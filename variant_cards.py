#!/usr/bin/env python3
import argparse
import os
import re
import sys
import time
from pathlib import Path
from typing import Optional, Tuple, List

import pandas as pd
import requests

# --- Constants / heuristics ---
POSSIBLE_GENE_COLS = ["GeneSymbol", "GeneSymbol;HGNC_ID", "GENE", "gene", "SYMBOL", "Symbol"]
POSSIBLE_DISEASE_COLS = [
    "PhenotypeList", "Phenotype", "Condition(s)", "ConditionList", "DiseaseName",
    "Condition", "Disease", "PHENOTYPE"
]
POSSIBLE_PMID_COLS = ["PubMedIDs", "PUBMED_IDS", "PMIDs", "PMID", "pubmed_id"]
POSSIBLE_NAME_COLS = ["Name", "VariantName", "VARIANT_NAME"]
POSSIBLE_HGVSC_COLS = ["HGVSc", "HGVS_cDNA", "HGVS_c", "hgvs_c"]
POSSIBLE_HGVSP_COLS = ["HGVSp", "Protein_change", "HGVS_p", "hgvs_p", "ProteinChange"]

NCBI_TOOL = "gene-browser-demo"
NCBI_EMAIL = os.environ.get("NCBI_EMAIL", "")  # set this if you have one (recommended by NCBI)


def pick_first_existing(df: pd.DataFrame, candidates: List[str]) -> Optional[str]:
    for c in candidates:
        if c in df.columns:
            return c
    return None


def normalize_gene_value(x: str) -> str:
    return str(x).strip().upper()


def extract_mutation_codes(row) -> Tuple[Optional[str], Optional[str]]:
    """Return (cDNA_code, protein_code) if present, using multiple column fallbacks and regex parsing."""
    # 1) direct columns
    for c in POSSIBLE_HGVSC_COLS:
        if c in row and pd.notna(row[c]) and str(row[c]).strip():
            cdna = str(row[c]).strip()
            break
    else:
        cdna = None

    for c in POSSIBLE_HGVSP_COLS:
        if c in row and pd.notna(row[c]) and str(row[c]).strip():
            prot = str(row[c]).strip()
            break
    else:
        prot = None

    # 2) parse from Name-like column if missing
    if (cdna is None or prot is None):
        name_col = None
        for c in POSSIBLE_NAME_COLS:
            if c in row:
                name_col = c
                break
        name_val = str(row[name_col]).strip() if (name_col and pd.notna(row[name_col])) else ""

        if cdna is None:
            m = re.search(r"(c\.[^ \),;]+)", name_val)
            if m:
                cdna = m.group(1)

        if prot is None:
            # typical formats: "(p.Glu23Valfs)", "p.Glu23Valfs", "(p.Trp24*)"
            m = re.search(r"\(?(p\.[A-Za-z][^)\s]+)\)?", name_val)
            if m:
                prot = m.group(1)

    return cdna, prot


def extract_first_disease(row) -> Optional[str]:
    for c in POSSIBLE_DISEASE_COLS:
        if c in row and pd.notna(row[c]):
            # ClinVar often has semicolon- or comma-separated lists
            txt = str(row[c]).strip()
            if not txt:
                continue
            # Split on semicolon or comma, pick the first non-empty chunk
            parts = re.split(r"[;|,]", txt)
            for p in parts:
                p2 = p.strip()
                if p2 and p2.lower() not in {"not provided", "not specified"}:
                    return p2
    return None


def parse_first_pmid(row) -> Optional[str]:
    for c in POSSIBLE_PMID_COLS:
        if c in row and pd.notna(row[c]):
            raw = str(row[c]).strip()
            if not raw:
                continue
            # split on ; or , or whitespace
            ids = re.split(r"[;,|\s]+", raw)
            for pm in ids:
                pm2 = pm.strip()
                if pm2.isdigit():
                    return pm2
    return None


def fetch_pubmed_one_sentence(pmid: str, timeout: float = 10.0) -> Optional[str]:
    """
    Fetch the PubMed abstract text (first sentence) for a PMID via NCBI E-utilities.
    Returns None if any error occurs.
    """
    if not pmid:
        return None

    # Use efetch with rettype=abstract, retmode=text for simplicity.
    # Respect NCBI etiquette (1–3 req/sec), include tool/email if available.
    params = {
        "db": "pubmed",
        "id": pmid,
        "rettype": "abstract",
        "retmode": "text",
    }
    if NCBI_TOOL:
        params["tool"] = NCBI_TOOL
    if NCBI_EMAIL:
        params["email"] = NCBI_EMAIL

    try:
        resp = requests.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
                            params=params, timeout=timeout)
        resp.raise_for_status()
        text = resp.text.strip()
        if not text:
            return None
        # Grab the first non-empty sentence-ish line.
        # Abstracts come with headers; pick the first line that looks like content.
        lines = [ln.strip() for ln in text.splitlines() if ln.strip()]
        if not lines:
            return None
        # Join the first paragraph then split into sentence(s)
        para = lines[0]
        # Sometimes the first non-empty line is a title; try the next non-empty line if too short
        if len(para.split()) < 6 and len(lines) > 1:
            para = lines[1]
        # naive sentence split
        sent = re.split(r"(?<=[.!?])\s+", para)
        first = sent[0].strip()
        # keep it brief
        return first if first else None
    except Exception:
        return None


def build_one_sentence_fallback(gene: str, cdna: Optional[str], prot: Optional[str],
                                disease: Optional[str], significance: Optional[str]) -> str:
    parts = []
    if gene:
        parts.append(gene)
    if cdna or prot:
        mut = cdna or prot
        parts.append(f"variant {mut}")
    if disease:
        parts.append(f"is associated with {disease}")
    if significance:
        parts.append(f"({significance})")
    txt = " ".join(parts).strip()
    if not txt:
        txt = "Variant associated information not available."
    if not txt.endswith("."):
        txt += "."
    return txt


def main():
    ap = argparse.ArgumentParser(description="Print 2–3 variant cards from a ClinVar file for a given gene.")
    ap.add_argument("--gene", required=True, help="Gene symbol, e.g. BRCA1")
    ap.add_argument("--clinvar", required=True, help="Path to ClinVar TSV/CSV (BRCA1 subset to start)")
    ap.add_argument("-n", "--num", type=int, default=3, help="Number of variant cards to print (default: 3)")
    args = ap.parse_args()

    path = Path(args.clinvar)
    if not path.exists():
        print(f"ERROR: File not found: {path}", file=sys.stderr)
        sys.exit(1)

    # Load with pandas. Try TSV then CSV automatically.
    try:
        if path.suffix.lower() in {".tsv", ".txt"}:
            df = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
        else:
            df = pd.read_csv(path, dtype=str, low_memory=False)
    except Exception as e:
        print(f"ERROR: Could not read {path}: {e}", file=sys.stderr)
        sys.exit(1)

    # Normalize columns
    df.columns = [c.strip() for c in df.columns]

    # Identify the gene column
    gene_col = pick_first_existing(df, POSSIBLE_GENE_COLS)
    if gene_col is None:
        # Try to derive gene from Name column later; for now, warn
        print("WARN: No obvious gene column found; will try Name parsing only.", file=sys.stderr)

    # Filter by gene
    gene = normalize_gene_value(args.gene)
    if gene_col:
        df["_gene_norm"] = df[gene_col].fillna("").map(normalize_gene_value)
        hits = df[df["_gene_norm"] == gene].copy()
    else:
        hits = df.copy()

    # If we still have no rows, try parsing gene symbol out of Name column (e.g., NM_007294.3(BRCA1):c.68_69delAG ...)
    if hits.empty:
        name_col = pick_first_existing(df, POSSIBLE_NAME_COLS)
        if name_col:
            mask = df[name_col].fillna("").str.contains(rf"\({re.escape(gene)}\)", na=False, regex=True)
            hits = df[mask].copy()

    if hits.empty:
        print(f"No variants found for gene '{gene}' in {path}.")
        sys.exit(0)

    # Pick some columns we might want for display
    # Clinical significance is often useful if present
    significance_col = None
    for c in ["ClinicalSignificance", "clinical_significance", "CLNSIG", "Significance", "CLIN_SIG"]:
        if c in df.columns:
            significance_col = c
            break

    # Limit to N
    hits = hits.head(max(1, args.num))

    # Print cards
    cards = []
    for _, row in hits.iterrows():
        cdna, prot = extract_mutation_codes(row)
        disease = extract_first_disease(row)
        pmid = parse_first_pmid(row)
        significance = row.get(significance_col, None)

        summary = None
        # Try PubMed if we have a PMID
        if pmid:
            summary = fetch_pubmed_one_sentence(pmid)
            # Be polite to NCBI
            time.sleep(0.35)

        if not summary or len(summary.split()) < 5:
            summary = build_one_sentence_fallback(gene, cdna, prot, disease, significance)

        # Pick a variant label
        mut_label = cdna or prot or row.get(pick_first_existing(df, POSSIBLE_NAME_COLS), "Unspecified variant")

        cards.append({
            "Gene": gene,
            "Mutation": mut_label,
            "Disease/Phenotype": disease or "Not specified",
            "Clinical significance": significance if pd.notna(significance) else "Not specified",
            "PMID": pmid or "—",
            "Summary": summary
        })

    # Pretty print
    for i, c in enumerate(cards, 1):
        print("=" * 60)
        print(f"Variant {i}: {c['Mutation']}")
        print("- Gene:                 ", c["Gene"])
        print("- Disease/Phenotype:    ", c["Disease/Phenotype"])
        print("- Clinical significance:", c["Clinical significance"])
        print("- PMID:                 ", c["PMID"])
        print("- Summary:              ", c["Summary"])
    print("=" * 60)


if __name__ == "__main__":
    main()

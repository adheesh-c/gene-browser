from pathlib import Path
import pandas as pd
import streamlit as st
try:
    from rapidfuzz import process, fuzz
    def fuzzy_gene_candidates(q, choices, limit=5):
        matches = process.extract(q, choices, scorer=fuzz.WRatio, limit=limit)
        return [m[0] for m in matches if m[1] > 60]
except Exception:
    # built-in fallback (slower)
    from difflib import get_close_matches
    def fuzzy_gene_candidates(q, choices, limit=5):
        return get_close_matches(q, choices, n=limit, cutoff=0.6)

# --- NEW IMPORTS ---
import time, re, os
from typing import Optional, Tuple, List
import requests

# Identify likely ClinVar column names
POSSIBLE_PMID_COLS = ["PubMedIDs", "PUBMED_IDS", "PMIDs", "PMID", "pubmed_id"]
POSSIBLE_NAME_COLS = ["Name", "VariantName", "VARIANT_NAME"]  # optional for mutation label

# NCBI etiquette
NCBI_TOOL = "gene-browser"
NCBI_EMAIL = st.secrets.get("NCBI_EMAIL", "")

# --- Column name candidates seen in ClinVar dumps ---
POSSIBLE_GENE_COLS   = ["GeneSymbol","GeneSymbol;HGNC_ID","GENE","gene","SYMBOL","Symbol"]
POSSIBLE_DISEASE_COLS= ["PhenotypeList","Phenotype","Condition(s)","ConditionList","DiseaseName","Condition","Disease","PHENOTYPE"]
POSSIBLE_PMID_COLS   = ["PubMedIDs","PUBMED_IDS","PMIDs","PMID","pubmed_id"]
POSSIBLE_NAME_COLS   = ["Name","VariantName","VARIANT_NAME"]
POSSIBLE_HGVSC_COLS  = ["HGVSc","HGVS_cDNA","HGVS_c","hgvs_c"]
POSSIBLE_HGVSP_COLS  = ["HGVSp","Protein_change","HGVS_p","hgvs_p","ProteinChange"]

def _pick_first(df: pd.DataFrame, candidates: List[str]) -> Optional[str]:
    for c in candidates:
        if c in df.columns: return c
    return None

def _normalize_gene(x: str) -> str:
    return str(x).strip().upper()

def _extract_mutation_codes(row) -> Tuple[Optional[str], Optional[str]]:
    cdna = None
    prot = None
    for c in POSSIBLE_HGVSC_COLS:
        if c in row and pd.notna(row[c]) and str(row[c]).strip():
            cdna = str(row[c]).strip(); break
    for c in POSSIBLE_HGVSP_COLS:
        if c in row and pd.notna(row[c]) and str(row[c]).strip():
            prot = str(row[c]).strip(); break
    if cdna is None or prot is None:
        name_col = next((c for c in POSSIBLE_NAME_COLS if c in row), None)
        name_val = str(row[name_col]).strip() if (name_col and pd.notna(row[name_col])) else ""
        if cdna is None:
            m = re.search(r"(c\.[^ \),;]+)", name_val)
            if m: cdna = m.group(1)
        if prot is None:
            m = re.search(r"\(?(p\.[A-Za-z][^)\s]+)\)?", name_val)
            if m: prot = m.group(1)
    return cdna, prot

def _extract_first_disease(row) -> Optional[str]:
    for c in POSSIBLE_DISEASE_COLS:
        if c in row and pd.notna(row[c]):
            txt = str(row[c]).strip()
            if not txt: continue
            parts = re.split(r"[;|,]", txt)
            for p in parts:
                p2 = p.strip()
                if p2 and p2.lower() not in {"not provided","not specified"}:
                    return p2
    return None

def _first_pmid(row) -> Optional[str]:
    for c in POSSIBLE_PMID_COLS:
        if c in row and pd.notna(row[c]):
            raw = str(row[c]).strip()
            ids = re.split(r"[;,|\s]+", raw)
            for pm in ids:
                pm2 = pm.strip()
                if pm2.isdigit(): return pm2
    return None

@st.cache_data(ttl=60*60, show_spinner=False)
def _pubmed_one_sentence(pmid: str) -> Optional[str]:
    if not pmid: return None
    params = {"db":"pubmed","id":pmid,"rettype":"abstract","retmode":"text"}
    if NCBI_TOOL: params["tool"] = NCBI_TOOL
    if NCBI_EMAIL: params["email"] = NCBI_EMAIL
    try:
        r = requests.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
                         params=params, timeout=10)
        r.raise_for_status()
        text = r.text.strip()
        if not text: return None
        lines = [ln.strip() for ln in text.splitlines() if ln.strip()]
        if not lines: return None
        para = lines[0] if len(lines[0].split())>=6 else (lines[1] if len(lines)>1 else lines[0])
        sent = re.split(r"(?<=[.!?])\s+", para)[0].strip()
        return sent or None
    except Exception:
        return None

def _fallback_sentence(gene: str, cdna: Optional[str], prot: Optional[str],
                       disease: Optional[str], significance: Optional[str]) -> str:
    bits = []
    if gene: bits.append(gene)
    if cdna or prot: bits.append(f"variant {(cdna or prot)}")
    if disease: bits.append(f"is associated with {disease}")
    if significance and str(significance).strip(): bits.append(f"({significance})")
    txt = " ".join(bits).strip()
    if not txt: txt = "Variant associated information not available."
    if not txt.endswith("."): txt += "."
    return txt

def build_variant_cards(df: pd.DataFrame, gene: str, n: int = 3) -> list[dict]:
    """Return up to n cards with mutation label, disease, significance, PMID, and one-sentence summary."""
    gene_norm = str(gene).strip().upper()
    df = df.copy()
    df.columns = [c.strip() for c in df.columns]

    # Columns we might use for display
    sig_col  = next((c for c in ["ClinicalSignificance","clinical_significance","CLNSIG","Significance","CLIN_SIG"] if c in df.columns), None)
    name_col = next((c for c in ["HGVSc","HGVS_cDNA","HGVS_c","hgvs_c","HGVSp","Protein_change","HGVS_p","hgvs_p","ProteinChange"] if c in df.columns), None)
    if not name_col:
        name_col = next((c for c in POSSIBLE_NAME_COLS if c in df.columns), None)

    # Preferred gene column(s)
    gene_col = next((c for c in ["GeneSymbol","GENE","gene","SYMBOL","Symbol","GeneSymbol;HGNC_ID"] if c in df.columns), None)

    # Filter by gene
    if gene_col:
        df["_gene_norm"] = df[gene_col].fillna("").map(lambda x: str(x).strip().upper())
        hits = df[df["_gene_norm"] == gene_norm]
    else:
        hits = df

    # Fallback: match (GENE) in "Name" field if no direct gene column matched
    if hits.empty and "Name" in df.columns:
        mask = df["Name"].fillna("").str.contains(rf"\({re.escape(gene_norm)}\)", regex=True, na=False)
        hits = df[mask]

    # Take top n rows (you can sort by significance/date if you prefer)
    hits = hits.head(max(1, n))

    cards = []
    for _, row in hits.iterrows():
        # Mutation label preference: cDNA, then protein, then Name as fallback
        mut_label = None
        for c in ["HGVSc","HGVS_cDNA","HGVS_c","hgvs_c","HGVSp","Protein_change","HGVS_p","hgvs_p","ProteinChange"]:
            if c in row and pd.notna(row[c]) and str(row[c]).strip():
                mut_label = str(row[c]).strip()
                break
        if not mut_label and name_col and pd.notna(row.get(name_col)):
            mut_label = str(row[name_col]).strip()

        # Disease/phenotype (first from a list-like column set)
        disease = None
        for c in ["PhenotypeList","Phenotype","Condition(s)","ConditionList","DiseaseName","Condition","Disease","PHENOTYPE"]:
            if c in row and pd.notna(row[c]):
                raw = str(row[c]).strip()
                if raw:
                    first = re.split(r"[;,\|]", raw)[0].strip()
                    if first and first.lower() not in {"not provided","not specified"}:
                        disease = first
                        break

        significance = str(row.get(sig_col, "") or "").strip() if sig_col else ""

        # --- NEW: PubMed one-sentence summary ---
        pmid = _parse_first_pmid(row)
        summary = _pubmed_first_sentence(pmid) if pmid else None
        # Be polite to NCBI if we actually hit the API (cache prevents most calls)
        if pmid and summary is None:
            time.sleep(0.35)

        if not summary or len(summary.split()) < 5:
            summary = _fallback_sentence(gene_norm, mut_label, disease, significance)

        cards.append({
            "Mutation": mut_label or "Unspecified variant",
            "Disease/Phenotype": disease or "Not specified",
            "Clinical significance": significance or "Not specified",
            "PMID": pmid or "—",
            "Summary": summary,
        })

    return cards


DATA_PATH = Path(__file__).parent / "clinvar_sample.csv"



# st.write("CWD:", Path.cwd())
# st.write("Expected data path:", DATA_PATH)
# st.write("File exists?", DATA_PATH.exists())

# ---------- Page config ----------
st.set_page_config(page_title="Gene → Mutations", page_icon="🧬", layout="wide")

# ---------- Caching ----------
@st.cache_data(show_spinner=False)
def load_variants(csv_path: Path) -> pd.DataFrame:
    if not csv_path.exists():
        st.error(f"Data file not found: {csv_path}")
        st.stop()
        # Fallback empty frame with expected columns
        cols = ["gene","variant_id","protein_change","cdna_change",
                "clinical_significance","condition","source"]
        return pd.DataFrame(columns=cols)
    if csv_path.suffix.lower() in {".tsv", ".txt"}:
        return pd.read_csv(csv_path, sep="\t", dtype=str, low_memory=False)
    return pd.read_csv(csv_path, dtype=str, low_memory=False)
    df = pd.read_csv(csv_path)
    # Normalize column names for safety
    df.columns = [c.strip().lower() for c in df.columns]
    # Optional: strip whitespace, uppercase gene symbols
    if "gene" in df.columns:
        df["gene"] = df["gene"].astype(str).str.strip().str.upper()
    return df

variants_df = load_variants(DATA_PATH)

def canonicalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Create stable, lower-case columns your UI expects, regardless of ClinVar header variants."""
    df = df.copy()
    df.columns = [c.strip() for c in df.columns]

    # gene → create 'gene'
    gene_col = _pick_first(df, ["gene", "GENE", "GeneSymbol", "SYMBOL", "Symbol", "GeneSymbol;HGNC_ID"])
    if gene_col:
        df["gene"] = df[gene_col].fillna("").map(_normalize_gene)

    # clinical_significance → create 'clinical_significance'
    sig_col = _pick_first(df, ["clinical_significance","ClinicalSignificance","CLNSIG","Significance","CLIN_SIG"])
    if sig_col and "clinical_significance" not in df.columns:
        df["clinical_significance"] = df[sig_col]

    # condition/phenotype → create 'condition'
    cond_col = _pick_first(df, POSSIBLE_DISEASE_COLS)
    if cond_col and "condition" not in df.columns:
        df["condition"] = df[cond_col]

    # variant_id (optional)
    vid_col = _pick_first(df, ["VariantID","variant_id","VariationID","VCV","VCV_ID"])
    if vid_col and "variant_id" not in df.columns:
        df["variant_id"] = df[vid_col]

    # protein change → create 'protein_change'
    prot_col = _pick_first(df, ["protein_change","HGVSp","Protein_change","HGVS_p","hgvs_p","ProteinChange"])
    if prot_col and "protein_change" not in df.columns:
        df["protein_change"] = df[prot_col]

    # cDNA change → create 'cdna_change'
    cdna_col = _pick_first(df, ["cdna_change","HGVSc","HGVS_cDNA","HGVS_c","hgvs_c"])
    if cdna_col and "cdna_change" not in df.columns:
        df["cdna_change"] = df[cdna_col]

    # source (optional)
    if "source" not in df.columns:
        df["source"] = "ClinVar"

    return df

variants_df = canonicalize_columns(variants_df)


def _parse_first_pmid(row) -> str | None:
    """Extract the first numeric PMID from any of the usual ClinVar PMID columns."""
    for col in POSSIBLE_PMID_COLS:
        if col in row and pd.notna(row[col]):
            raw = str(row[col]).strip()
            if not raw:
                continue
            # Split on semicolons, commas, pipes, whitespace
            for pm in re.split(r"[;,\|\s]+", raw):
                pm = pm.strip()
                if pm.isdigit():
                    return pm
    return None

@st.cache_data(ttl=60*60, show_spinner=False)
def _pubmed_first_sentence(pmid: str) -> str | None:
    """Fetch PubMed abstract (text mode) and return a single first sentence."""
    if not pmid:
        return None
    params = {"db": "pubmed", "id": pmid, "rettype": "abstract", "retmode": "text"}
    if NCBI_TOOL:
        params["tool"] = NCBI_TOOL
    if NCBI_EMAIL:
        params["email"] = NCBI_EMAIL

    try:
        r = requests.get(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
            params=params,
            timeout=10,
        )
        r.raise_for_status()
        text = r.text.strip()
        if not text:
            return None

        # Find a content-y line (skip headers). If first line is very short, try next.
        lines = [ln.strip() for ln in text.splitlines() if ln.strip()]
        if not lines:
            return None
        para = lines[0] if len(lines[0].split()) >= 6 else (lines[1] if len(lines) > 1 else lines[0])

        # Naive sentence split; good enough for a one-liner
        sent = re.split(r"(?<=[.!?])\s+", para)[0].strip()
        return sent or None

    except Exception:
        # Keep the UI resilient if PubMed fails
        return None

def _fallback_sentence(gene: str, mutation_label: str | None, disease: str | None, significance: str | None) -> str:
    bits = []
    if gene: bits.append(gene)
    if mutation_label: bits.append(f"variant {mutation_label}")
    if disease: bits.append(f"is associated with {disease}")
    if significance and str(significance).strip(): bits.append(f"({significance})")
    txt = " ".join(bits).strip()
    if not txt:
        txt = "Variant associated information not available."
    if not txt.endswith("."):
        txt += "."
    return txt


# ---------- Sidebar (filters & info) ----------
with st.sidebar:
    st.header("Filters")
    gene_options = sorted(variants_df["gene"].dropna().unique().tolist()) if "gene" in variants_df.columns else []
    sig_options = sorted(variants_df["clinical_significance"].dropna().unique().tolist()) if "clinical_significance" in variants_df.columns else []
    cond_options = sorted(variants_df["condition"].dropna().unique().tolist()) if "condition" in variants_df.columns else []

    selected_sigs = st.multiselect("Clinical significance", options=sig_options, default=[])
    selected_conditions = st.multiselect("Condition", options=cond_options, default=[])

    st.markdown("---")
    st.caption("Data source: ClinVar subset. PubMed summaries via NCBI E-utilities.")


# ---------- Title & search ----------
st.title("🧬 Gene → Related Mutations")
st.write("Type a gene symbol (e.g., **BRCA1**) to see related variants. Use the sidebar to filter.")

col1, col2 = st.columns([2,1])
with col1:
    query = st.text_input("Gene symbol", value="", placeholder="e.g., BRCA1, TP53, BRCA2").strip().upper()
with col2:
    topk = st.number_input("Max results", min_value=5, max_value=5000, value=200, step=5)

# ---------- Autocomplete / fuzzy help ----------
# def fuzzy_gene_candidates(user_text: str, choices: list[str], limit: int = 5):
#     if not user_text or not choices:
#         return []
#     # RapidFuzz returns (match, score, index)
#     matches = process.extract(user_text, choices, scorer=fuzz.WRatio, limit=limit)
#     return [m[0] for m in matches if m[1] > 60]  # simple threshold

if query and (query not in variants_df["gene"].unique()):
    suggestions = fuzzy_gene_candidates(query, gene_options, limit=5)
    if suggestions:
        st.info("Did you mean:")
        cols = st.columns(len(suggestions))
        for i, s in enumerate(suggestions):
            if cols[i].button(s):
                st.session_state["__query_override"] = s

# apply override if user clicked a suggestion
query = st.session_state.get("__query_override", query)

# ---------- Filtering logic ----------
def filter_variants(df: pd.DataFrame, gene_text: str) -> pd.DataFrame:
    if not gene_text:
        return df.iloc[0:0]  # show nothing until user types
    # Exact gene match first
    exact = df[df["gene"] == gene_text]
    if len(exact) > 0:
        return exact
    # Fall back to fuzzy gene filtering (startswith or contains)
    starts = df[df["gene"].str.startswith(gene_text, na=False)]
    if len(starts) > 0:
        return starts
    contains = df[df["gene"].str.contains(gene_text, na=False)]
    return contains

results = filter_variants(variants_df, query)

# Apply sidebar filters
if selected_sigs:
    results = results[results["clinical_significance"].isin(selected_sigs)]
if selected_conditions:
    results = results[results["condition"].isin(selected_conditions)]

# Limit and display
results = results.head(int(topk))

# ---------- Summary chips ----------
left, right = st.columns(2)
with left:
    st.metric("Matches", len(results))
with right:
    st.metric("Unique conditions", results["condition"].nunique() if not results.empty else 0)

# ---------- Table ----------
if results.empty and query:
    st.warning("No variants found for that input and filters. Try removing filters or check the gene symbol.")
elif query:
    # Order columns nicely if present
    preferred_cols = ["gene", "variant_id", "protein_change", "cdna_change",
                      "clinical_significance", "condition", "source"]
    ordered_cols = [c for c in preferred_cols if c in results.columns] + \
                   [c for c in results.columns if c not in preferred_cols]
    st.dataframe(results[ordered_cols], use_container_width=True, hide_index=True)

    # Download button
    csv_bytes = results[ordered_cols].to_csv(index=False).encode("utf-8")
    st.download_button(
        label="Download results (CSV)",
        data=csv_bytes,
        file_name=f"{query}_variants.csv",
        mime="text/csv",
    )


if query and not results.empty:
    with st.expander("🧬 Variant cards with PubMed summary", expanded=True):
        left, right = st.columns([1, 5])
        with left:
            n_cards = st.number_input("How many", min_value=1, max_value=5, value=3, step=1)
        with right:
            st.caption("Each card includes a one-sentence summary extracted from PubMed when available.")
        if st.button("Generate cards"):
            cards = build_variant_cards(variants_df, query, n=int(n_cards))
            for i, c in enumerate(cards, 1):
                st.markdown("---")
                st.markdown(f"**Variant {i}: {c['Mutation']}**")
                st.markdown(f"- **Disease/Phenotype:** {c['Disease/Phenotype']}")
                st.markdown(f"- **Clinical significance:** {c['Clinical significance']}")
                if c["PMID"] != "—":
                    st.markdown(f"- **PMID:** [{c['PMID']}](https://pubmed.ncbi.nlm.nih.gov/{c['PMID']}/)")
                st.markdown(f"> {c['Summary']}")

# ---------- Footer helpers ----------
# with st.expander("How to plug in real ClinVar / NCBI"):
#     st.markdown("""
# - Replace `data/clinvar_sample.csv` with an export from ClinVar or your ETL.
# - Keep column names similar (gene, variant_id, protein_change, cdna_change, clinical_significance, condition, source).
# - For performance on large data:
#   - Use **parquet** instead of CSV and load with `pd.read_parquet`.
#   - Add `@st.cache_data` to any expensive loaders or API calls.
#   - Consider pre-indexing by gene into per-gene Parquet files.
# """)

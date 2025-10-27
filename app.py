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
# --- Column name candidates seen in ClinVar dumps + your CSV ---
POSSIBLE_GENE_COLS    = ["GeneSymbol","GeneSymbol;HGNC_ID","GENE","gene","SYMBOL","Symbol"]
POSSIBLE_DISEASE_COLS = ["PhenotypeList","Phenotype","Condition(s)","ConditionList","DiseaseName",
                         "Condition","Disease","PHENOTYPE","condition","disease","phenotype"]
POSSIBLE_PMID_COLS    = ["PubMedIDs","PUBMED_IDS","PMIDs","PMID","pubmed_id","pmid"]
POSSIBLE_NAME_COLS    = ["Name","VariantName","VARIANT_NAME"]
# Include your CSV's field names here ↓↓↓
POSSIBLE_HGVSC_COLS   = ["HGVSc","HGVS_cDNA","HGVS_c","hgvs_c","cdna_change"]
POSSIBLE_HGVSP_COLS   = ["HGVSp","Protein_change","HGVS_p","hgvs_p","ProteinChange","protein_change"]


# ---------- PubMed helpers (robust) ----------
@st.cache_data(ttl=60*60, show_spinner=False)
def _pubmed_one_sentence(pmid: str) -> str | None:
    """
    Return first sentence from PubMed abstract using XML (more reliable than text mode).
    """
    if not pmid:
        return None
    params = {"db": "pubmed", "id": pmid, "retmode": "xml"}
    if NCBI_TOOL: params["tool"] = NCBI_TOOL
    if NCBI_EMAIL: params["email"] = NCBI_EMAIL

    try:
        r = requests.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
                         params=params, timeout=10)
        r.raise_for_status()
        from xml.etree import ElementTree as ET
        root = ET.fromstring(r.text)
        texts = []
        for node in root.findall(".//Abstract/AbstractText"):
            part = (node.text or "").strip()
            if part:
                texts.append(part)
        abstract = " ".join(texts).strip()
        if not abstract:
            return None
        import re as _re
        sent = _re.split(r"(?<=[.!?])\s+", abstract)[0].strip()
        return sent or None
    except Exception:
        return None

@st.cache_data(ttl=60*60, show_spinner=False)
def _guess_pmid(gene: str | None, mutation: str | None, disease: str | None) -> str | None:
    """
    Try multiple PubMed ESearch queries to find a relevant PMID.
    Priority: gene + mutation, then gene + disease, then gene only.
    """
    g = (gene or "").strip()
    m = (mutation or "").strip()
    d = (disease or "").strip()

    queries = []
    if g and m:
        queries.append(f'{g}[Title/Abstract] AND ("{m}"[All Fields])')
    if g and d:
        queries.append(f'{g}[Title/Abstract] AND ("{d}"[Title/Abstract])')
    if g:
        queries.append(f'{g}[Title/Abstract]')

    for q in queries:
        params = {"db": "pubmed", "retmode": "json", "retmax": "1", "sort": "relevance", "term": q}
        if NCBI_TOOL: params["tool"] = NCBI_TOOL
        if NCBI_EMAIL: params["email"] = NCBI_EMAIL
        try:
            r = requests.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
                             params=params, timeout=10)
            r.raise_for_status()
            ids = r.json().get("esearchresult", {}).get("idlist", [])
            if ids:
                return ids[0]
        except Exception:
            pass
    return None


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
            if not txt:
                continue
            parts = re.split(r"[;|,]", txt)
            for p in parts:
                p2 = p.strip()
                if p2 and p2.lower() not in {"not provided","not specified","na","n/a"}:
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
        # Mutation label preference: cDNA, then protein, then Name as fallback (includes your CSV fields)
        mut_label = None
        for c in ["cdna_change","HGVSc","HGVS_cDNA","HGVS_c","hgvs_c",
                  "protein_change","HGVSp","Protein_change","HGVS_p","hgvs_p","ProteinChange"]:
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
        
        if not pmid:
            pmid = _guess_pmid(gene_norm, mut_label, disease)

        summary = _pubmed_one_sentence(pmid) if pmid else None

        # Be polite to NCBI if we actually hit the API and got nothing (will be cached on future calls)
        if pmid and summary is None:
            time.sleep(0.35)

        if not summary or len(summary.split()) < 5:
            summary = _fallback_sentence(gene_norm, mut_label, disease, significance)




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

st.markdown("""
<style>
.hero {
  padding: 2rem; border-radius: 16px;
  background: linear-gradient(180deg, #f5f6ff 0%, #ffffff 100%);
  border: 1px solid #e9ecff;
}
.hero h1 { margin: 0 0 0.5rem 0; }
.badge { display:inline-block; padding:0.25rem 0.6rem; border-radius:999px; background:#eef2ff; margin-right:0.4rem; }
.kidcard { border:1px solid #e6e8ff; border-radius:14px; padding:1rem; background:#fafbff; }
</style>
""", unsafe_allow_html=True)


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

# --- PMID parsing from dataset ---
def _parse_first_pmid(row) -> str | None:
    for col in POSSIBLE_PMID_COLS:
        if col in row and pd.notna(row[col]):
            raw = str(row[col]).strip()
            if not raw:
                continue
            for pm in re.split(r"[;,\|\s]+", raw):
                pm = pm.strip()
                if pm.isdigit():
                    return pm
    return None


# --- NEW: Attach a PMID column to your dataframe (manual or guessed) ---
@st.cache_data(ttl=60*60, show_spinner=False)
def attach_pmids(df: pd.DataFrame) -> pd.DataFrame:
    """
    Returns a copy of df with a 'PMID' column added if missing.
    Uses any existing PMID in the row; otherwise tries to guess via PubMed.
    """
    if "PMID" in df.columns:  # already there
        return df.copy()

    out_rows = []
    for _, r in df.iterrows():
        pmid = _parse_first_pmid(r)
        if not pmid:
            pmid = _guess_pmid(
                gene=str(r.get("gene","")),
                mutation=str(r.get("cdna_change") or r.get("protein_change") or ""),
                disease=str(r.get("condition",""))
            )
        rr = r.copy()
        rr["PMID"] = pmid or ""
        out_rows.append(rr)
    return pd.DataFrame(out_rows)


variants_df = load_variants(DATA_PATH)
variants_df = attach_pmids(variants_df)

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

# Teen-friendly landing when no query yet
if not query:
    st.markdown('<div class="hero">', unsafe_allow_html=True)
    st.markdown("## 🧬 Welcome to the Gene Variant Explorer")
    st.write("Type a gene (like **BRCA1**) to see real variants, the conditions they’re linked to, and a one-sentence research summary from PubMed.")
    st.write("Use these quick-picks to try it:")
    cols = st.columns(5)
    picks = ["BRCA1", "BRCA2", "TP53", "APC", "MLH1"]
    for i, g in enumerate(picks):
        if cols[i].button(g):
            st.session_state["__query_override"] = g
            st.rerun()
    st.markdown('</div>', unsafe_allow_html=True)

    with st.expander("What’s a gene? What’s a variant? (Simple explanations)"):
        st.markdown("""
- **Gene:** an instruction in your DNA that helps your cells work (like a recipe).
- **Variant (mutation):** a small change in that instruction. Some are harmless, some matter for health.
- **Why look?** Understanding common variants can guide better health conversations with doctors.
*This app is for learning only — not medical advice.*
        """)

# Apply quick-pick override
query = st.session_state.get("__query_override", query)


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
    st.dataframe(results[ordered_cols], width="stretch", hide_index=True)

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
                else:
                    st.caption("No PMID in dataset; attempted PubMed search.")
                st.markdown(f"> {c['Summary']}")

            with st.expander("Diagnostics (why you might not see a one-liner)"):
                st.write("Showing up to the first 5 candidate rows with their derived PubMed info.")
                import itertools
                diag_rows = []
                # build again but keep pmid + raw query ingredients
                sample_hits = results.head(min(len(results), 5))
                for _, row in sample_hits.iterrows():
                    gene_norm = (query or "").strip().upper()
                    mut = str(row.get("cdna_change") or row.get("protein_change") or row.get("Name") or "")
                    dis = str(row.get("condition") or row.get("Phenotype") or "")
                    pmid = row.get("PMID", "") or _parse_first_pmid(row) or _guess_pmid(gene_norm, mut, dis) or ""
                    summary = _pubmed_one_sentence(pmid) if pmid else ""
                    diag_rows.append({
                        "gene": gene_norm,
                        "mutation": mut,
                        "disease": dis,
                        "pmid": pmid,
                        "summary_len": len(summary.split()) if summary else 0
                    })
                st.dataframe(pd.DataFrame(diag_rows), width="stretch", hide_index=True)


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

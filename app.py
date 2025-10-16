import os
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

DATA_PATH = Path(__file__).parent / "data" / "clinvar_sample.csv"

# ---------- Page config ----------
st.set_page_config(page_title="Gene → Mutations", page_icon="🧬", layout="wide")

# ---------- Caching ----------
@st.cache_data(show_spinner=False)
def load_variants(csv_path: Path) -> pd.DataFrame:
    if not csv_path.exists():
        # Fallback empty frame with expected columns
        cols = ["gene","variant_id","protein_change","cdna_change",
                "clinical_significance","condition","source"]
        return pd.DataFrame(columns=cols)
    df = pd.read_csv(csv_path)
    # Normalize column names for safety
    df.columns = [c.strip().lower() for c in df.columns]
    # Optional: strip whitespace, uppercase gene symbols
    if "gene" in df.columns:
        df["gene"] = df["gene"].astype(str).str.strip().str.upper()
    return df

variants_df = load_variants(DATA_PATH)

# ---------- Sidebar (filters & info) ----------
with st.sidebar:
    st.header("Filters")
    # Precompute options
    gene_options = sorted(variants_df["gene"].dropna().unique().tolist())
    sig_options = sorted(variants_df["clinical_significance"].dropna().unique().tolist())
    cond_options = sorted(variants_df["condition"].dropna().unique().tolist())

    selected_sigs = st.multiselect("Clinical significance", options=sig_options, default=[])
    selected_conditions = st.multiselect("Condition", options=cond_options, default=[])

    st.markdown("---")
    st.caption("Data source: sample CSV (swap with ClinVar export/API later).")

# ---------- Title & search ----------
st.title("🧬 Gene → Related Mutations")
st.write("Type a gene symbol (e.g., **BRCA1**) to see related variants. Use the sidebar to filter.")

col1, col2 = st.columns([2,1])
with col1:
    query = st.text_input("Gene symbol", value="", placeholder="e.g., BRCA1, TP53, BRCA2").strip().upper()
with col2:
    topk = st.number_input("Max results", min_value=5, max_value=5000, value=200, step=5)

# ---------- Autocomplete / fuzzy help ----------
def fuzzy_gene_candidates(user_text: str, choices: list[str], limit: int = 5):
    if not user_text or not choices:
        return []
    # RapidFuzz returns (match, score, index)
    matches = process.extract(user_text, choices, scorer=fuzz.WRatio, limit=limit)
    return [m[0] for m in matches if m[1] > 60]  # simple threshold

if query and (query not in variants_df["gene"].unique()):
    suggestions = fuzzy_gene_candidates(query, gene_options, limit=5)
    if suggestions:
        st.info(f"Did you mean: {', '.join(suggestions)}")

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

# ---------- Footer helpers ----------
with st.expander("How to plug in real ClinVar / NCBI"):
    st.markdown("""
- Replace `data/clinvar_sample.csv` with an export from ClinVar or your ETL.
- Keep column names similar (gene, variant_id, protein_change, cdna_change, clinical_significance, condition, source).
- For performance on large data:
  - Use **parquet** instead of CSV and load with `pd.read_parquet`.
  - Add `@st.cache_data` to any expensive loaders or API calls.
  - Consider pre-indexing by gene into per-gene Parquet files.
""")

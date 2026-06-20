from pathlib import Path
import pandas as pd
import streamlit as st

import time, re, os
from typing import Optional, Tuple, List
import requests
from datetime import datetime

# ----------------------------
# Fuzzy matching (RapidFuzz preferred)
# ----------------------------
try:
    from rapidfuzz import process, fuzz
    def fuzzy_gene_candidates(q, choices, limit=5):
        matches = process.extract(q, choices, scorer=fuzz.WRatio, limit=limit)
        return [m[0] for m in matches if m[1] > 60]
except Exception:
    from difflib import get_close_matches
    def fuzzy_gene_candidates(q, choices, limit=5):
        return get_close_matches(q, choices, n=limit, cutoff=0.6)

# ----------------------------
# Query params helpers
# ----------------------------
def _get_param_list(key):
    val = st.query_params.get(key, "")
    if isinstance(val, list):
        return val
    if not val:
        return []
    return val.split("|")

def _set_params_keep_page(gene: str, sigs: list[str], conds: list[str]):
    """
    Update URL params used for sharing searches, WITHOUT overwriting the 'page' param.
    """
    page = st.query_params.get("page", "home")
    if isinstance(page, list):
        page = page[0] if page else "home"

    st.query_params["page"] = page
    st.query_params["gene"] = gene or ""
    st.query_params["sig"]  = "|".join(sigs) if sigs else ""
    st.query_params["cond"] = "|".join(conds) if conds else ""

# ----------------------------
# ClinVar column candidates
# ----------------------------
POSSIBLE_GENE_COLS    = ["GeneSymbol","GeneSymbol;HGNC_ID","GENE","gene","SYMBOL","Symbol"]
POSSIBLE_DISEASE_COLS = ["PhenotypeList","Phenotype","Condition(s)","ConditionList","DiseaseName",
                         "Condition","Disease","PHENOTYPE","condition","disease","phenotype"]
POSSIBLE_PMID_COLS    = ["PubMedIDs","PUBMED_IDS","PMIDs","PMID","pubmed_id","pmid","PMID"]
POSSIBLE_NAME_COLS    = ["Name","VariantName","VARIANT_NAME"]

POSSIBLE_HGVSC_COLS   = ["HGVSc","HGVS_cDNA","HGVS_c","hgvs_c","cdna_change"]
POSSIBLE_HGVSP_COLS   = ["HGVSp","Protein_change","HGVS_p","hgvs_p","ProteinChange","protein_change"]
POSSIBLE_CHROM_COLS   = ["Chromosome","chromosome","CHROM","chrom","Chr","chr"]
POSSIBLE_POS_COLS     = ["Start","start","POS","pos","Position","position"]

# NCBI etiquette
NCBI_TOOL = "gene-browser"
NCBI_EMAIL = st.secrets.get("NCBI_EMAIL", "")

# ----------------------------
# Gene summaries (hardcoded for common educational genes; NCBI API fallback for others)
# ----------------------------
GENE_SUMMARIES = {
    "BRCA1": "BRCA1 encodes a tumor suppressor that helps repair damaged DNA through homologous recombination. Pathogenic variants significantly increase lifetime risk of breast and ovarian cancer.",
    "BRCA2": "BRCA2 plays a central role in DNA repair alongside BRCA1. Variants are associated with elevated risk of breast, ovarian, prostate, and pancreatic cancers.",
    "TP53": "TP53 encodes p53, often called the 'guardian of the genome,' which halts cell division or triggers cell death when DNA is damaged. It is the most frequently mutated gene across all human cancers.",
    "CFTR": "CFTR encodes a chloride channel that regulates fluid balance in the lungs and digestive tract. Loss-of-function variants cause cystic fibrosis, a recessive disease affecting mucus, sweat, and digestion.",
    "HTT": "HTT encodes huntingtin, a protein essential for normal brain development. An expanded CAG repeat causes Huntington's disease, a progressive and currently incurable neurodegenerative disorder.",
    "APOE": "APOE encodes apolipoprotein E, a key regulator of cholesterol and lipid transport in blood and brain. The ε4 allele is the strongest known genetic risk factor for late-onset Alzheimer's disease.",
    "LDLR": "LDLR encodes the LDL receptor, which clears cholesterol-carrying particles from the bloodstream. Variants that impair the receptor cause familial hypercholesterolemia and early-onset heart disease.",
    "RB1": "RB1 encodes the retinoblastoma protein, a critical brake on the cell cycle that prevents uncontrolled growth. Biallelic loss causes retinoblastoma, a childhood eye tumor, and predisposes to other cancers.",
    "APC": "APC encodes a tumor suppressor that regulates the Wnt cell-growth signaling pathway. Germline variants cause familial adenomatous polyposis (FAP), with hundreds of colon polyps and near-certain colorectal cancer risk.",
    "MLH1": "MLH1 encodes a DNA mismatch repair protein that corrects replication errors before they become mutations. Variants cause Lynch syndrome, the most common hereditary colorectal cancer predisposition syndrome.",
    "PTEN": "PTEN encodes a phosphatase that suppresses the PI3K/AKT cell-growth signaling pathway. Loss of function is linked to Cowden syndrome and is among the most frequent alterations in sporadic cancers.",
    "VHL": "VHL encodes a protein that regulates HIF-1α and the cellular response to low oxygen levels. Germline variants cause von Hippel-Lindau disease, predisposing to kidney cancer and brain hemangioblastomas.",
    "NF1": "NF1 encodes neurofibromin, which suppresses RAS-driven cell proliferation signals. Variants cause neurofibromatosis type 1, characterized by benign nerve-sheath tumors, skin markings, and learning differences.",
    "NF2": "NF2 encodes merlin, a tumor suppressor in Schwann and meningeal cells. Variants cause neurofibromatosis type 2, presenting with bilateral vestibular schwannomas and progressive hearing loss.",
    "TSC1": "TSC1 encodes hamartin, which together with TSC2 suppresses mTOR-driven cell growth. Variants cause tuberous sclerosis complex, with benign tumors forming in the brain, skin, kidneys, and lungs.",
    "TSC2": "TSC2 encodes tuberin, the catalytic partner of hamartin in the mTOR-suppressing complex. Variants cause tuberous sclerosis complex with multi-organ benign tumor formation.",
    "PKD1": "PKD1 encodes polycystin-1, a membrane protein involved in cell-to-cell signaling in kidney tubules. Variants account for ~85% of autosomal dominant polycystic kidney disease (ADPKD), the most common inherited kidney disorder.",
    "PKD2": "PKD2 encodes polycystin-2, an ion channel that partners with polycystin-1. Variants cause the remaining ~15% of ADPKD, generally with a milder course.",
    "HEXA": "HEXA encodes the alpha subunit of beta-hexosaminidase A, which breaks down GM2 gangliosides in lysosomes. Biallelic variants cause Tay-Sachs disease, a fatal infantile neurological disorder.",
    "PAH": "PAH encodes phenylalanine hydroxylase, which converts phenylalanine to tyrosine in the liver. Biallelic variants cause phenylketonuria (PKU); untreated, excess phenylalanine causes severe intellectual disability.",
    "MUTYH": "MUTYH encodes a base-excision repair enzyme that corrects oxidative DNA damage. Biallelic variants cause MUTYH-associated polyposis, increasing colorectal cancer risk.",
    "KRAS": "KRAS encodes a GTPase central to growth factor signaling from cell-surface receptors. Somatic KRAS mutations are among the most frequent in human cancers, including pancreatic, colorectal, and lung adenocarcinoma.",
}

@st.cache_data(ttl=60*60*24, show_spinner=False)
def fetch_ncbi_gene_summary(gene_symbol: str) -> str | None:
    try:
        search_params = {
            "db": "gene",
            "term": f"{gene_symbol}[Gene Name] AND Homo sapiens[Organism]",
            "retmode": "json", "retmax": "1",
        }
        if NCBI_TOOL: search_params["tool"] = NCBI_TOOL
        if NCBI_EMAIL: search_params["email"] = NCBI_EMAIL
        r = requests.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
                         params=search_params, timeout=10)
        r.raise_for_status()
        ids = r.json().get("esearchresult", {}).get("idlist", [])
        if not ids:
            return None
        fetch_params = {"db": "gene", "id": ids[0], "retmode": "json"}
        if NCBI_TOOL: fetch_params["tool"] = NCBI_TOOL
        if NCBI_EMAIL: fetch_params["email"] = NCBI_EMAIL
        r2 = requests.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",
                          params=fetch_params, timeout=10)
        r2.raise_for_status()
        result = r2.json().get("result", {})
        summary = result.get(ids[0], {}).get("summary", "")
        if summary:
            sents = re.split(r"(?<=[.!?])\s+", summary.strip())
            return " ".join(sents[:2])
        return None
    except Exception:
        return None

# ----------------------------
# PubMed helpers
# ----------------------------
@st.cache_data(ttl=60*60, show_spinner=False)
def _pubmed_one_sentence(pmid: str) -> str | None:
    if not pmid:
        return None

    params = {"db": "pubmed", "id": pmid, "retmode": "xml"}
    if NCBI_TOOL: params["tool"] = NCBI_TOOL
    if NCBI_EMAIL: params["email"] = NCBI_EMAIL

    try:
        r = requests.get(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
            params=params,
            timeout=10
        )
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

        sent = re.split(r"(?<=[.!?])\s+", abstract)[0].strip()
        return sent or None

    except Exception:
        return None

@st.cache_data(ttl=60*60, show_spinner=False)
def _guess_pmid(gene: str | None, mutation: str | None, disease: str | None) -> str | None:
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
            r = requests.get(
                "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
                params=params,
                timeout=10
            )
            r.raise_for_status()
            ids = r.json().get("esearchresult", {}).get("idlist", [])
            if ids:
                return ids[0]
        except Exception:
            pass

    return None

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

def _pick_first(df: pd.DataFrame, candidates: List[str]) -> Optional[str]:
    for c in candidates:
        if c in df.columns:
            return c
    return None

def _normalize_gene(x: str) -> str:
    return str(x).strip().upper()

def _fallback_sentence(gene: str, mutation_label: Optional[str], disease: Optional[str], significance: Optional[str]) -> str:
    bits = []
    if gene:
        bits.append(gene)
    if mutation_label:
        bits.append(f"variant {mutation_label}")
    if disease:
        bits.append(f"is associated with {disease}")
    if significance and str(significance).strip():
        bits.append(f"({significance})")

    txt = " ".join(bits).strip()
    if not txt:
        txt = "Variant associated information not available."
    if not txt.endswith("."):
        txt += "."
    return txt

# ----------------------------
# Variant cards builder
# ----------------------------
def build_variant_cards(df: pd.DataFrame, gene: str, n: int = 3) -> list[dict]:
    gene_norm = str(gene).strip().upper()
    df = df.copy()
    df.columns = [c.strip() for c in df.columns]

    sig_col  = next((c for c in ["ClinicalSignificance","clinical_significance","CLNSIG","Significance","CLIN_SIG"] if c in df.columns), None)
    name_col = next((c for c in POSSIBLE_NAME_COLS if c in df.columns), None)
    gene_col = next((c for c in POSSIBLE_GENE_COLS if c in df.columns), None)

    # filter by gene
    if gene_col:
        df["_gene_norm"] = df[gene_col].fillna("").map(lambda x: str(x).strip().upper())
        hits = df[df["_gene_norm"] == gene_norm]
    else:
        hits = df

    # fallback: parse gene in Name column
    if hits.empty and name_col:
        mask = df[name_col].fillna("").str.contains(rf"\({re.escape(gene_norm)}\)", regex=True, na=False)
        hits = df[mask]

    hits = hits.head(max(1, n))

    cards = []
    for _, row in hits.iterrows():
        # mutation label
        mut_label = None
        for c in POSSIBLE_HGVSC_COLS + POSSIBLE_HGVSP_COLS:
            if c in row and pd.notna(row[c]) and str(row[c]).strip():
                mut_label = str(row[c]).strip()
                break
        if not mut_label and name_col and pd.notna(row.get(name_col)):
            mut_label = str(row[name_col]).strip()

        # disease
        disease = None
        for c in POSSIBLE_DISEASE_COLS:
            if c in row and pd.notna(row[c]):
                raw = str(row[c]).strip()
                if raw:
                    first = re.split(r"[;,\|]", raw)[0].strip()
                    if first and first.lower() not in {"not provided","not specified"}:
                        disease = first
                        break

        significance = str(row.get(sig_col, "") or "").strip() if sig_col else ""

        pmid = _parse_first_pmid(row) or ""
        if not pmid:
            pmid = _guess_pmid(gene_norm, mut_label, disease) or ""

        summary = _pubmed_one_sentence(pmid) if pmid else None
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

# ----------------------------
# Data loading + normalization
# ----------------------------
DATA_PATH = Path(__file__).parent / "clinvar_sample.csv"

@st.cache_data(show_spinner=False)
def load_variants(csv_path: Path) -> pd.DataFrame:
    if not csv_path.exists():
        return pd.DataFrame(columns=["gene","variant_id","protein_change","cdna_change",
                                     "clinical_significance","condition","source","PMID"])
    if csv_path.suffix.lower() in {".tsv", ".txt"}:
        return pd.read_csv(csv_path, sep="\t", dtype=str, low_memory=False)
    return pd.read_csv(csv_path, dtype=str, low_memory=False)

def canonicalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df.columns = [c.strip() for c in df.columns]

    gene_col = _pick_first(df, ["gene"] + POSSIBLE_GENE_COLS)
    if gene_col:
        df["gene"] = df[gene_col].fillna("").map(_normalize_gene)

    sig_col = _pick_first(df, ["clinical_significance","ClinicalSignificance","CLNSIG","Significance","CLIN_SIG"])
    if sig_col and "clinical_significance" not in df.columns:
        df["clinical_significance"] = df[sig_col]

    cond_col = _pick_first(df, ["condition"] + POSSIBLE_DISEASE_COLS)
    if cond_col and "condition" not in df.columns:
        df["condition"] = df[cond_col]

    vid_col = _pick_first(df, ["variant_id","VariantID","VariationID","VCV","VCV_ID"])
    if vid_col and "variant_id" not in df.columns:
        df["variant_id"] = df[vid_col]

    prot_col = _pick_first(df, ["protein_change","HGVSp","Protein_change","HGVS_p","hgvs_p","ProteinChange"])
    if prot_col and "protein_change" not in df.columns:
        df["protein_change"] = df[prot_col]

    cdna_col = _pick_first(df, ["cdna_change","HGVSc","HGVS_cDNA","HGVS_c","hgvs_c"])
    if cdna_col and "cdna_change" not in df.columns:
        df["cdna_change"] = df[cdna_col]

    if "source" not in df.columns:
        df["source"] = "ClinVar"

    if "PMID" not in df.columns:
        df["PMID"] = ""

    return df

# ----------------------------
# Page config + styles
# ----------------------------
st.set_page_config(page_title="Gene Variant Explorer", page_icon="🧬", layout="wide")

st.markdown("""
<style>
.hero {
  padding: 2rem; border-radius: 16px;
  background: linear-gradient(180deg, #f5f6ff 0%, #ffffff 100%);
  border: 1px solid #e9ecff;
}
.kidcard { border:1px solid #e6e8ff; border-radius:14px; padding:1rem; background:#fafbff; }
.smallmuted { color: rgba(0,0,0,0.6); font-size: 0.95rem; }
.gene-info-box {
  padding: 1rem 1.2rem; border-radius: 12px;
  background: #f0f4ff; border-left: 4px solid #4f46e5;
  margin-bottom: 1rem;
}
</style>
""", unsafe_allow_html=True)

# ----------------------------
# Session state init (for gene + filters)
# ----------------------------
if "loaded_params" not in st.session_state:
    st.session_state["loaded_params"] = True
    st.session_state["query"] = (st.query_params.get("gene", "") or "").upper()
    st.session_state["selected_sigs"] = _get_param_list("sig")
    st.session_state["selected_conditions"] = _get_param_list("cond")

if "search_history" not in st.session_state:
    st.session_state["search_history"] = []

# ----------------------------
# TOP NAV + PAGE ROUTING
# ----------------------------
def _get_page_from_url() -> str:
    p = st.query_params.get("page", "home")
    if isinstance(p, list):
        p = p[0] if p else "home"
    return (p or "home").lower()

def _set_page_in_url(page: str):
    st.query_params["page"] = page

PAGE_LABELS = {
    "home": "Home",
    "tool": "Tool",
    "impact": "Impact",
    "about": "About",
    "feedback": "Feedback",
}

current_page = _get_page_from_url()
if current_page not in PAGE_LABELS:
    current_page = "home"

# top nav
nav_cols = st.columns([1,1,1,1,1,6])
nav_keys = ["home", "tool", "impact", "about", "feedback"]
for i, key in enumerate(nav_keys):
    label = PAGE_LABELS[key]
    if key == current_page:
        nav_cols[i].markdown(f"**{label}**")
    else:
        if nav_cols[i].button(label, key=f"nav_{key}"):
            _set_page_in_url(key)
            st.rerun()
st.markdown("---")

# ----------------------------
# Load dataframe once
# ----------------------------
variants_df = canonicalize_columns(load_variants(DATA_PATH))

# ----------------------------
# Pages
# ----------------------------
def render_home():
    st.title("🧬 Gene Variant Explorer")
    st.markdown("<div class='smallmuted'>A student-built tool to make genetics understandable using real public scientific data.</div>", unsafe_allow_html=True)

    # Add quick start guide
    with st.expander("Quick Start Guide (First Time Users)", expanded=True):
        st.markdown("""
        1. **Click "Tool" in the navigation bar above**
        2. **Try searching for "BRCA1"** - a well-studied gene
        3. **Look at the results table** - each row is a different variant
        4. **Click "Generate cards"** at the bottom to see easy-to-read summaries
        5. **Use the sidebar filters** to narrow down results by disease or significance

        **Don't understand the terms?** Check out the "About" page for explanations!
        """)

    st.markdown("""
<div class="hero">
<h3>What this tool does</h3>
<ul>
  <li>Search a gene (example: <b>BRCA1</b>, <b>TP53</b>) and view real variants.</li>
  <li>See associated conditions from ClinVar (NIH).</li>
  <li>Read a <b>one-sentence research summary</b> pulled from PubMed abstracts.</li>
  <li>Built for <b>learning</b> — not diagnosis.</li>
</ul>

<h3>How to use it</h3>
<ol>
  <li>Click <b>Tool</b> in the top bar.</li>
  <li>Type a gene like <b>BRCA1</b>.</li>
  <li>Filter by condition or clinical significance.</li>
  <li>Generate cards for readable summaries.</li>
</ol>

<p><b>Disclaimer:</b> Educational only. Not medical advice.</p>
</div>
""", unsafe_allow_html=True)

    st.markdown("### Data sources")
    st.write("- **ClinVar (NIH)** for variants + conditions")
    st.write("- **PubMed (NCBI)** for research abstracts")

    st.info("Ready? Click **Tool** in the top navigation bar.")


def render_about():
    st.title("About")
    st.markdown("""
### Purpose
Genetics information is often written for professionals and can be confusing for people with little background in genetics. This app helps students and families explore variant evidence responsibly.

### What's inside
- Variant rows from a ClinVar-style file
- PubMed summaries fetched via NCBI E-utilities
- Shareable links that preserve filters

### How to read genetics information

#### Understanding the table columns:
- **Gene**: The gene symbol (e.g., BRCA1, TP53) - a section of DNA that contains instructions for making a protein
- **Variant ID**: A unique identifier for this specific variant (like a serial number)
- **Protein change** (p.Val1736Ala): Describes how the protein is changed
  - Format: `p.OriginalAminoAcidPositionNewAminoAcid`
  - Example: `p.Val1736Ala` means at position 1736, Valine (Val) was replaced with Alanine (Ala)
  - Think of it like changing one letter in a word: "cat" → "bat"
- **cDNA change** (c.5207T>C): Describes the DNA change in the coding sequence
  - Format: `c.PositionOriginalBase>NewBase`
  - Example: `c.5207T>C` means at position 5207, Thymine (T) was replaced with Cytosine (C)
  - This is the DNA-level change that causes the protein change
- **Clinical significance**: What this variant means for health
  - **Pathogenic**: ⚠️ Known to cause disease
  - **Likely Pathogenic**: ⚠️ Strongly suspected to cause disease
  - **Benign**: ✅ Does not cause disease
  - **Likely Benign**: ✅ Unlikely to cause disease
  - **VUS (Variant of Uncertain Significance)**: ❓ Unknown effect - needs more research
- **Condition**: The disease or health condition associated with this variant
- **Source**: Where the data comes from (usually ClinVar, a NIH database)
- **PMID**: Link to the research paper (PubMed ID) - click to read the original study

#### Understanding variant notation:
- **c.** = coding DNA sequence (the part that makes proteins)
- **p.** = protein sequence (the actual protein building blocks)
- **fs** = frameshift (shifts how the code is read, often causing major problems)
- **del** = deletion (DNA was removed)
- **ins** = insertion (DNA was added)
- *** or Ter** = stop codon (protein stops being made early, usually problematic)

#### Example: Understanding a variant
**Example variant:** `p.Val1736Ala` in BRCA1

**What this means:**
- In the BRCA1 gene, at position 1736
- The amino acid **Valine (Val)** was changed to **Alanine (Ala)**
- This is like changing one letter in a word: "cat" → "bat"

**Why it matters:**
- This small change can affect how the protein works
- If the protein doesn't work correctly, it may not protect against cancer
- That's why this variant might be "Pathogenic" (disease-causing)

### Glossary
- **Gene**: A section of DNA that contains instructions for making a protein. Genes are like recipes in a cookbook.
- **Variant**: A change in the DNA sequence. Like a typo in a recipe - sometimes it matters, sometimes it doesn't.
- **Pathogenic**: A variant that causes disease. Think of it as a 'broken' gene that doesn't work correctly.
- **Benign**: A variant that does NOT cause disease. The gene still works fine despite the change.
- **VUS (Variant of Uncertain Significance)**: We don't know yet if this variant causes disease. More research is needed.
- **BRCA1/BRCA2**: Genes that help repair DNA damage. Variants in these genes increase cancer risk.
- **TP53**: Often called the "guardian of the genome" - helps prevent cancer by stopping damaged cells.

### Privacy
This app does not collect personal health information.
""")


def render_impact():
    st.title("Impact")
    st.markdown("""
Show proof that your project helped people.

### What will be put
- Number of users (visits)
- Number of gene searches
- CSV downloads
- Feedback submissions
- Workshops / presentations delivered
- People reached

### Plan for impact
1. Present it to one real group (club/class/library).
2. Collect feedback with the Feedback tab.
3. Improve the tool based on feedback.
4. Report numbers here.
""")


def render_feedback():
    st.title("Feedback")
    st.write("Help improve the tool. What should we add or fix?")

    with st.form("feedback_form", clear_on_submit=True):
        name = st.text_input("Your name (optional)")
        role = st.selectbox("I am a...", ["Student", "Parent", "Teacher", "Other"])
        msg  = st.text_area("What should we improve?")
        submitted = st.form_submit_button("Send")

    if submitted:
        ts = datetime.utcnow().isoformat()
        row = pd.DataFrame([{"timestamp": ts, "name": name, "role": role, "message": msg}])
        try:
            log_path = Path("feedback_log.csv")
            if log_path.exists():
                row.to_csv(log_path, mode="a", header=False, index=False)
            else:
                row.to_csv(log_path, index=False)
            st.success("Thanks! Your feedback was recorded.")
        except Exception:
            st.info("Thanks! (Logging failed on this host.)")


def filter_variants(df: pd.DataFrame, gene_text: str) -> pd.DataFrame:
    if not gene_text:
        return df.iloc[0:0]
    exact = df[df["gene"] == gene_text]
    if len(exact) > 0:
        return exact
    starts = df[df["gene"].str.startswith(gene_text, na=False)]
    if len(starts) > 0:
        return starts
    contains = df[df["gene"].str.contains(gene_text, na=False)]
    return contains

def color_significance(val):
    if pd.isna(val):
        return ''
    val_str = str(val).lower()
    if 'pathogenic' in val_str and 'likely' not in val_str:
        return 'background-color: #ffebee'
    elif 'likely pathogenic' in val_str:
        return 'background-color: #fff3e0'
    elif 'benign' in val_str and 'likely' not in val_str:
        return 'background-color: #e8f5e9'
    elif 'likely benign' in val_str:
        return 'background-color: #f1f8e9'
    elif 'uncertain' in val_str or 'vus' in val_str:
        return 'background-color: #f3e5f5'
    return ''

def _normalize_sig_category(val: str) -> str:
    v = str(val).lower()
    if "likely pathogenic" in v:
        return "Likely Pathogenic"
    if "pathogenic" in v and "conflicting" not in v:
        return "Pathogenic"
    if "likely benign" in v:
        return "Likely Benign"
    if "benign" in v:
        return "Benign"
    if "uncertain" in v or "vus" in v:
        return "VUS"
    if "conflicting" in v:
        return "Conflicting"
    return "Other"


def render_tool():
    # Sidebar only on tool page
    with st.sidebar:
        st.header("Filters")
        st.caption("Use these to narrow down your search results")

        sig_options  = sorted(variants_df["clinical_significance"].dropna().unique().tolist()) if "clinical_significance" in variants_df.columns else []
        cond_options = sorted(variants_df["condition"].dropna().unique().tolist()) if "condition" in variants_df.columns else []

        selected_sigs = st.multiselect(
            "Clinical significance",
            options=sig_options,
            default=[s for s in st.session_state.get("selected_sigs", []) if s in sig_options],
            help="Filter by whether variants are disease-causing (Pathogenic) or not (Benign)"
        )
        selected_conditions = st.multiselect(
            "Condition",
            options=cond_options,
            default=[c for c in st.session_state.get("selected_conditions", []) if c in cond_options],
            help="Filter by specific diseases or health conditions"
        )

        st.session_state["selected_sigs"] = selected_sigs
        st.session_state["selected_conditions"] = selected_conditions

        st.markdown("---")

        # Recent searches
        history = st.session_state.get("search_history", [])
        if history:
            st.caption("Recent searches")
            for past_gene in history:
                if st.button(past_gene, key=f"hist_{past_gene}"):
                    st.session_state["query"] = past_gene
                    _set_params_keep_page(past_gene, st.session_state.get("selected_sigs", []), st.session_state.get("selected_conditions", []))
                    st.rerun()
            st.markdown("---")

        # Teacher mode
        teacher_mode = st.checkbox(
            "Teacher mode",
            value=False,
            help="Shows iframe embed code so you can embed this search in slides or a class website"
        )

        st.markdown("---")
        st.caption("Data source: ClinVar subset. PubMed summaries via NCBI E-utilities.")

    st.title("🧬 Gene Variant Explorer — Tool")
    st.write("Type a gene symbol (e.g., **BRCA1**) to see related variants. Use the sidebar to filter.")

    # Add help section
    with st.expander("📚 What do these terms mean?", expanded=False):
        st.markdown("""
        - **Protein change (p.Val1736Ala)**: How the protein building blocks are altered. The format shows which amino acid was changed at which position.
        - **cDNA change (c.5207T>C)**: The DNA sequence change in the gene. Shows which DNA base was changed.
        - **Clinical significance**: Whether this variant is known to cause disease
          - **Pathogenic**: ⚠️ Disease-causing
          - **Likely Pathogenic**: ⚠️ Strongly suspected to cause disease
          - **Benign**: ✅ Not disease-causing
          - **Likely Benign**: ✅ Unlikely to cause disease
          - **VUS**: ❓ Unknown significance (needs more research)
        - **PMID**: Research paper identifier - click to read the original study
        """)
        st.info("💡 For more detailed explanations, check out the **About** page!")

    # Add example section
    with st.expander("💡 Example: Understanding a variant", expanded=False):
        st.markdown("""
        **Example variant:** `p.Val1736Ala` in BRCA1

        **What this means:**
        - In the BRCA1 gene, at position 1736
        - The amino acid **Valine (Val)** was changed to **Alanine (Ala)**
        - This is like changing one letter in a word: "cat" → "bat"

        **Why it matters:**
        - This small change can affect how the protein works
        - If the protein doesn't work correctly, it may not protect against cancer
        - That's why this variant might be "Pathogenic" (disease-causing)
        """)

    # File uploader (optional)
    st.markdown("#### Try your own file")
    up = st.file_uploader(
        "Upload a ClinVar-style CSV/TSV",
        type=["csv","tsv","txt"],
        help="We look for gene, cdna_change/protein_change, clinical_significance, condition, PMID (optional)"
    )

    # Use uploaded data if provided
    df_live = variants_df
    if up is not None:
        sep = "\t" if (up.name.endswith(".tsv") or up.name.endswith(".txt")) else ","
        user_df = pd.read_csv(up, sep=sep, dtype=str, low_memory=False)
        df_live = canonicalize_columns(user_df)
        st.success(f"Loaded {len(df_live)} rows from your file.")

    # Search input
    col1, col2 = st.columns([2,1])
    with col1:
        query = st.text_input(
            "Gene symbol",
            value=st.session_state.get("query", ""),
            placeholder="e.g., BRCA1, TP53, BRCA2"
        ).strip().upper()
        st.session_state["query"] = query
    with col2:
        topk = st.number_input("Max results", min_value=5, max_value=5000, value=200, step=5)

    # Keep URL share params updated
    _set_params_keep_page(
        st.session_state.get("query",""),
        st.session_state.get("selected_sigs", []),
        st.session_state.get("selected_conditions", [])
    )

    # Landing
    if not query:
        st.markdown('<div class="hero">', unsafe_allow_html=True)
        st.markdown("## Welcome!")
        st.write("**New to genetics?** Start by searching for a well-known gene:")

        gene_descriptions = {
            "BRCA1": "Breast cancer gene 1 - helps repair DNA",
            "BRCA2": "Breast cancer gene 2 - also helps repair DNA",
            "TP53": "Tumor protein 53 - 'guardian of the genome'",
            "APC": "Adenomatous polyposis coli - involved in colon cancer",
            "MLH1": "DNA mismatch repair gene - helps fix DNA errors"
        }

        cols = st.columns(5)
        picks = ["BRCA1", "BRCA2", "TP53", "APC", "MLH1"]
        for i, g in enumerate(picks):
            with cols[i]:
                if st.button(g, key=f"pick_{g}"):
                    st.session_state["query"] = g
                    _set_params_keep_page(g, st.session_state.get("selected_sigs", []), st.session_state.get("selected_conditions", []))
                    st.rerun()
                st.caption(gene_descriptions[g])

        st.markdown("---")
        st.markdown("** Tip:** After searching, click 'Generate cards' to see easy-to-read summaries!")
        st.markdown("</div>", unsafe_allow_html=True)
        return

    # Fuzzy suggestions
    gene_options = sorted(df_live["gene"].dropna().unique().tolist()) if "gene" in df_live.columns else []
    if query and (query not in df_live["gene"].unique()):
        suggestions = fuzzy_gene_candidates(query, gene_options, limit=5)
        if suggestions:
            st.info("Did you mean:")
            cols = st.columns(len(suggestions))
            for i, s in enumerate(suggestions):
                if cols[i].button(s, key=f"sug_{s}"):
                    st.session_state["query"] = s
                    _set_params_keep_page(s, st.session_state.get("selected_sigs", []), st.session_state.get("selected_conditions", []))
                    st.rerun()

    # Gene summary card
    gene_blurb = GENE_SUMMARIES.get(query)
    if gene_blurb is None:
        with st.spinner(f"Looking up {query}..."):
            gene_blurb = fetch_ncbi_gene_summary(query)
    if gene_blurb:
        st.markdown(
            f'<div class="gene-info-box"><b>{query}</b> — {gene_blurb}</div>',
            unsafe_allow_html=True
        )

    # Filter + apply sidebar filters
    results = filter_variants(df_live, query)
    if st.session_state.get("selected_sigs"):
        results = results[results["clinical_significance"].isin(st.session_state["selected_sigs"])]
    if st.session_state.get("selected_conditions"):
        results = results[results["condition"].isin(st.session_state["selected_conditions"])]

    results = results.head(int(topk))

    # Update search history
    if query and not results.empty:
        history = st.session_state.get("search_history", [])
        if query not in history:
            history = [query] + history
        st.session_state["search_history"] = history[:10]

    # Metrics
    left, right = st.columns(2)
    with left:
        st.metric("Matches", len(results))
    with right:
        st.metric("Unique conditions", results["condition"].nunique() if not results.empty else 0)

    # Pathogenicity breakdown chart
    if not results.empty and "clinical_significance" in results.columns:
        cats = results["clinical_significance"].dropna().map(_normalize_sig_category)
        order = ["Pathogenic", "Likely Pathogenic", "VUS", "Conflicting", "Likely Benign", "Benign", "Other"]
        counts = cats.value_counts().reindex(order, fill_value=0)
        counts = counts[counts > 0]
        if not counts.empty:
            with st.expander("Pathogenicity breakdown", expanded=True):
                st.bar_chart(counts)
                st.caption("Distribution of clinical significance categories across the variants shown.")

    # Chromosome location banner (only when all variants are on one chromosome)
    chrom_col = _pick_first(df_live, POSSIBLE_CHROM_COLS)
    pos_col   = _pick_first(df_live, POSSIBLE_POS_COLS)
    if not results.empty and chrom_col and pos_col:
        chroms = results[chrom_col].dropna().unique()
        if len(chroms) == 1:
            pos_vals = results[pos_col].dropna()
            if not pos_vals.empty:
                try:
                    st.caption(f"📍 Chromosome {chroms[0]}, position {int(float(pos_vals.iloc[0])):,}")
                except (ValueError, TypeError):
                    pass

    # Share section
    with st.expander("Share"):
        st.write("Copy the link from your browser address bar to share this search + filters.")

    # Teacher mode embed code
    if teacher_mode:
        with st.expander("Embed this search (Teacher mode)", expanded=True):
            embed_params = f"?embed=true&page=tool&gene={query}"
            if st.session_state.get("selected_sigs"):
                embed_params += "&sig=" + "|".join(st.session_state["selected_sigs"])
            iframe_code = (
                f'<iframe src="https://your-app.streamlit.app{embed_params}" '
                f'width="100%" height="700" frameborder="0"></iframe>'
            )
            st.code(iframe_code, language="html")
            st.caption(
                "Replace `your-app.streamlit.app` with your deployment URL. "
                "The `?embed=true` parameter hides the Streamlit toolbar for a cleaner classroom view."
            )

    # Table
    if results.empty:
        st.warning("No variants found. Try removing filters or check the gene symbol.")
        return

    st.markdown("""
    **Understanding the results table:**
    - Each row represents a different variant (change) in the gene
    - Color coding: Red/Orange = disease-causing, Green = not disease-causing, Purple = unknown
    - Click on PMID numbers to read the original research paper
    - Use filters in the sidebar to narrow results
    """)

    preferred_cols = ["gene", "variant_id", "protein_change", "cdna_change",
                      "clinical_significance", "condition", "source", "PMID"]
    ordered_cols = [c for c in preferred_cols if c in results.columns] + [c for c in results.columns if c not in preferred_cols]

    if "clinical_significance" in ordered_cols:
        styled_df = results[ordered_cols].style.map(color_significance, subset=['clinical_significance'])
        st.dataframe(styled_df, width="stretch", hide_index=True)
    else:
        st.dataframe(results[ordered_cols], width="stretch", hide_index=True)

    csv_bytes = results[ordered_cols].to_csv(index=False).encode("utf-8")
    st.download_button(
        label="Download results (CSV)",
        data=csv_bytes,
        file_name=f"{query}_variants.csv",
        mime="text/csv",
    )

    # Cards
    with st.expander("🧬 Variant cards with PubMed summary", expanded=True):
        c1, c2 = st.columns([1,5])
        with c1:
            n_cards = st.number_input("How many", min_value=1, max_value=5, value=3, step=1)
        with c2:
            st.caption("Each card includes a one-sentence summary extracted from PubMed when available.")

        if st.button("Generate cards"):
            cards = build_variant_cards(df_live, query, n=int(n_cards))
            gene_context = GENE_SUMMARIES.get(query, "")
            for i, c in enumerate(cards, 1):
                st.markdown("---")
                st.markdown(f"**Variant {i}: {c['Mutation']}**")

                if c['Mutation'].startswith('p.'):
                    st.caption("💡 'p.' means this describes a change in the protein (the building blocks that do the work)")
                elif c['Mutation'].startswith('c.'):
                    st.caption("💡 'c.' means this describes a change in the DNA code")

                st.markdown(f"- **Disease/Phenotype:** {c['Disease/Phenotype']}")

                sig = c['Clinical significance']
                sig_explanation = {
                    'Pathogenic': '⚠️ This variant is known to cause disease',
                    'Likely Pathogenic': '⚠️ This variant is strongly suspected to cause disease',
                    'Benign': '✅ This variant does not cause disease',
                    'Likely Benign': '✅ This variant is unlikely to cause disease',
                    'VUS': '❓ The significance of this variant is uncertain - more research is needed',
                    'Variant of Uncertain Significance': '❓ The significance of this variant is uncertain - more research is needed'
                }.get(sig, '')

                if sig_explanation:
                    st.markdown(f"- **Clinical significance:** {sig} — {sig_explanation}")
                else:
                    st.markdown(f"- **Clinical significance:** {sig}")

                if c["PMID"] != "—":
                    st.markdown(f"- **Research paper:** [PMID {c['PMID']}](https://pubmed.ncbi.nlm.nih.gov/{c['PMID']}/) (click to read the original study)")
                else:
                    st.caption("No research paper found for this variant.")

                st.markdown(f"**Summary:** {c['Summary']}")

                # Why does this matter?
                sig_lower = sig.lower()
                disease = c["Disease/Phenotype"]
                why_parts = []
                if "pathogenic" in sig_lower and "benign" not in sig_lower:
                    if disease and disease != "Not specified":
                        why_parts.append(f"This variant is classified as **{sig}**, meaning it is linked to an increased risk of **{disease}**.")
                    else:
                        why_parts.append(f"This variant is classified as **{sig}**, meaning it may contribute to disease.")
                    if gene_context:
                        why_parts.append(f"**About {query}:** {gene_context}")
                elif "benign" in sig_lower:
                    why_parts.append(f"This variant is classified as **{sig}** — the {query} gene still functions normally despite this DNA change.")
                elif "uncertain" in sig_lower or "vus" in sig_lower:
                    why_parts.append("This variant's impact is **uncertain** (Variant of Uncertain Significance). Scientists need more evidence to determine whether it causes disease.")
                    if disease and disease != "Not specified":
                        why_parts.append(f"It has been observed in people with {disease}, but this association alone does not confirm it as the cause.")

                if why_parts:
                    with st.expander("Why does this matter?"):
                        for part in why_parts:
                            st.markdown(part)

# ----------------------------
# Render selected page
# ----------------------------
if current_page == "home":
    render_home()
elif current_page == "tool":
    render_tool()
elif current_page == "impact":
    render_impact()
elif current_page == "about":
    render_about()
elif current_page == "feedback":
    render_feedback()

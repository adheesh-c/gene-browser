# Data Sources & Provenance

**Educational use only. Not medical advice.** The bundled dataset is a small,
hand-curated *teaching* set. It is **not** a complete or clinical-grade copy of
ClinVar and must never be used for diagnosis, personal risk assessment, or any
clinical decision.

## What's in `clinvar_sample.csv`

A curated teaching dataset of **218 variants across 25 genes**, built to give the
GeneScope tool a realistic feel while staying small and correct. Every gene has a
**spread of clinical significance** (at least one Pathogenic, one Uncertain, and
one Benign row) so the pathogenicity chart and filters are meaningful.

| Column | Meaning |
|---|---|
| `gene` | HGNC gene symbol |
| `variant_id` | ClinVar Variation accession (`VCV…`) |
| `protein_change` | HGVS protein change (`p.…`), blank for intronic/non-coding variants |
| `cdna_change` | HGVS coding-DNA change (`c.…`) |
| `clinical_significance` | ClinVar aggregate germline classification |
| `condition` | Associated condition/trait (ClinVar `trait_name`) |
| `source` | Always `ClinVar` for this set |
| `PMID` | A real PubMed ID associated with the variant/gene (see below) |
| `collection` | Themed grouping for one-click starting points |

### Collections
- **brain_tumor** — IDH1, EGFR, ATRX, TP53 (ties the tool to the aptamer / brain-tumor theme)
- **hereditary_cancer** — BRCA1, BRCA2, TP53*, RB1, APC, MLH1, PTEN, VHL, MUTYH, KRAS
- **neuro** — HTT, APOE, NF1, NF2, TSC1, TSC2
- **metabolic** — CFTR, LDLR, HEXA, PAH
- **renal** — PKD1, PKD2

\* `collection` holds one value per row. **TP53** is placed in `brain_tumor` (it is
central to gliomas and Li-Fraumeni syndrome, which prominently includes brain
tumors). It remains fully searchable by gene symbol and is a hereditary-cancer
gene as well; the single-valued `collection` column just controls which themed
starting button surfaces it.

## Retrieval

- **Source database:** NCBI **ClinVar**, queried live via the E-utilities API
  (`esearch` + `esummary`, `db=clinvar`).
- **ClinVar retrieval date:** **2026-08-14**
- **Classifications** come straight from each record's aggregate
  `germline_classification` (ClinVar's own review-status-weighted call). Records
  with *"no assertion criteria provided"* review status, and conditions of
  *"not provided" / "not specified"*, were filtered out; better-reviewed records
  (expert panel / multiple submitters) were preferred. Nothing was relabeled by hand.
- For each gene, records were pulled per significance bucket using the ClinVar
  property filters `"clinsig pathogenic"`, `"clinsig benign"`, and `"clinsig vus"`
  to guarantee a spread.
- **Biologically honest gaps:** a few genes have fewer rows in a bucket because
  ClinVar genuinely has few well-reviewed coding variants there — e.g. **IDH1**
  (an oncogene hotspot; its benign entries here are real intronic variants) and
  **APOE**. **HTT**'s hallmark pathogenic mechanism is a *CAG repeat expansion*,
  not a single-nucleotide variant, so it is not represented as a discrete SNV
  row; a real pathogenic HTT coding variant is included instead, and the repeat
  expansion is described in the app's gene summary.

## PMIDs — how they were chosen and verified

- **Primary:** each variant's PMID is a real PubMed citation **linked to that
  ClinVar record** (`elink`, `clinvar → pubmed`). Because ClinVar often links the
  same high-throughput functional study (e.g. saturation-genome-editing papers)
  to many variants of a gene, the 218 rows reference **46 distinct PMIDs** — a
  PMID may repeat across several variants of the same gene. Each is genuinely
  associated with the variant in ClinVar.
- **Fallback (5 rows):** five variants (BRCA2 ×2, APOE ×1, LDLR ×2) had **no**
  ClinVar-linked citation. Each was given a real, gene-level PubMed reference
  found via `esearch` (`db=pubmed`): BRCA2 → PMID `29446198`, APOE → `37957317`,
  LDLR → `31491741`. These are gene references, not necessarily variant-specific.
- **Verification:** **every** PMID in the file was re-checked against PubMed
  `esummary` on **2026-08-14** and confirmed to resolve (return a title). No PMID
  was invented; any row whose PMID could not be verified was dropped.

## Reproducing / updating

The set was generated with throttled E-utilities calls (~3 requests/sec, sending
`tool` and `email` parameters per NCBI etiquette). To refresh it, re-run the same
ClinVar `esearch`/`esummary` + PubMed verification pipeline and update the
retrieval date above. Because ClinVar classifications change over time, always
re-verify labels against the live record before relying on them.

## Licensing / attribution

ClinVar and PubMed are public resources from the U.S. National Center for
Biotechnology Information (NCBI/NLM). Data is used here for non-commercial
educational purposes. See:
- ClinVar: https://www.ncbi.nlm.nih.gov/clinvar/
- PubMed: https://pubmed.ncbi.nlm.nih.gov/
- NCBI E-utilities usage policy: https://www.ncbi.nlm.nih.gov/books/NBK25497/

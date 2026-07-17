# 🧬 Gene Variant Explorer

🔗 **Live app:** [genescope.streamlit.app](https://genescope.streamlit.app/)

A student-built, education-first web app that makes genetics understandable
using real public scientific data. It's meant to help students, families, and
teachers explore how specific gene variants relate to health conditions —
without needing a background in genetics.

Built with [Streamlit](https://streamlit.io) and Python.

## What it does

- **Search a gene** (e.g. `BRCA1`, `TP53`) with fuzzy matching, so small typos
  still find the right gene.
- **Browse known variants** for that gene and the conditions they're
  associated with, pulled from a ClinVar-style dataset.
- **Filter results** by clinical significance (Pathogenic, Benign, VUS, etc.)
  and by condition.
- **Read plain-language summaries** — a one-sentence, human-readable summary
  of the relevant research, generated from PubMed abstracts via the NCBI
  E-utilities API (with a small hand-written glossary for well-known genes
  like BRCA1, TP53, and CFTR).
- **Generate shareable "cards"** for individual variants, and export search
  results as CSV.
- **Learn the vocabulary** — an About page explains variant notation
  (`c.`, `p.`, `fs`, `del`, `ins`, etc.) and clinical significance terms in
  plain English.
- **Share searches via URL** — filters are stored in query params so a
  specific search can be linked to someone else.
- **Collect feedback** from users through an in-app form that submits
  directly to a Google Form.

The app is organized into five pages: **Home**, **Tool**, **Impact**,
**About**, and **Feedback**. The Impact page shows real usage metrics (each
shows "Coming soon..." until it's actually being tracked); the developer's
internal notes on proving impact live in `impact_notes.md`, not in the app.

### Data sources
- [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) (NIH) — variants and
  associated conditions (a sample dataset is bundled as `clinvar_sample.csv`)
- [PubMed](https://pubmed.ncbi.nlm.nih.gov/) abstracts via NCBI E-utilities —
  research summaries

⚠️ **Educational use only. Not medical advice — please consult a doctor.**
No personal health information is collected. Feedback submissions (name,
role, message) are sent to a Google Form connected to the developer's
account (see `GOOGLE_FORM_ID` and the `GOOGLE_FORM_ENTRY_*` constants near
`render_feedback()` in `app.py`).

## Running it locally

```bash
pip install -r requirements.txt
streamlit run app.py
```

The app also includes a `.devcontainer` config, so it can be opened directly
in a GitHub Codespace with everything pre-installed.

## Deployment

This app deploys on **[Streamlit Community Cloud](https://streamlit.io/cloud)**,
straight from this GitHub repo's `main` branch (`app.py` as the entry point).
Streamlit Cloud auto-redeploys on every push to `main`.

- **Live URL:** [genescope.streamlit.app](https://genescope.streamlit.app/)
- **App name on Streamlit Cloud:** `genescope` (differs from the repo name,
  `gene-browser` — worth remembering when looking it up in the Streamlit
  Cloud dashboard)
- **Repo/branch:** `adheesh-c/gene-browser` @ `main`
- **Secrets configured on Streamlit Cloud:** `NCBI_EMAIL` (see
  `.streamlit/secrets.toml.example` for the format) — optional, but polite
  to NCBI's API etiquette rules.

### How to (re)deploy

1. Go to [share.streamlit.io](https://share.streamlit.io) and sign in with
   the GitHub account that owns this repo.
2. Click **New app** → pick repo `adheesh-c/gene-browser`, branch `main`,
   main file path `app.py`.
3. Under **Advanced settings → Secrets**, paste in the contents of
   `.streamlit/secrets.toml.example` with real values filled in.
4. Click **Deploy**. Streamlit Cloud builds from `requirements.txt` and gives
   you a public `*.streamlit.app` URL.
5. Copy that URL into the **Live app** line at the top of this README (and
   into the placeholder in the Deployment section above) so it's easy to
   find later.

Free-tier apps go to sleep after a period of inactivity and take a few
seconds to wake up on the next visit — that's expected, not a bug.

## Hosting options (for future toolkits)

The current app is Streamlit-specific, so Streamlit Community Cloud is the
right fit today. For reference, here's how it compares to other options if
future education toolkits in this series need something different:

| Option | Why it fits | Notes |
|---|---|---|
| **[Streamlit Community Cloud](https://streamlit.io/cloud)** (recommended) | Free, zero-config, built specifically for Streamlit apps. Connect your GitHub repo and it deploys automatically on every push. | Public URL like `your-app.streamlit.app`. Manage secrets (e.g. `NCBI_EMAIL`) in its dashboard. Sleeps after inactivity on the free tier, wakes on next visit. |
| **[Hugging Face Spaces](https://huggingface.co/spaces)** | Also free and supports Streamlit natively; good if you want your work discoverable alongside a portfolio of other projects. | Similar free-tier sleep behavior. |
| **[Render](https://render.com)** | General-purpose free web hosting; not Streamlit-specific, works for any app you deploy from GitHub. | Free tier spins down when idle and takes ~30–60s to wake up. |
| **[Railway](https://railway.app)** / **[Fly.io](https://fly.io)** | General-purpose hosting with small free/low-cost tiers. Good once you outgrow simple Streamlit apps (e.g. Flask/FastAPI/Node backends). | More setup than Streamlit Cloud; useful long-term skill to learn. |

**For this app specifically**, Streamlit Community Cloud is the easiest
starting point — sign in with GitHub, point it at this repo, and it's live.

**Longer-term**, since this is meant to be the first of a series of
education tools, it's worth picking one general-purpose platform (Render,
Railway, or Fly.io) to standardize on once apps stop being pure Streamlit —
that way each new toolkit doesn't need a brand-new hosting setup, and you
can put them all under one account/dashboard. Vercel or Netlify are strong
options too if any future project ends up being a JavaScript/React front end
instead of Python.

The Feedback page submits directly to a Google Form (see `app.py`), so
responses persist in the connected Google Sheet/Form regardless of the
host's filesystem — no local, ephemeral file involved.

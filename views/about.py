"""Page 4: About."""
import streamlit as st

from data_loader import AD_MIN_LEN, AD_THRESHOLD, MORFFY_CITATION, MORFFY_URL
from views.common import FULL_DATASETS, cards, full_dataset_getter, load_data, page_header

# ---- fill these in when available (empty string -> shown as "coming soon") ----
DOI = ""
PREPRINT_URL = ""
LAB_URL = ""
GITHUB_URL = "https://github.com/Leedog159/Poplar_TF_Viewer"
DATA_URL = ""          # Zenodo / supplementary tables

CONTACTS = [  # (name, role, email)
    ("Liam D. Kirkpatrick", "", "LDKirkpatrick@berkeley.edu"),
    ("Anthony E. Sarkiss", "", "Asarkiss@berkeley.edu"),
    ("Patrick M. Shih", "Corresponding Author", "pmshih@berkeley.edu"),
]

# Markdown: the co-first-author asterisks are escaped (\*) so they render literally.
CITATION_MD = (
    r"Kirkpatrick, L. D.\*, Sarkiss, A. E.\*, Waldburger, L., Kesselhaut, M. A., Zintel, M., Liu, S., "
    r"Keasling, J. D., Staller, M. V., & Shih, P. M. (2026). A genome-wide map of poplar transcription "
    r"factor activation domains uncovers pervasive turnover and tuning. "
    r"\[Manuscript in preparation / bioRxiv forthcoming\]. "
    + (f"DOI: [{DOI}](https://doi.org/{DOI})" if DOI else r"DOI: \[Pending submission\]")
)


def _link(label: str, url: str) -> str:
    return f"[{label}]({url})" if url else f"{label} — *coming soon*"


def render():
    # ---- citation (top) ------------------------------------------------------
    with st.container(border=True):
        page_header("📄", "How to Cite")
        st.info(CITATION_MD, icon=":material/format_quote:")
        st.caption(r"\* These authors contributed equally.")

    # ---- about the project --------------------------------------------------
    page_header("🌳", "About the Project")
    st.markdown(
        "Interactive, genome-wide access to **activation domains (ADs)** identified across "
        "**2,466 *Populus trichocarpa* transcription factors**, with comparative analyses against "
        "*Arabidopsis thaliana* orthologs. Browse per-residue activity predictions and empirical "
        "tile measurements for single TFs, compare ADs across orthogroups, and export AD calls in bulk."
    )
    data = load_data()
    sp = data.preds.species.value_counts()
    cards([
        ("Poplar TFs in viewer", f"{sp.get('Poplar', 0):,}", ""),
        ("Arabidopsis TFs in viewer", f"{sp.get('Arabidopsis', 0):,}", ""),
        ("Tiles designed / recovered", "87,798 / 75,533", ""),
        ("Orthogroups with MSA", f"{len(data._msa):,}", ""),
    ])

    # ---- methodology ---------------------------------------------------------
    page_header("🧪", "Methodology & Data Overview")
    with st.expander("Experimental screen: yeast sort-seq of 40-aa tiles", expanded=True):
        st.markdown(
            "- Every TF was tiled into **40-amino-acid fragments** (10-aa step), giving "
            "**87,798 designed tiles**, of which **75,533** were recovered in the screen.\n"
            "- Tiles were fused to a DNA-binding domain in yeast, sorted by reporter expression "
            "(sort-seq), and scored by their distribution across sort bins.\n"
            "- Empirical scores can be shown as grey 40-aa bars in the TF Viewer "
            "(*Show empirical yeast tiles*)."
        )
    with st.expander("Computational model: ADHunter deep ensemble", expanded=True):
        st.markdown(
            "- An **ensemble of 20 ResNet models (ADHunter)** was trained on the tile measurements.\n"
            "- Predictions are made in sliding windows and smoothed at **1-aa resolution**, giving a "
            "per-residue activity score that reduces measurement noise and resolves AD boundaries "
            "beyond the 10-aa tiling step.\n"
            "- The sky-blue curve in the viewers is this per-residue ensemble prediction."
        )
    with st.expander("In planta validation: *Nicotiana benthamiana*", expanded=True):
        st.markdown(
            "- A subset of predicted ADs was validated by transient expression in "
            "***N. benthamiana*** leaves, confirming that yeast-derived activity calls carry over "
            "to a plant cellular context."
        )

    # ---- external data ---------------------------------------------------------
    page_header("📚", "External Datasets & Attribution")
    with st.container(border=True):
        st.markdown(
            "**Arabidopsis empirical yeast tile activity** — the *Arabidopsis thaliana* sort-seq "
            "tiling data shown for Arabidopsis TFs is from:\n\n"
            f"> {MORFFY_CITATION} [{MORFFY_URL}]({MORFFY_URL})"
        )

    # ---- interpretation ------------------------------------------------------
    page_header("📏", "Score Interpretation & Thresholds")
    with st.container(border=True):
        st.markdown(
            f"""
- **Activity score** — sublibrary-standardized **z-score** (0 = sublibrary mean; units of SD).
- **Activation domain (AD) call** — a contiguous run with predicted activity **≥ {AD_THRESHOLD}**
  for at least **{AD_MIN_LEN} consecutive amino acids**. Coordinates are 1-based and inclusive.
- **Disorder** — Metapredict score; residues with score **> 0.5** are considered disordered.
- **Domains** — Pfam hits at gathering (GA) cutoffs; DNA-binding domains (DBDs) are shown in amber, other domains in slate.
"""
        )
    with st.expander("Notes on coordinates in this viewer"):
        st.markdown(
            "- Empirical tiles are placed on the protein by matching their sequence; the C-terminal tile "
            "covers the last 40 aa. Tiles screened from a different isoform keep their nominal "
            "position and are flagged *isoform mismatch* on hover.\n"
            "- The Orthogroup Viewer plots on multiple-sequence-alignment columns where a "
            "precomputed alignment exists; otherwise on residue positions.\n"
            "- Orthogroups are connected components of the poplar–Arabidopsis ortholog mapping, "
            "restricted to TF–TF pairs."
        )

    # ---- full dataset downloads --------------------------------------------
    page_header("💾", "Download Full Datasets",
                subtitle="Complete tables behind the viewer. Each file is generated on first click "
                         "(the full activity matrix and AD FASTA take ~30 s the very first time) and then cached.")
    n_ptr, n_ath = int(sp.get("Poplar", 0)), int(sp.get("Arabidopsis", 0))
    items = [
        ("matrix", "Predicted & Empirical Activity Matrix",
         f"All {n_ptr + n_ath:,} TFs ({n_ptr:,} Ptr, {n_ath:,} Ath): activator call, peak z-score, AD intervals "
         "and peak scores, active empirical tiles, Pfam domains, and full protein/CDS + AD/Pfam slice sequences."),
        ("orthology", "Orthology Mapping Table",
         f"All {len(data.mapping):,} poplar–Arabidopsis mapping rows with relationship type, TF flags, "
         "orthogroup ID, group composition (Ptr:Ath) and MSA availability."),
        ("pfam", "Pfam Domain Annotations",
         "Merged Pfam domain intervals for every TF (same merging as the viewers): domain, 1-based "
         "coordinates, DNA-binding flag, contributing Pfam models/accessions, best i-E-value."),
        ("ad_fasta", "All Identified ADs (FASTA)",
         "Every called AD peptide across both species; headers `>GeneID_ADn_start-end peak_score=X.XX`."),
    ]
    cols = st.columns(2)
    for i, (kind, title, desc) in enumerate(items):
        fname, mime = FULL_DATASETS[kind]
        with cols[i % 2].container(border=True):
            st.markdown(f"**{title}**")
            st.caption(desc)
            st.download_button(f"Download {fname}", data=full_dataset_getter(kind), file_name=fname,
                               mime=mime, on_click="ignore", icon=":material/download:",
                               key=f"full_{kind}", width="stretch")

    # ---- contact & links -----------------------------------------------------
    page_header("🔗", "Contact & Links")
    st.markdown("\n".join(
        f"- **{name}**" + (f" ({role})" if role else "") + f": [{email}](mailto:{email})"
        for name, role, email in CONTACTS))
    st.markdown(
        f"- **Lab website:** {_link('Shih Lab', LAB_URL)}\n"
        f"- **Code:** {_link('GitHub repository', GITHUB_URL)}\n"
        f"- **Raw data downloads:** {_link('Zenodo / supplementary tables', DATA_URL)}\n"
        f"- **Preprint:** {_link('bioRxiv', PREPRINT_URL)}"
    )

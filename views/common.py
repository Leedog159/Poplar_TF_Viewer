"""Shared state and cached data access for all pages."""
import hashlib
from pathlib import Path

import streamlit as st

import data_loader
from data_loader import DATA_DIR, FASTA_KINDS, TFData

ACTIVE_GENE = "selected_gene"  # gene ID shared across pages
DEFAULT_GENE = "Potri.001G258700"
SPECIES_LATIN = {"Poplar": "Populus trichocarpa", "Arabidopsis": "Arabidopsis thaliana"}

# ---- global look ---------------------------------------------------------------
_CSS = """
<style>
/* Clear Streamlit's fixed top toolbar so the atlas header is never covered.
   .block-container is the legacy class; data-testid is what Streamlit 1.4x+ renders. */
.block-container, [data-testid="stMainBlockContainer"] {
  padding-top: 4rem !important; padding-bottom: 2rem !important; }
.stDeployButton, [data-testid="stAppDeployButton"] { display: none !important; }

/* Theme-agnostic: text inherits Streamlit's theme color; muted text uses opacity; fills and
   borders are translucent neutrals, so the same CSS reads in both dark and light mode. */
.atlas-header { border-bottom: 1px solid rgba(128,128,128,.25); padding: 0 0 .8rem; margin-bottom: 1.1rem; }
.atlas-header .t { font-size: 1.6rem; font-weight: 700; letter-spacing: -.01em; line-height: 1.25; }
.atlas-header .s { opacity: .68; font-size: .93rem; margin-top: .25rem; line-height: 1.4; }
.page-h { display: flex; flex-wrap: wrap; align-items: baseline; gap: .6rem; margin: .2rem 0 .15rem; }
.page-h .t { font-size: 1.25rem; font-weight: 650; }
.page-h .b { font-size: .72rem; font-weight: 600; color: #0ea5e9; background: rgba(56,189,248,.12);
             border: 1px solid rgba(56,189,248,.3); border-radius: 999px; padding: .1rem .55rem; }
.page-sub { opacity: .68; font-size: .88rem; margin-bottom: .6rem; }
.cards { display: grid; grid-template-columns: repeat(auto-fit, minmax(150px, 1fr)); gap: .6rem;
         margin: .3rem 0 .9rem; }
.card { border: 1px solid rgba(128,128,128,.25); border-radius: 10px; padding: .5rem .8rem .55rem;
        background: rgba(128,128,128,.06); min-width: 0; }
.card .l { font-size: .7rem; text-transform: uppercase; letter-spacing: .05em; opacity: .68; }
.card .v { font-size: 1.05rem; font-weight: 600; overflow-wrap: anywhere; line-height: 1.35; }
.card .u { font-size: .78rem; opacity: .68; font-weight: 400; }
.tag { display: inline-block; border-radius: 999px; padding: .08rem .6rem; font-size: .85rem; font-weight: 600;
       vertical-align: middle; }
.tag.on { background: rgba(34,197,94,.16); color: #22c55e; border: 1px solid rgba(34,197,94,.4); }
.tag.off { background: rgba(128,128,128,.12); color: inherit; opacity: .8; border: 1px solid rgba(128,128,128,.35); }
.tag.sp { font-size: .75rem; padding: .02rem .5rem; margin-right: .35rem; }
.tag.ptr { background: rgba(56,189,248,.14); color: #38bdf8; border: 1px solid rgba(56,189,248,.45); }
.tag.ath { background: rgba(192,132,252,.14); color: #c084fc; border: 1px solid rgba(192,132,252,.45); }
.seq-h { font-size: .9rem; margin-bottom: .2rem; }
.seq-h .m { opacity: .68; }
/* per-TF header rows on the Orthogroup Viewer page */
.tf-row { display: flex; flex-wrap: wrap; align-items: center; gap: .35rem 1.1rem; padding: .45rem .1rem .1rem;
          border-top: 1px solid rgba(128,128,128,.22); }
.tf-row .id { font-weight: 650; font-size: 1rem; }
.tf-row .gene { opacity: .68; font-size: .82rem; }
.tf-row .k { opacity: .68; font-size: .74rem; text-transform: uppercase; letter-spacing: .05em; margin-right: .3rem; }
.tf-row .v { font-size: .9rem; font-weight: 600; }
.legend-row { display: flex; flex-wrap: wrap; gap: .3rem 1.1rem; font-size: .82rem; opacity: .85; margin-bottom: .2rem; }
.legend-row .sw { display: inline-block; vertical-align: middle; margin-right: .35rem; }
.atlas-footer { color: #888; font-size: .75rem; margin-top: 2.5rem; padding-top: .6rem;
                border-top: 1px solid rgba(128,128,128,.18); }
</style>
"""


def inject_css() -> None:
    st.markdown(_CSS, unsafe_allow_html=True)


def render_header() -> None:
    st.markdown(
        '<div class="atlas-header"><div class="t">🌳 Poplar TF Activation Domain Atlas</div>'
        '<div class="s">Genome-wide activation domain identification across <i>Populus trichocarpa</i> '
        'and <i>Arabidopsis thaliana</i> transcription factors.</div></div>',
        unsafe_allow_html=True)


FOOTER_CITATION = ("Citation: Kirkpatrick, L. D.*, Sarkiss, A. E.*, et al. (2026) — "
                   "Poplar TF Activation Domain Atlas.")


def render_footer() -> None:
    """One-line muted citation at the bottom of every page (called from app.py after the page)."""
    text = FOOTER_CITATION.replace("*", "&#42;")  # literal asterisks, never Markdown emphasis
    st.markdown(f'<div class="atlas-footer">{text}</div>', unsafe_allow_html=True)


def page_header(icon: str, title: str, badge: str | None = None, subtitle: str | None = None) -> None:
    b = f'<span class="b">{badge}</span>' if badge else ""
    st.markdown(f'<div class="page-h"><span class="t">{icon} {title}</span>{b}</div>'
                + (f'<div class="page-sub">{subtitle}</div>' if subtitle else ""), unsafe_allow_html=True)


def cards(items: list[tuple[str, str, str]]) -> None:
    """Row of metric cards; items = (label, value_html, unit/sub-text)."""
    html = "".join(f'<div class="card"><div class="l">{l}</div><div class="v">{v}'
                   + (f' <span class="u">{u}</span>' if u else "") + "</div></div>" for l, v, u in items)
    st.markdown(f'<div class="cards">{html}</div>', unsafe_allow_html=True)


def is_dark_theme() -> bool:
    """Active theme for chart colors. Streamlit infers it from the app background; it can lag
    one rerun behind a theme switch in the settings menu (CSS above is theme-agnostic)."""
    try:
        return (st.context.theme.type or "dark") != "light"
    except Exception:
        return True


def status_tag(rec) -> str:
    return ('<span class="tag on">Active</span>' if len(rec.ads)
            else '<span class="tag off" title="No AD detected">Inactive</span>')


def species_tag(rec) -> str:
    return ('<span class="tag sp ptr">Ptr</span>' if rec.species == "Poplar"
            else '<span class="tag sp ath">Ath</span>')


def gene_cards(rec) -> None:
    cards([
        (f"<i>{SPECIES_LATIN[rec.species]}</i>",
         f'{species_tag(rec)}{rec.protein_id}<br><span class="u">gene {rec.gene_id}</span>', ""),
        ("Activator status", status_tag(rec), ""),
        ("TF family", rec.family, ""),
        ("Protein length", f"{rec.length:,}", "aa"),
        ("Peak model score", f"{float(rec.preds.max()):.2f}", "z"),
        ("Identified ADs", str(len(rec.ads)), ""),
    ])


def _fingerprint() -> str:
    """Changes whenever data_loader.py or any input data file changes."""
    h = hashlib.sha1(Path(data_loader.__file__).read_bytes())
    for f in sorted(p for d in (DATA_DIR, DATA_DIR.parent) for p in d.iterdir() if p.is_file()):
        st_ = f.stat()
        h.update(f"{f.name}:{st_.st_size}:{st_.st_mtime_ns}".encode())
    return h.hexdigest()


@st.cache_resource(show_spinner="Loading data…", max_entries=1)
def _load_data(fingerprint: str) -> TFData:
    # cache_resource (not cache_data): one shared in-memory object; cache_data
    # would pickle/copy ~100 MB of arrays on every access.
    return TFData()


@st.cache_resource
def _cache_state() -> dict:
    return {"fingerprint": None}


def load_data() -> TFData:
    """Cached TFData, rebuilt automatically when the loader code or data files change.

    Without this, a running server keeps serving a TFData built from an older
    data_loader.py (stale methods/fields) after the code is edited.
    """
    fp = _fingerprint()
    state = _cache_state()
    if state["fingerprint"] != fp:
        if state["fingerprint"] is not None:
            st.cache_data.clear()  # derived tables were computed from the old object
        state["fingerprint"] = fp
    return _load_data(fp)


@st.cache_data(show_spinner=False, max_entries=256)
def orthogroup_summary(gene_id: str):
    data = load_data()
    return data.orthogroup_summary(data.orthogroup(gene_id))


@st.cache_data(show_spinner=False, max_entries=64)
def bulk_results(tokens: tuple[str, ...]):
    """(summary DataFrame, CSV bytes, {fasta kind: text}, not-found IDs) for raw tokens.

    The CSV is written here (cached) with full, untruncated sequence strings."""
    data = load_data()
    ids, missing = data.resolve(list(tokens))
    df = data.bulk_summary(ids)
    fastas = {k: data.fasta(ids, what=k) for k in FASTA_KINDS}
    return df, df.to_csv(index=False).encode("utf-8"), fastas, missing


@st.cache_data(show_spinner=False)
def family_counts(species: str | None, activator: bool | None = None) -> dict[str, int]:
    return load_data().family_counts(species, activator)


@st.cache_data(show_spinner=False)
def gene_ids(species: str | None, family: str | None, activator: bool | None = None) -> list[str]:
    return load_data().gene_ids(species, family, activator)


# ---- "Filter by Activity" (TF Viewer + Orthogroup Viewer) -----------------------
ACTIVITY_ALL = "All"
ACTIVITY_FILTERS = {ACTIVITY_ALL: None, "Activators only": True, "Non-activators only": False}


def activity_matches(choice: str, rec) -> bool:
    want = ACTIVITY_FILTERS.get(choice)
    return want is None or bool(len(rec.ads)) == want


# ---- CSV / FASTA export block (Bulk Lookup + Orthogroup Viewer) -----------------
FASTA_BUTTONS = [  # (kind, label, filename)
    ("ad_aa", "Download AD Peptides (FASTA)", "ad_peptides.fasta"),
    ("ad_cds", "Download AD CDS (FASTA)", "ad_cds.fasta"),
    ("protein", "Download Full TF Proteins (FASTA)", "tf_proteins.fasta"),
    ("cds", "Download Full TF CDS (FASTA)", "tf_cds.fasta"),
]


def export_section(df, csv_bytes: bytes, fastas: dict[str, str], file_prefix: str = "") -> None:
    """Summary CSV button + FASTA downloads for a bulk_results() payload.

    `file_prefix` is prepended to every download file name (contents are unaffected)."""
    st.download_button("📥 Download Summary CSV", csv_bytes, file_name=f"{file_prefix}tf_ad_summary.csv",
                       mime="text/csv", type="secondary",
                       help="All columns, including full untruncated sequence columns.")

    with st.expander("📥 FASTA Sequence Downloads", expanded=True):
        cols = st.columns(2)
        for i, (kind, label, fname) in enumerate(FASTA_BUTTONS):
            text_ = fastas[kind]
            n = text_.count(">")
            cols[i % 2].download_button(f"{label} · {n}", text_.encode(), file_name=f"{file_prefix}{fname}",
                                        mime="text/plain", disabled=not text_, width="stretch",
                                        key=f"dl_{kind}")
        n_ads, n_ad_cds = fastas["ad_aa"].count(">"), fastas["ad_cds"].count(">")
        n_cds_missing = int((~df.CDS_Available).sum())
        notes = []
        if not n_ads:
            notes.append("No ADs were called in the matched genes, so the AD exports are empty.")
        if n_cds_missing:
            notes.append(f"{n_cds_missing} gene(s) have no CDS on file and are omitted from the CDS exports"
                         + (f" ({n_ads - n_ad_cds} AD(s))." if n_ads > n_ad_cds else "."))
        notes.append("Headers: `>GeneID_AD1_start-end peak_score=X.XX` (peptides), "
                     "`>GeneID_AD1_CDS_start-end …` (CDS); AA coordinates, 1-based inclusive. "
                     "Full CDS include the stop codon where present.")
        st.caption(" ".join(notes))


# ---- widget-state persistence across pages -------------------------------------
# Streamlit discards a widget's session_state entry when its page is not rendered,
# so each page restores its widgets from a shadow key and saves them back.
def restore(key: str, default):
    if key not in st.session_state:
        st.session_state[key] = st.session_state.get(f"_keep_{key}", default)


def save(*keys: str) -> None:
    for k in keys:
        if k in st.session_state:
            st.session_state[f"_keep_{k}"] = st.session_state[k]


def kept(key: str, default=None):
    return st.session_state.get(f"_keep_{key}", default)


FULL_DATASETS = {  # kind -> (file name, mime)
    "matrix": ("tf_activation_domain_matrix.csv", "text/csv"),
    "orthology": ("poplar_arabidopsis_orthology.csv", "text/csv"),
    "pfam": ("pfam_domains_merged.csv", "text/csv"),
    "ad_fasta": ("all_activation_domains.fasta", "text/plain"),
}


@st.cache_data(show_spinner=False, max_entries=len(FULL_DATASETS), persist="disk")
def _full_dataset(kind: str, fingerprint: str) -> bytes:
    """Build one full-dataset export. Persisted to disk (keyed by the data fingerprint),
    so the ~30 s cold build of the all-TF tables happens once, not per session/restart."""
    data = load_data()
    if kind == "matrix":
        return data.bulk_summary(data.gene_ids()).to_csv(index=False).encode("utf-8")
    if kind == "orthology":
        return data.orthology_table().to_csv(index=False).encode("utf-8")
    if kind == "pfam":
        return data.pfam_table().to_csv(index=False).encode("utf-8")
    if kind == "ad_fasta":
        return data.fasta(data.gene_ids(), what="ad_aa").encode("utf-8")
    raise ValueError(kind)


def full_dataset_getter(kind: str):
    """Zero-argument callable for st.download_button(data=...): runs only on click."""
    return lambda: _full_dataset(kind, _fingerprint())


def active_gene(data: TFData) -> str:
    g = st.session_state.get(ACTIVE_GENE, DEFAULT_GENE)
    return g if data.has(g) else DEFAULT_GENE


def set_active_gene(gene_id: str) -> None:
    st.session_state[ACTIVE_GENE] = gene_id

"""Page 3: Bulk Lookup & Export."""
import streamlit as st

from data_loader import (AD_MIN_LEN, AD_THRESHOLD, BULK_SEQUENCE_COLS, BULK_SUMMARY_COLS,
                         parse_id_file, parse_id_text)
from views.common import bulk_results, cards, export_section, page_header, restore, save

EXAMPLE = "Potri.001G258700\nPotri.013G113100\nPotri.009G053900\nAT5G12870"

# Display labels only; the CSV keeps the raw column names.
_TXT = st.column_config.TextColumn
COLUMN_CONFIG = {
    "Gene_ID": _TXT("Gene ID", pinned=True),
    "Family": _TXT("Family"),
    "Has_AD": st.column_config.CheckboxColumn("Activator"),
    "Num_ADs": st.column_config.NumberColumn("# ADs", format="%d"),
    "Peak_Predicted_Score": st.column_config.NumberColumn("Peak score (z)", format="%.2f"),
    "Species": _TXT("Species"),
    "AD_Coordinates": _TXT("AD coordinates"),
    "AD_Peak_Scores": _TXT("AD peak scores"),
    "Pfam_Domains": _TXT("Pfam domains"),
    "Empirical_Active_Tiles": _TXT("Active tiles (≥ 1.0)"),
    "Protein_ID": _TXT("Protein ID"),
    "Length": st.column_config.NumberColumn("Length (aa)", format="%d"),
    "CDS_Available": st.column_config.CheckboxColumn("CDS"),
    "Full_Protein_Seq": _TXT("Protein sequence", width="medium"),
    "Full_CDS_Seq": _TXT("CDS sequence", width="medium"),
    "AD_Peptide_Seqs": _TXT("AD peptides", width="medium"),
    "AD_CDS_Seqs": _TXT("AD CDS", width="medium"),
    "Pfam_Peptide_Seqs": _TXT("Pfam peptides", width="medium"),
    "Pfam_CDS_Seqs": _TXT("Pfam CDS", width="medium"),
}


SUBMITTED = "bulk_submitted"  # tokens of the last Lookup (plain key: survives page switches)


def _load_example():
    st.session_state["bulk_text"] = EXAMPLE  # fills the box only; the user still presses Lookup


def _current_tokens(text: str, up) -> list[str]:
    tokens = parse_id_text(text)
    if up is not None:
        tokens += parse_id_file(up.getvalue(), up.name)
    return tokens


def render():
    page_header("📦", "Bulk Lookup & Export", "CSV · FASTA",
                "Poplar (Potri.…) and Arabidopsis (AT…G…) gene or transcript IDs, separated by any mix of "
                "newlines, commas, tabs, semicolons or spaces; case-insensitive (AT1G71930 and AT1G71930.1 "
                f"both match). ADs = predicted activity ≥ {AD_THRESHOLD} for ≥ {AD_MIN_LEN} aa.")

    restore("bulk_text", "")
    c1, c2 = st.columns([3, 2])
    with c1:
        text = st.text_area("Paste gene IDs", key="bulk_text", height=160, placeholder=EXAMPLE)
        b1, b2, _ = st.columns([1.3, 1.3, 3])
        lookup = b1.button("🔍 Lookup", type="secondary", width="stretch")
        b2.button("Load example", type="secondary", on_click=_load_example, width="stretch",
                  help="Fills the box with example IDs (replaces its contents); press Lookup to run.")
    with c2:
        up = st.file_uploader("…or upload a .txt / .csv", type=["txt", "csv"])
    save("bulk_text")

    tokens = _current_tokens(text, up)
    if lookup:
        if not tokens:
            st.session_state.pop(SUBMITTED, None)
            st.warning("Enter or upload gene IDs, then press Lookup.")
            return
        st.session_state[SUBMITTED] = tuple(tokens)
    submitted = st.session_state.get(SUBMITTED)
    if not submitted:
        st.info("Paste or upload gene IDs, then press **🔍 Lookup**.")
        return
    if tuple(tokens) != submitted:
        st.caption("ℹ️ The input has changed since the last lookup. Press **🔍 Lookup** to update the results.")
    tokens = list(submitted)

    with st.spinner("Looking up genes…"):
        df, csv_bytes, fastas, missing = bulk_results(tuple(tokens))

    if missing:
        st.warning(f"{len(missing)} ID(s) not found in the dataset: " + ", ".join(missing[:50])
                   + (" …" if len(missing) > 50 else ""))
    if df.empty:
        st.error("None of the input IDs matched the dataset.")
        return

    cards([
        ("Genes matched", f"{len(df):,}", ""),
        ("Activators", f"{int(df.Has_AD.sum()):,}", "≥ 1 AD"),
        ("Total ADs", f"{int(df.Num_ADs.sum()):,}", ""),
        ("Not found", f"{len(missing):,}", ""),
    ])

    restore("bulk_show_seqs", False)
    show_seqs = st.toggle("Show sequence columns", key="bulk_show_seqs",
                          help="Full protein/CDS and AD/Pfam slice sequences. Always included in the CSV.")
    save("bulk_show_seqs")
    st.dataframe(df, hide_index=True, width="stretch", column_config=COLUMN_CONFIG,
                 column_order=BULK_SUMMARY_COLS + (BULK_SEQUENCE_COLS if show_seqs else []))

    export_section(df, csv_bytes, fastas)

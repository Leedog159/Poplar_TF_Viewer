"""Page 1: TF Viewer."""
import numpy as np
import plotly.graph_objects as go
import streamlit as st
from plotly.subplots import make_subplots

from data_loader import AD_MIN_LEN, AD_THRESHOLD, MORFFY_URL
from plots import (BAR_TEXT, COLORS, DISORDER_CUTOFF, FONT, PLOTLY_CONFIG, color, region_labels, runs,
                   style_figure)
from views.common import (active_gene, family_counts, gene_cards, gene_ids, is_dark_theme, kept, load_data,
                          page_header, restore, save, set_active_gene)

ALL = "All Families"


def _select_gene(data) -> str:
    current = data.get(active_gene(data))
    restore("single_species", current.species)
    restore("single_family", ALL)
    # A gene selected on another page (e.g. "View in TF Viewer") takes over:
    # switch species if needed and drop a family filter that would hide it.
    if kept("single_gene") != current.gene_id:
        st.session_state["single_species"] = current.species
        if st.session_state["single_family"] not in (ALL, current.family):
            st.session_state["single_family"] = ALL
        st.session_state["single_gene"] = current.gene_id
    species = st.radio("Species", ["Poplar", "Arabidopsis"], horizontal=True, key="single_species")

    counts = family_counts(species)
    fam_options = [ALL] + list(counts)
    if st.session_state.get("single_family") not in fam_options:
        st.session_state["single_family"] = ALL
    total = sum(counts.values())
    c1, c2 = st.columns([1, 2])
    family = c1.selectbox("Filter by Family", fam_options, key="single_family",
                          format_func=lambda f: f"{ALL} ({total})" if f == ALL else f"{f} ({counts[f]})")

    ids = gene_ids(species, None if family == ALL else family)
    if st.session_state.get("single_gene") not in ids:
        st.session_state["single_gene"] = current.gene_id if current.gene_id in ids else ids[0]
    gene_id = c2.selectbox(f"Gene ID ({len(ids)} shown, type to search)", ids, key="single_gene")
    set_active_gene(gene_id)
    save("single_species", "single_family", "single_gene")
    return gene_id


def build_figure(rec, show_tiles: bool = False, dark: bool = True) -> go.Figure:
    """Activity line (row 1), AD/Pfam interval lanes directly beneath it (row 2), disorder (row 3)."""
    L = rec.length
    x = np.arange(1, L + 1)
    emp_label = f"Empirical Yeast Activity ({rec.empirical_source})"
    fig = make_subplots(rows=3, cols=1, shared_xaxes=True, vertical_spacing=0.04,
                        row_heights=[0.6, 0.18, 0.22])
    shapes = []  # assigned once at the end (per-call add_shape is slow)

    # ---- Row 1: continuous activity line (+ optional empirical tiles) -------
    if show_tiles and len(rec.tiles):
        tx, ty, tt = [], [], []
        for s, e, a, pl in zip(rec.tiles.start, rec.tiles.end, rec.tiles.activity, rec.tiles.placement):
            label = (f"Empirical tile {int(s)}–{int(e)}: <b>{a:.2f}</b>"
                     + (" (isoform mismatch)" if pl == "nominal" else ""))
            tx += [int(s) - 0.5, int(e) + 0.5, None]
            ty += [a, a, None]
            tt += [label, label, None]
        fig.add_trace(go.Scatter(x=tx, y=ty, mode="lines", name=emp_label, text=tt,
                                 hovertemplate="%{text}<extra></extra>",
                                 line=dict(color=COLORS["tile"], width=4)), row=1, col=1)
    regions = region_labels(rec)
    custom = np.stack([np.array(list(rec.seq), dtype=object),
                       np.where(regions != "", "<br>" + regions, "")], axis=1)
    fig.add_trace(go.Scatter(x=x, y=rec.preds, mode="lines", name="Predicted activity (ADHunter)",
                             line=dict(color=color("pred", dark), width=2.5), customdata=custom,
                             hovertemplate="Residue %{x} <b>%{customdata[0]}</b> · predicted <b>%{y:.2f}</b>"
                                           "%{customdata[1]}<extra></extra>"), row=1, col=1)
    shapes.append(dict(type="line", xref="x domain", yref="y", x0=0, x1=1, y0=AD_THRESHOLD, y1=AD_THRESHOLD,
                       line=dict(color=COLORS["threshold"], dash="dash", width=1.5), opacity=0.9))
    fig.add_annotation(x=1, xref="x domain", y=AD_THRESHOLD, yref="y", text="AD threshold (1.0)",
                       showarrow=False, xanchor="right", yanchor="bottom",
                       font=dict(size=10, color=COLORS["threshold"]))

    # ---- Row 2: dedicated interval lanes for ADs and Pfam domains -----------
    def add_bars(starts, ends, lane, kinds, names, details):
        if not len(starts):
            return
        starts, ends = np.asarray(starts, int), np.asarray(ends, int)
        hover = [f"<b>{n}</b>: residues {s}–{e} ({e - s + 1} aa){d}"
                 for n, s, e, d in zip(names, starts, ends, details)]
        fig.add_trace(go.Bar(
            y=[lane] * len(starts), x=(ends - starts + 1).tolist(), base=(starts - 0.5).tolist(),
            orientation="h", width=0.66, showlegend=False, hovertext=hover, hoverinfo="text",
            marker=dict(color=[COLORS[k] for k in kinds], cornerradius=4, line_width=0),
            text=names, textposition="inside", insidetextanchor="middle", constraintext="none",
            textfont=dict(color=[BAR_TEXT[k] for k in kinds], size=11, family=FONT)), row=2, col=1)

    add_bars(rec.ads.Start, rec.ads.End, "AD", ["ad"] * len(rec.ads),
             [a.replace("AD", "AD ") for a in rec.ads.AD],
             [f"<br>peak {float(p):.2f}" for p in rec.ads["Peak Score"]])
    add_bars(rec.pfam.start, rec.pfam.end, "Pfam", ["dbd" if d else "pfam" for d in rec.pfam.is_dbd],
             rec.pfam.pfam_name.tolist(),
             [("<br>DNA-binding" if d else "") + f"<br>{m}" for d, m in zip(rec.pfam.is_dbd, rec.pfam.members)])
    fig.update_yaxes(categoryorder="array", categoryarray=["Pfam", "AD"], type="category",
                     range=[-0.6, 1.6], showgrid=False, row=2, col=1)

    # ---- Row 3: disorder (slate bands where >= 0.5) -------------------------
    for s, e in runs(rec.disorder >= DISORDER_CUTOFF):
        shapes.append(dict(type="rect", xref="x3", yref="y3 domain", x0=s - 0.5, x1=e + 0.5, y0=0, y1=1,
                           fillcolor=COLORS["disorder_fill"], line_width=0, layer="below"))
    shapes.append(dict(type="line", xref="x3 domain", yref="y3", x0=0, x1=1, y0=DISORDER_CUTOFF,
                       y1=DISORDER_CUTOFF, line=dict(color=color("axis", dark), dash="dot", width=1)))
    fig.add_trace(go.Scatter(x=x, y=rec.disorder, mode="lines", name="Disorder (Metapredict)",
                             line=dict(color=color("disorder_line", dark), width=1.4), showlegend=False,
                             hovertemplate="Residue %{x} · disorder <b>%{y:.2f}</b><extra></extra>"),
                  row=3, col=1)
    fig.add_trace(go.Scatter(x=[None], y=[None], mode="markers", name="Disordered (≥ 0.5)",
                             marker=dict(symbol="square", size=11, color=COLORS["disorder_fill"],
                                         line=dict(color=color("axis", dark), width=1))), row=3, col=1)

    style_figure(fig, height=540, dark=dark)
    fig.update_layout(shapes=shapes, barmode="overlay")
    fig.update_xaxes(range=[0.5, L + 0.5], showspikes=False)
    fig.update_xaxes(title_text="Amino acid position", row=3, col=1)
    fig.update_yaxes(title_text="Activity (z)", row=1, col=1)
    fig.update_yaxes(title_text="Disorder", range=[0, 1], dtick=0.5, row=3, col=1)
    return fig


def _seq_code(title: str, seq: str | None):
    st.markdown(f'<div class="seq-h">{title}</div>', unsafe_allow_html=True)
    if seq is None:
        st.caption("Not available")
    else:
        st.code(seq, language="text", wrap_lines=True)


def _interval_card(rec, title: str, start, end, extra: str = ""):
    s, e = int(start), int(end)
    cds_rng = f"CDS: {(s - 1) * 3 + 1}–{e * 3} bp" if rec.cds is not None else "CDS: n/a"
    with st.container(border=True):
        st.markdown(f'<div class="seq-h"><b>{title}</b>: Residues {s}–{e} ({e - s + 1} aa) '
                    f'<span class="m">| {cds_rng}{extra}</span></div>', unsafe_allow_html=True)
        c1, c2 = st.columns(2)
        with c1:
            _seq_code("Peptide Sequence (aa)", rec.subseq(s, e))
        with c2:
            _seq_code("Coding Sequence (CDS)", rec.cds_slice(s, e))


def _sequences_section(rec):
    with st.expander("🧬 Sequences & Annotated Domains", expanded=True):
        t_full, t_ad, t_pf = st.tabs([":material/description: Full sequences",
                                      f":material/bolt: Activation domains ({len(rec.ads)})",
                                      f":material/category: Pfam domains ({len(rec.pfam)})"])
        with t_full:
            c1, c2 = st.columns(2)
            with c1:
                _seq_code(f"Protein · {rec.length:,} aa", rec.seq)
            with c2:
                _seq_code(f"CDS · {len(rec.cds):,} nt" if rec.cds else "CDS", rec.cds)
        with t_ad:
            if rec.ads.empty:
                st.caption(f"No ADs called (predicted ≥ {AD_THRESHOLD} for ≥ {AD_MIN_LEN} aa).")
            for name, s, e, pk in zip(rec.ads.AD, rec.ads.Start, rec.ads.End, rec.ads["Peak Score"]):
                _interval_card(rec, name.replace("AD", "AD "), s, e, f" · peak {float(pk):.2f}")
        with t_pf:
            if rec.pfam.empty:
                st.caption("No Pfam domains.")
            for p in rec.pfam.itertuples(index=False):
                extra = (" · DNA-binding" if p.is_dbd else "") + f" · {p.members}"
                _interval_card(rec, p.pfam_name, p.start, p.end, extra)


def render():
    data = load_data()
    page_header("🔬", "TF Viewer", "Per-residue activity",
                f"ADHunter predictions, empirical yeast tiles, disorder and domains. ADs are called where "
                f"predicted activity ≥ {AD_THRESHOLD} for ≥ {AD_MIN_LEN} consecutive residues.")
    rec = data.get(_select_gene(data))
    gene_cards(rec)

    # off by default for both species, so the model prediction leads; remembered per species
    tiles_key = f"show_tiles_{rec.species}"
    restore(tiles_key, False)
    show_tiles = st.checkbox("Show empirical yeast tiles", key=tiles_key,
                             help=f"Empirical Yeast Activity ({rec.empirical_source})")
    save(tiles_key)

    st.plotly_chart(build_figure(rec, show_tiles, dark=is_dark_theme()), width="stretch", theme=None,
                    config=PLOTLY_CONFIG)
    if show_tiles and rec.species == "Arabidopsis":
        st.caption(f"Arabidopsis empirical sort-seq tile data sourced from "
                   f"[Morffy et al. (Nature 2024)]({MORFFY_URL}).")
    st.caption("Interval lanes below the activity line: cyan = AD, amber = Pfam DNA-binding domain, slate = other "
               "Pfam domain (overlapping hits of the same domain merged). Grey bands = disordered (≥ 0.5)"
               + ("; grey 40-aa bars = empirical yeast tiles." if show_tiles else "."))

    if not rec.ads.empty:
        st.dataframe(rec.ads[["AD", "Start", "End", "Length", "Peak Score", "Mean Score"]],
                     hide_index=True, width="content", column_config={
                         "Peak Score": st.column_config.NumberColumn(format="%.2f"),
                         "Mean Score": st.column_config.NumberColumn(format="%.2f")})

    _sequences_section(rec)
    st.caption("Tip: open **Orthogroup Viewer** in the sidebar to compare this gene's orthogroup.")

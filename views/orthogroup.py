"""Page 2: Orthogroup Viewer."""
import numpy as np
import plotly.graph_objects as go
import streamlit as st
from plotly.subplots import make_subplots

from data_loader import AD_THRESHOLD
from plots import BAR_TEXT, COLORS, FONT, PLOTLY_CONFIG, SPECIES_TAG, color, region_labels, style_figure
from views.common import (ACTIVITY_ALL, ACTIVITY_FILTERS, active_gene, activity_matches, bulk_results, cards,
                          export_section, family_counts, gene_ids, is_dark_theme, kept, load_data,
                          orthogroup_summary, page_header, restore, save, set_active_gene, species_tag,
                          status_tag)

TRACK_PX = 190
ALL = "All Families"
MSA = "Aligned Coordinates (MSA)"
RESIDUE = "Residue Coordinates (Ungapped)"


def _select_gene(data) -> str:
    current = data.get(active_gene(data))
    restore("og_family", ALL)
    restore("og_activity", ACTIVITY_ALL)
    if kept("og_gene") != current.gene_id:  # gene chosen on another page takes over
        if st.session_state["og_family"] not in (ALL, current.family):
            st.session_state["og_family"] = ALL
        if not activity_matches(st.session_state["og_activity"], current):
            st.session_state["og_activity"] = ACTIVITY_ALL
        st.session_state["og_gene"] = current.gene_id

    # Family | Activity | Gene ID on one row; activity is read first (it sets the family options).
    c_fam, c_act, c_gene = st.columns([1, 1, 2])
    activity = c_act.selectbox("Filter by Activity", list(ACTIVITY_FILTERS), key="og_activity")
    activator = ACTIVITY_FILTERS[activity]

    counts = family_counts(None, activator)  # both species combined
    fam_options = [ALL] + list(counts)
    if st.session_state["og_family"] not in fam_options:
        st.session_state["og_family"] = ALL
    total = sum(counts.values())
    family = c_fam.selectbox("Filter by Family", fam_options, key="og_family",
                             format_func=lambda f: f"{ALL} ({total})" if f == ALL else f"{f} ({counts[f]})")
    ids = gene_ids(None, None if family == ALL else family, activator)
    if st.session_state.get("og_gene") not in ids:
        st.session_state["og_gene"] = current.gene_id if current.gene_id in ids else ids[0]
    gene_id = c_gene.selectbox(f"Poplar or Arabidopsis gene ID ({len(ids)} shown, type to search)", ids,
                               key="og_gene")
    set_active_gene(gene_id)
    save("og_family", "og_activity", "og_gene")
    return gene_id


def _open_in_single(gene_id: str) -> None:
    set_active_gene(gene_id)
    from views.nav import SINGLE  # local import: nav imports this module
    st.switch_page(SINGLE)


def _shared_scale(recs, og, use_msa):
    """Common y-range and x-extent so the separate per-TF tracks stay comparable and aligned."""
    ymin = min(float(r.preds.min()) for r in recs)
    ymax = max(float(r.preds.max()) for r in recs)
    yrange = [min(ymin, 0) - 0.3, max(ymax, AD_THRESHOLD) + 0.3]
    xmax = og.aln_len if use_msa else max(r.length for r in recs)
    return yrange, xmax


def build_track(rec, og, use_msa: bool, yrange, xmax: int, first: bool, last: bool,
                dark: bool = True) -> go.Figure:
    """One TF: activity line (row 1) with AD / Pfam interval lanes directly beneath it (row 2)."""
    fig = make_subplots(rows=2, cols=1, shared_xaxes=True, vertical_spacing=0.06, row_heights=[0.72, 0.28])
    resid = np.arange(1, rec.length + 1)
    cols = og.columns(rec.gene_id) if use_msa else resid
    to_x = (lambda p, c=cols: int(c[int(p) - 1]))  # residue (1-based) -> x coordinate
    regions = region_labels(rec)
    regions = np.where(regions != "", "<br><b>" + regions + "</b>", "")
    who = f"<b>{SPECIES_TAG[rec.species]} {rec.protein_id}</b>"

    # row 1: predicted activity; on the MSA scale gap columns are NaN, so the line breaks at indels
    if use_msa:
        x = np.arange(1, xmax + 1)
        y = np.full(xmax, np.nan)
        y[cols - 1] = rec.preds
        custom = np.full((xmax, 3), "", dtype=object)
        custom[cols - 1, 0] = resid.astype(str)
        custom[cols - 1, 1] = list(rec.seq)
        custom[cols - 1, 2] = regions
        hover = (who + "<br>column %{x} · residue %{customdata[0]} <b>%{customdata[1]}</b> · "
                 "predicted <b>%{y:.2f}</b>%{customdata[2]}<extra></extra>")
    else:
        x, y = resid, rec.preds
        custom = np.stack([np.array(list(rec.seq), dtype=object), regions], axis=1)
        hover = (who + "<br>residue %{x} <b>%{customdata[0]}</b> · predicted <b>%{y:.2f}</b>"
                 "%{customdata[1]}<extra></extra>")
    fig.add_trace(go.Scatter(x=x, y=y, mode="lines", line=dict(color=color("pred", dark), width=2.5),
                             customdata=custom, hovertemplate=hover, connectgaps=False,
                             name="Predicted activity (ADHunter)"), row=1, col=1)
    shapes = [dict(type="line", xref="x domain", yref="y", x0=0, x1=1, y0=AD_THRESHOLD, y1=AD_THRESHOLD,
                   opacity=0.9, line=dict(color=COLORS["threshold"], dash="dash", width=1.5))]
    if first:
        fig.add_annotation(x=1, xref="x domain", y=AD_THRESHOLD, yref="y", text="AD threshold (1.0)",
                           showarrow=False, xanchor="right", yanchor="bottom",
                           font=dict(size=10, color=COLORS["threshold"]))

    # row 2: interval lanes (bars span the mapped columns; in MSA mode they include internal gaps)
    def lane(starts, ends, name_lane, kinds, names, details):
        if not len(starts):
            return
        x0 = np.array([to_x(s) for s in starts]); x1 = np.array([to_x(e) for e in ends])
        hov = [f"{who}<br><b>{n}</b>: residues {int(s)}–{int(e)}{d}"
               for n, s, e, d in zip(names, starts, ends, details)]
        fig.add_trace(go.Bar(
            y=[name_lane] * len(x0), x=(x1 - x0 + 1).tolist(), base=(x0 - 0.5).tolist(), orientation="h",
            width=0.7, showlegend=False, hovertext=hov, hoverinfo="text",
            marker=dict(color=[COLORS[k] for k in kinds], cornerradius=3, line_width=0),
            text=names, textposition="inside", insidetextanchor="middle", constraintext="none",
            textfont=dict(color=[BAR_TEXT[k] for k in kinds], size=10, family=FONT)), row=2, col=1)

    lane(rec.ads.Start, rec.ads.End, "AD", ["ad"] * len(rec.ads),
         [a.replace("AD", "AD ") for a in rec.ads.AD], [f" · peak {float(p):.2f}" for p in rec.ads["Peak Score"]])
    lane(rec.pfam.start, rec.pfam.end, "Pfam", ["dbd" if d else "pfam" for d in rec.pfam.is_dbd],
         rec.pfam.pfam_name.tolist(), [" · DNA-binding" if d else "" for d in rec.pfam.is_dbd])

    style_figure(fig, height=TRACK_PX + (28 if last else 0), dark=dark)
    fig.update_layout(shapes=shapes, showlegend=False, barmode="overlay",
                      margin=dict(l=40, r=40, t=6, b=34 if last else 6))
    fig.update_xaxes(range=[0.5, xmax + 0.5])
    if last:
        fig.update_xaxes(title_text="Alignment column (MSA)" if use_msa else "Residue position (ungapped)",
                         row=2, col=1)
    fig.update_yaxes(range=yrange, title_text="z", title_standoff=2, row=1, col=1)
    fig.update_yaxes(type="category", categoryorder="array", categoryarray=["Pfam", "AD"],
                     range=[-0.6, 1.6], showgrid=False, tickfont=dict(size=10), row=2, col=1)
    return fig


def build_tracks(data, og, use_msa: bool = True, dark: bool = True) -> list[tuple]:
    """[(TFRecord, figure)] for each orthogroup member, on a shared scale."""
    use_msa = use_msa and og.aligned
    recs = [data.get(g) for g in og.members]
    yrange, xmax = _shared_scale(recs, og, use_msa)
    n = len(recs)
    return [(r, build_track(r, og, use_msa, yrange, xmax, i == 0, i == n - 1, dark)) for i, r in enumerate(recs)]


def _legend(dark: bool) -> str:
    sw = '<span class="sw" style="{}"></span>'
    items = [
        (f"width:18px;height:3px;background:{color('pred', dark)}", "Predicted activity (ADHunter)"),
        (f"width:18px;border-top:2px dashed {COLORS['threshold']}", "AD threshold (1.0)"),
        (f"width:12px;height:12px;border-radius:3px;background:{COLORS['ad']}", "Called AD"),
        (f"width:12px;height:12px;border-radius:3px;background:{COLORS['dbd']}", "Pfam DBD"),
        (f"width:12px;height:12px;border-radius:3px;background:{COLORS['pfam']}", "Other Pfam"),
    ]
    return '<div class="legend-row">' + "".join(f"<span>{sw.format(s)}{t}</span>" for s, t in items) + "</div>"


def _tf_header(rec, is_query: bool) -> str:
    q = ' <span class="gene">(query)</span>' if is_query else ""
    return (f'<div class="tf-row"><span>{species_tag(rec)}<span class="id">{rec.protein_id}</span></span>'
            f'<span>{status_tag(rec)}</span>'
            f'<span class="gene">gene {rec.gene_id}{q}</span>'
            f'<span><span class="k">Length</span><span class="v">{rec.length:,} aa</span></span>'
            f'<span><span class="k">ADs</span><span class="v">{len(rec.ads)}</span>'
            f' <span class="gene">(peak {float(rec.preds.max()):.2f} z)</span></span></div>')


def render():
    data = load_data()
    dark = is_dark_theme()
    page_header("🧬", "Orthogroup Viewer", "Orthologs & paralogs",
                "All TF orthologs and paralogs connected to the query in the poplar–Arabidopsis "
                "ortholog mapping, on a shared activity scale.")
    gene_id = _select_gene(data)
    og = data.orthogroup(gene_id)
    recs = [data.get(g) for g in og.members]

    families = sorted({r.family for r in recs})
    cards([
        ("TF family", " / ".join(families), ""),
        ("Poplar (Ptr) genes", str(len(og.ptr)), ""),
        ("Arabidopsis (Ath) genes", str(len(og.ath)), ""),
        ("Copy ratio", og.ratio, "Ptr:Ath"),
    ])
    if len(og.members) == 1:
        st.info("No TF orthologs found for this gene in mapping.pkl.")

    restore("og_coords", MSA)
    coords = st.radio("Coordinate space", [MSA, RESIDUE], horizontal=True, key="og_coords",
                      help="MSA: homologous positions line up vertically (curves break at gaps). "
                           "Ungapped: each protein on its own 1…length axis.")
    save("og_coords")
    use_msa = coords == MSA and og.aligned
    if coords == MSA and not og.aligned and len(og.members) > 1:
        st.warning("No precomputed alignment for this orthogroup, so tracks are shown on residue "
                   "coordinates. Compare using the domain markers.")

    st.markdown(_legend(dark), unsafe_allow_html=True)
    for rec, fig in build_tracks(data, og, use_msa, dark):
        c1, c2 = st.columns([5, 1.4], vertical_alignment="center")
        c1.markdown(_tf_header(rec, rec.gene_id == og.query), unsafe_allow_html=True)
        if c2.button("🔍 View in TF Viewer", key=f"open_{rec.gene_id}", width="stretch"):
            _open_in_single(rec.gene_id)
        st.plotly_chart(fig, width="stretch", theme=None, config=PLOTLY_CONFIG, key=f"track_{rec.gene_id}")
    st.caption("All tracks share the same y-axis (z-score) and x-extent. Interval lanes under each curve: "
               "cyan = called AD, amber = Pfam DNA-binding domain, slate = other Pfam domain. "
               + ("Breaks in a curve are alignment gaps." if use_msa else ""))

    page_header("📋", "Orthogroup summary")
    st.dataframe(orthogroup_summary(gene_id), hide_index=True, width="stretch", column_config={
        "Max pred": st.column_config.NumberColumn("Peak score (z)", format="%.2f"),
        "# ADs": st.column_config.NumberColumn("# ADs"),
    })

    # Same cached pipeline as Bulk Lookup, so the files match a Bulk Lookup of these members exactly.
    df, csv_bytes, fastas, _ = bulk_results(tuple(og.members))
    export_section(df, csv_bytes, fastas, file_prefix=f"orthogroup_{og.query}_")

"""Shared Plotly styling and per-residue helpers for the viewer's charts."""
import numpy as np

FONT = "Inter, 'Helvetica Neue', Helvetica, Arial, sans-serif"

# Data colors (both themes). Domain-bar trio (ad / dbd / pfam) validated against a dark surface:
# all >= 3:1 contrast; CVD and normal-vision separation pass.
COLORS = {
    "pred": "#38bdf8",                     # predicted activity (sky blue; darker in light mode)
    "tile": "rgba(148,163,184,0.45)",      # empirical yeast tiles: muted slate, recessive
    "threshold": "#f43f5e",                # AD threshold line (rose)
    "disorder_line": "#a5b4fc",            # disorder profile (soft periwinkle)
    "disorder_fill": "rgba(148,163,184,0.2)",  # disordered regions (>= 0.5)
    "ad": "#06b6d4",                       # activation domains (cyan)
    "dbd": "#f59e0b",                      # Pfam DNA-binding domains (amber)
    "pfam": "#64748b",                     # other Pfam domains (slate)
    "plot_bg": "rgba(0,0,0,0)",
}
# Theme-dependent chart chrome (+ a darker prediction line so it keeps contrast on white).
CHROME = {
    True: {"ink": "#e6e6e6", "muted": "#a3a8b4", "grid": "rgba(255,255,255,0.08)",
           "axis": "rgba(255,255,255,0.22)", "hover_bg": "#262730", "hover_ink": "#fafafa",
           "template": "plotly_dark", "disorder_line": "#a5b4fc"},
    False: {"ink": "#31333f", "muted": "#6b7280", "grid": "rgba(0,0,0,0.07)",
            "axis": "rgba(0,0,0,0.25)", "hover_bg": "#ffffff", "hover_ink": "#31333f",
            "template": "plotly_white", "pred": "#0284c7", "disorder_line": "#6366f1"},
}


def color(key: str, dark: bool = True) -> str:
    return CHROME[dark].get(key, COLORS.get(key))
# label text on bars: dark ink on bright cyan/amber, light ink on slate
BAR_TEXT = {"ad": "#0f172a", "dbd": "#0f172a", "pfam": "#f8fafc"}
SPECIES_TAG = {"Poplar": "Ptr", "Arabidopsis": "Ath"}
DISORDER_CUTOFF = 0.5

# Toolbar: keep zoom, pan, reset/home and image export; drop select/lasso and other clutter.
PLOTLY_CONFIG = {
    "displaylogo": False,
    "displayModeBar": True,
    "modeBarButtonsToRemove": ["select2d", "lasso2d", "zoomIn2d", "zoomOut2d", "autoScale2d",
                               "toggleSpikelines", "hoverClosestCartesian", "hoverCompareCartesian"],
    "toImageButtonOptions": {"format": "svg", "filename": "tf_activation_domains"},
}


def pfam_color(is_dbd) -> str:
    return COLORS["dbd"] if is_dbd else COLORS["pfam"]


def style_figure(fig, height: int, dark: bool = True):
    """Transparent canvas, faint grid, sans-serif type; plotly_dark or plotly_white chrome."""
    ch = CHROME[dark]
    fig.update_layout(
        template=ch["template"], height=height, font=dict(family=FONT, size=12, color=ch["ink"]),
        margin=dict(l=40, r=40, t=30, b=30), paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor=COLORS["plot_bg"],
        hovermode="closest", hoverlabel=dict(font=dict(family=FONT, size=12, color=ch["hover_ink"]),
                                             bgcolor=ch["hover_bg"], bordercolor=ch["axis"]),
        legend=dict(orientation="h", yanchor="bottom", y=1.0, xanchor="right", x=1,
                    font=dict(size=11), bgcolor="rgba(0,0,0,0)"),
    )
    axis = dict(showline=True, linecolor=ch["axis"], linewidth=1, mirror=True, gridcolor=ch["grid"],
                zeroline=False, ticks="outside", tickcolor=ch["axis"], ticklen=4,
                title_font=dict(size=12, color=ch["muted"]), tickfont=dict(size=11))
    fig.update_xaxes(**axis)
    fig.update_yaxes(**axis)
    return fig


def runs(mask: np.ndarray) -> list[tuple[int, int]]:
    """1-based inclusive (start, end) of each True run."""
    m = np.concatenate([[False], np.asarray(mask, bool), [False]])
    e = np.flatnonzero(np.diff(m.astype(np.int8)))
    return [(int(a) + 1, int(b)) for a, b in zip(e[::2], e[1::2])]


def region_labels(rec) -> np.ndarray:
    """Per-residue text naming the AD / Pfam features covering it (for unified hover)."""
    lab = [[] for _ in range(rec.length)]
    for name, s, e in zip(rec.ads.AD, rec.ads.Start, rec.ads.End):
        for j in range(int(s) - 1, int(e)):
            lab[j].append(name.replace("AD", "AD "))
    for name, s, e in zip(rec.pfam.pfam_name, rec.pfam.start, rec.pfam.end):
        for j in range(int(s) - 1, int(e)):
            lab[j].append(name)
    return np.array([" · ".join(x) if x else "" for x in lab], dtype=object)


def tile_mean(rec) -> np.ndarray:
    """Per-residue mean of the (sequence-placed) empirical tiles covering each residue."""
    tot = np.zeros(rec.length)
    n = np.zeros(rec.length)
    for s, e, a in zip(rec.tiles.start, rec.tiles.end, rec.tiles.activity):
        tot[int(s) - 1:int(e)] += a
        n[int(s) - 1:int(e)] += 1
    with np.errstate(invalid="ignore"):
        return np.where(n > 0, tot / n, np.nan)

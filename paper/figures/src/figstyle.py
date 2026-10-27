"""Shared style and geometry helpers for the figures of the S4 paper.

Every figure script inserts this directory into sys.path and imports this
module, so all figures use one palette, one font and one set of line widths.
"""
import itertools
import os

import matplotlib

matplotlib.use("pdf")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
FIGDIR = os.path.dirname(HERE)
DATADIR = os.path.join(FIGDIR, "data")

# The paper's constants (Section 3).
SBAR = 41821 / 50000          # s-bar
SBAR2 = SBAR * SBAR
OMEGA = 1.0 - SBAR2           # omega = 1 - s-bar^2

# Okabe-Ito palette (colour-blind safe).  ORANGE is the main highlight,
# BLUE the secondary highlight, GREEN a third series where one is needed.
ORANGE = "#D55E00"
BLUE = "#0072B2"
GREEN = "#009E73"
INK = "#1A1A1A"               # points and text
MUTED = "#555555"             # secondary text
CONSTR = "#9A9A9A"            # thin construction lines
LIGHT = "#C8C8C8"             # very light lines (full circle, grid)
RING_FILL = "#DCE3EB"         # annulus fill: light neutral blue-grey
RING_EDGE = "#8E9DAD"         # annulus boundary circles
ORANGE_TINT = "#F6D5BF"
BLUE_TINT = "#C6DCEC"
GREEN_TINT = "#C4E6DA"

FULL_WIDTH = 6.0              # inches

LW_THIN = 0.5
LW = 0.9
LW_THICK = 2.4


def setup():
    plt.rcParams.update({
        "text.usetex": True,
        "text.latex.preamble": r"\usepackage{amsmath}\usepackage{amssymb}",
        "font.family": "serif",
        "font.serif": ["Computer Modern Roman"],
        "font.size": 9,
        "axes.labelsize": 9,
        "axes.titlesize": 9,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 8,
        "axes.linewidth": 0.5,
        "xtick.major.width": 0.5,
        "ytick.major.width": 0.5,
        "xtick.minor.width": 0.4,
        "ytick.minor.width": 0.4,
        "xtick.major.size": 3,
        "ytick.major.size": 3,
        "xtick.minor.size": 1.8,
        "ytick.minor.size": 1.8,
        "axes.edgecolor": "#444444",
        "xtick.color": "#444444",
        "ytick.color": "#444444",
        "axes.labelcolor": INK,
        "text.color": INK,
        "lines.linewidth": LW,
        "lines.solid_capstyle": "round",
        "lines.dash_capstyle": "butt",
        "legend.frameon": False,
        "legend.handlelength": 2.2,
        "hatch.linewidth": 0.5,
        "pdf.fonttype": 42,
        "savefig.dpi": 300,
    })


def geometry_axes(ax):
    """Equal aspect, no frame, no ticks."""
    ax.set_aspect("equal")
    ax.set_axis_off()


def panel_label(ax, label, x=0.0, y=1.0, **kw):
    ax.text(x, y, label, transform=ax.transAxes, ha="left", va="top",
            fontsize=9, **kw)


def save(fig, name, tight=True):
    path = os.path.join(FIGDIR, name + ".pdf")
    if tight:
        fig.savefig(path, bbox_inches="tight", pad_inches=0.02)
    else:
        fig.savefig(path)
    plt.close(fig)
    return path


def draw_annulus(ax, inner=SBAR, outer=1.0, center=(0.0, 0.0), zorder=0,
                 fill=RING_FILL, edge=RING_EDGE, lw=LW_THIN):
    """Shade the annulus inner <= |x - center| <= outer."""
    t = np.linspace(0.0, 2.0 * np.pi, 721)
    cx, cy = center
    xo, yo = cx + outer * np.cos(t), cy + outer * np.sin(t)
    xi, yi = cx + inner * np.cos(t[::-1]), cy + inner * np.sin(t[::-1])
    ax.fill(np.concatenate([xo, xi]), np.concatenate([yo, yi]),
            color=fill, lw=0, zorder=zorder)
    ax.plot(xo, yo, color=edge, lw=lw, zorder=zorder + 0.1)
    ax.plot(xi[::-1], yi[::-1], color=edge, lw=lw, zorder=zorder + 0.1)


def points(ax, pts, s=9, color=INK, zorder=5, **kw):
    pts = np.asarray(pts, dtype=float)
    ax.scatter(pts[:, 0], pts[:, 1], s=s, color=color, edgecolors="none",
               zorder=zorder, **kw)


# ---------------------------------------------------------------------
# General-position checks
# ---------------------------------------------------------------------

def orient(a, b, c):
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])


def incircle(a, b, c, d):
    """Zero exactly when a, b, c, d are concyclic (or three are collinear)."""
    ax, ay = a[0] - d[0], a[1] - d[1]
    bx, by = b[0] - d[0], b[1] - d[1]
    cx, cy = c[0] - d[0], c[1] - d[1]
    a2, b2, c2 = ax * ax + ay * ay, bx * bx + by * by, cx * cx + cy * cy
    return float(ax * (by * c2 - b2 * cy) - ay * (bx * c2 - b2 * cx)
                 + a2 * (bx * cy - by * cx))


def gp_margins(pts):
    """Return (min |orientation det| over triples,
               min |incircle det| over quadruples)."""
    pts = [tuple(map(float, p)) for p in pts]
    mo = min(abs(orient(*t)) for t in itertools.combinations(pts, 3))
    mc = min(abs(incircle(*q)) for q in itertools.combinations(pts, 4))
    return mo, mc


def distinct_sq_distances(pts, decimals=9):
    """Number of distinct squared distances (rounded) and the smallest gap
    between consecutive distinct values (to show the rounding is safe)."""
    pts = np.asarray(pts, dtype=float)
    vals = []
    for i in range(len(pts)):
        for j in range(i + 1, len(pts)):
            d = pts[i] - pts[j]
            vals.append(round(float(d @ d), decimals))
    u = sorted(set(vals))
    gap = min(b - a for a, b in zip(u, u[1:])) if len(u) > 1 else float("inf")
    return len(u), gap

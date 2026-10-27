"""fig_distance_lists: why the annulus.  The distances from two points of a
random set to all the other points: (a) uniform in the unit disk, where the
two lists differ; (b) uniform in the annulus sbar <= |p| <= 1, where they
nearly coincide."""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

import figstyle as fs  # noqa: E402

fs.setup()
N = 3000
BINS = np.linspace(0, 2, 41)          # width 0.05


def sample(rng, inner):
    r = np.sqrt(rng.uniform(inner ** 2, 1.0, N))   # uniform in area
    t = rng.uniform(0, 2 * np.pi, N)
    return np.c_[r * np.cos(t), r * np.sin(t)]


def nearest(P, target):
    return int(np.argmin(((P - np.asarray(target)) ** 2).sum(1)))


def polar(r, deg):
    a = np.radians(deg)
    return (r * np.cos(a), r * np.sin(a))


rows = [
    ("(a)", 0.0, 1, [(0.0, 0.0), polar(0.99, -35)]),
    ("(b)", fs.SBAR, 2, [polar(0.92, 40), polar(0.92, 215)]),
]

fig = plt.figure(figsize=(fs.FULL_WIDTH, 4.55))
H = 0.40                         # axes height (figure fraction) per row
bottoms = [0.565, 0.085]
styles = [dict(color=fs.ORANGE, ls="-", lw=1.2),
          dict(color=fs.BLUE, ls=(0, (4, 1.6)), lw=1.2)]
markers = [dict(marker="o", s=34), dict(marker="s", s=28)]
hist_axes = []
for (lab, inner, seed, targets), y0 in zip(rows, bottoms):
    rng = np.random.default_rng(seed)
    P = sample(rng, inner)
    idx = [nearest(P, t) for t in targets]
    # scatter panel
    axs = fig.add_axes([0.0, y0, 0.30, H])
    if inner > 0:
        fs.draw_annulus(axs, inner=inner)
    else:
        fs.draw_annulus(axs, inner=0.0)
    axs.scatter(P[:, 0], P[:, 1], s=0.35, color="#4A4A4A", edgecolors="none",
                zorder=2, rasterized=False)
    for k, (i, st, mk) in enumerate(zip(idx, styles, markers)):
        axs.scatter([P[i, 0]], [P[i, 1]], color=st["color"], edgecolors="white",
                    linewidths=0.9, zorder=6, **mk)
    fs.geometry_axes(axs)
    axs.set_xlim(-1.05, 1.05)
    axs.set_ylim(-1.05, 1.05)
    fs.panel_label(axs, lab, x=0.0, y=1.02)
    # labels for p1, p2
    for k, i in enumerate(idx):
        p = P[i]
        if np.hypot(*p) < 0.3:
            off = np.array([0.10, 0.10])
        else:
            off = -0.20 * p / np.hypot(*p)
        axs.text(p[0] + off[0], p[1] + off[1], rf"$p_{k + 1}$", ha="center",
                 va="center", fontsize=9, zorder=7,
                 bbox=dict(boxstyle="square,pad=0.08", fc="white", ec="none",
                           alpha=0.85))
    # histogram panel
    axh = fig.add_axes([0.42, y0 + 0.02, 0.57, H - 0.04])
    hist_axes.append(axh)
    for k, (i, st) in enumerate(zip(idx, styles)):
        d = np.hypot(*(P - P[i]).T)
        d = np.delete(d, i)
        axh.hist(d, bins=BINS, density=True, histtype="step", label=
                 rf"$p_{k + 1}$, \ $|p_{k + 1}| = {np.hypot(*P[i]):.2f}$", **st)
        msd = (d ** 2).mean()
        print(f"{lab} p{k + 1} = ({P[i, 0]:+.3f}, {P[i, 1]:+.3f}), |p| = "
              f"{np.hypot(*P[i]):.3f}, mean |p-q|^2 = {msd:.3f}, "
              f"max |p-q| = {d.max():.3f}")
    for side in ("top", "right"):
        axh.spines[side].set_visible(False)
    axh.set_xlim(0, 2)
    axh.set_xticks([0, 0.5, 1, 1.5, 2])
    axh.set_xticklabels([r"$0$", r"$0.5$", r"$1$", r"$1.5$", r"$2$"])
    axh.set_ylabel("density")
    axh.legend(loc="upper left", handlelength=2.4, borderaxespad=0.2)
    if y0 == bottoms[-1]:
        axh.set_xlabel(r"distance $|p_i q|$ to the other points $q$")
    else:
        axh.tick_params(labelbottom=True)

ymax = max(a.get_ylim()[1] for a in hist_axes)
for a in hist_axes:
    a.set_ylim(0, ymax)
fs.save(fig, "fig_distance_lists")

"""fig_annulus: (a) the annulus sbar <= |p| <= 1 to scale, with sample
points; (b) the inner radius s(delta) = sqrt(1 - 12 delta)."""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

import figstyle as fs  # noqa: E402

fs.setup()

fig = plt.figure(figsize=(fs.FULL_WIDTH, 2.75))
axa = fig.add_axes([0.0, 0.0, 0.44, 1.0])
axb = fig.add_axes([0.58, 0.16, 0.40, 0.80])

# ------------------------------------------------------------- panel (a)
ax = axa
fs.draw_annulus(ax)
rng = np.random.default_rng(40)
m = 60
rad = np.sqrt(rng.uniform(fs.SBAR2, 1.0, m))      # uniform in area
th = rng.uniform(0, 2 * np.pi, m)
pts = np.c_[rad * np.cos(th), rad * np.sin(th)]
assert np.all((pts ** 2).sum(1) >= fs.SBAR2) and np.all((pts ** 2).sum(1) <= 1)
fs.points(ax, pts, s=6)

# centroid
c = 0.035
ax.plot([-c, c], [-c, c], color=fs.INK, lw=0.8, zorder=6)
ax.plot([-c, c], [c, -c], color=fs.INK, lw=0.8, zorder=6)
ax.text(0.0, -0.08, r"centroid $0$", ha="center", va="top", fontsize=8)

# radius arrows
arrow = dict(arrowstyle="-|>", mutation_scale=7, color=fs.INK, lw=0.7,
             shrinkA=0, shrinkB=0)
a1 = np.radians(40)
ax.annotate("", xy=(np.cos(a1), np.sin(a1)), xytext=(0, 0), arrowprops=arrow,
            zorder=7)
ax.text(0.40 * np.cos(a1) - 0.07, 0.40 * np.sin(a1) + 0.07, r"$1$",
        ha="center", va="center")
a2 = np.radians(152)
tip2 = fs.SBAR * np.array([np.cos(a2), np.sin(a2)])
ax.annotate("", xy=tip2, xytext=(0, 0), arrowprops=arrow, zorder=7)
ax.text(0.47 * np.cos(a2) + 0.03, 0.47 * np.sin(a2) + 0.13,
        r"$\bar s \approx 0.836$", ha="center", va="center")
fs.geometry_axes(ax)
ax.set_xlim(-1.05, 1.05)
ax.set_ylim(-1.05, 1.05)
fs.panel_label(ax, "(a)", x=0.02, y=1.0)

# ------------------------------------------------------------- panel (b)
ax = axb
d = np.linspace(0, 1 / 12, 600)
s = np.sqrt(np.clip(1 - 12 * d, 0, None))
ax.fill_between(d, s, 1.0, color=fs.RING_FILL, lw=0, zorder=0)
ax.plot(d, s, color=fs.BLUE, lw=1.4, zorder=3)
ax.plot([0, 1 / 12], [1, 1], color=fs.RING_EDGE, lw=fs.LW, zorder=2)
d0 = 1 / 40
s0 = np.sqrt(1 - 12 * d0)
print(f"s(1/40) = sqrt(0.7) = {s0:.5f}; sbar = {fs.SBAR:.5f}; omega = {fs.OMEGA:.6f}")
ax.plot([d0, d0], [0, s0], color=fs.ORANGE, lw=0.8, ls=(0, (3, 2)), zorder=2)
ax.plot([0, d0], [s0, s0], color=fs.ORANGE, lw=0.8, ls=(0, (3, 2)), zorder=2)
ax.scatter([d0], [s0], s=22, color=fs.ORANGE, edgecolors="white",
           linewidths=0.8, zorder=5)
ax.text(d0 + 0.0025, s0 + 0.035, r"$s(\tfrac1{40}) = \sqrt{0.7} \approx 0.837$",
        ha="left", va="bottom", fontsize=8)
ax.text(0.0665, 0.70, r"points of $P$" + "\n" + r"lie in the band",
        ha="center", va="center", fontsize=8, color=fs.MUTED, linespacing=1.25)
ax.text(0.051, 0.30, r"$s(\delta) = \sqrt{1 - 12\delta}$", ha="center",
        va="center", fontsize=8)
ax.annotate(r"$\delta = 0$: the band is the unit circle", xy=(0.0003, 1.0),
            xytext=(0.0045, 1.085), ha="left", va="center", fontsize=8,
            arrowprops=dict(arrowstyle="-", color=fs.MUTED, lw=fs.LW_THIN,
                            shrinkA=1, shrinkB=1))
ax.set_xlim(0, 1 / 12)
ax.set_ylim(0, 1.13)
ax.set_xticks([0, 1 / 40, 1 / 20, 1 / 12])
ax.set_xticklabels([r"$0$", r"$\tfrac1{40}$", r"$\tfrac1{20}$", r"$\tfrac1{12}$"])
ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
ax.set_yticklabels([r"$0$", r"$0.25$", r"$0.5$", r"$0.75$", r"$1$"])
ax.set_xlabel(r"$\delta$ \ in \ $D \le (\tfrac13 + \delta)\,n$")
ax.set_ylabel(r"radius $|p|$")
for side in ("top", "right"):
    ax.spines[side].set_visible(False)
fs.panel_label(ax, "(b)", x=-0.20, y=1.06)

fs.save(fig, "fig_annulus")

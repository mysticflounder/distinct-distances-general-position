"""fig_farkas: a toy Farkas certificate.  The half-planes x >= 1, y >= 1 and
x + y <= 1 have no common point; the multipliers 1, 1, 1 prove it."""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Polygon  # noqa: E402

import figstyle as fs  # noqa: E402

fs.setup()
lo, hi = -0.5, 2.5
B = 3.0   # polygons extend past the window and are clipped

# A grid check that no point satisfies all three inequalities.
g = np.linspace(lo, hi, 601)
X, Y = np.meshgrid(g, g)
assert not np.any((X >= 1) & (Y >= 1) & (X + Y <= 1))

fig = plt.figure(figsize=(fs.FULL_WIDTH, 2.55))
ax = fig.add_axes([0.075, 0.13, 0.36, 0.83])

half_planes = [
    # polygon, colour, tint, hatch
    ([(1, -B), (B, -B), (B, B), (1, B)], fs.ORANGE, fs.ORANGE_TINT, "////"),
    ([(-B, 1), (B, 1), (B, B), (-B, B)], fs.BLUE, fs.BLUE_TINT, "\\\\\\\\"),
    ([(-B, -B), (-B, 1 + B), (1 + B, -B)], fs.GREEN, fs.GREEN_TINT, "...."),
]
for poly, col, tint, hatch in half_planes:
    ax.add_patch(Polygon(poly, closed=True, facecolor=tint, alpha=0.45,
                         edgecolor="none", zorder=1))
    ax.add_patch(Polygon(poly, closed=True, facecolor="none", hatch=hatch,
                         edgecolor=col, hatchcolor=col, lw=0, alpha=0.55,
                         zorder=2))
# boundary lines
ax.plot([1, 1], [lo, hi], color=fs.ORANGE, lw=1.3, zorder=3)
ax.plot([lo, hi], [1, 1], color=fs.BLUE, lw=1.3, ls=(0, (4, 1.6)), zorder=3)
ax.plot([lo, hi], [1 - lo, 1 - hi], color=fs.GREEN, lw=1.3, ls=(0, (1.2, 1.2)),
        zorder=3)
box = dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.9)
ax.text(1.75, 0.35, r"$x \ge 1$", ha="center", va="center", fontsize=9,
        bbox=box, zorder=5)
ax.text(0.15, 2.05, r"$y \ge 1$", ha="center", va="center", fontsize=9,
        bbox=box, zorder=5)
ax.text(-0.05, 0.25, r"$x + y \le 1$", ha="center", va="center", fontsize=9,
        bbox=box, zorder=5)
ax.set_xlim(lo, hi)
ax.set_ylim(lo, hi)
ax.set_aspect("equal")
ax.set_xticks([0, 1, 2])
ax.set_yticks([0, 1, 2])
ax.set_xticklabels([r"$0$", r"$1$", r"$2$"])
ax.set_yticklabels([r"$0$", r"$1$", r"$2$"])
ax.set_xlabel(r"$x$")
ax.set_ylabel(r"$y$", rotation=0, labelpad=8)
for side in ("top", "right"):
    ax.spines[side].set_visible(False)

# the certificate, in the margin
cert = (r"$\begin{array}{r@{\;}c@{\;}r@{\qquad}l}"
        r"-x & \le & -1 & \text{multiplier } 1\\"
        r"-y & \le & -1 & \text{multiplier } 1\\"
        r"x + y & \le & 1 & \text{multiplier } 1\\[2pt] \hline"
        r"\rule{0pt}{11pt}0 & \le & -1 & \text{the sum}"
        r"\end{array}$")
fig.text(0.52, 0.70, cert, ha="left", va="center", fontsize=9)
fig.text(0.52, 0.40, r"No point lies in all three half-planes.",
         ha="left", va="center", fontsize=9)

fs.save(fig, "fig_farkas")

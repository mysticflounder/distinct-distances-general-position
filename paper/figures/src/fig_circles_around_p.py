"""fig_circles_around_p: the trivial bound n - 1 <= 3D.  The point p and
nine other points on three circles centred at p, three on each circle."""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

import figstyle as fs  # noqa: E402

fs.setup()

radii = [1.0, 1.6, 2.3]
# Angles from a seeded search (rng seed 7, 6000 trials) that maximizes the
# smaller of the two general-position margins printed below.
angles_deg = [[95, 192, 304], [68, 158, 284], [20, 139, 263]]
pts = [np.zeros(2)]
for r, angs in zip(radii, angles_deg):
    for a in np.radians(angs):
        pts.append(r * np.array([np.cos(a), np.sin(a)]))
pts = np.array(pts)
mo, mc = fs.gp_margins(pts)
print(f"10 points (p = 0, radii {radii}): min |orientation det| = {mo:.4f}, "
      f"min |incircle det| = {mc:.4f}")
# Each circle centred at p holds exactly 3 points; a fourth would be
# concyclic.  The margins above check general position for all 10 points.
assert mo > 0.2 and mc > 0.2

fig, ax = plt.subplots(figsize=(2.9, 2.9))
tt = np.linspace(0, 2 * np.pi, 600)
for r in radii:
    ax.plot(r * np.cos(tt), r * np.sin(tt), color=fs.CONSTR, lw=fs.LW_THIN,
            zorder=1)
# one radius per circle, labelled
for i, (r, angs) in enumerate(zip(radii, angles_deg)):
    a = np.radians(angs[0])
    q = r * np.array([np.cos(a), np.sin(a)])
    ax.plot([0, q[0]], [0, q[1]], color=fs.CONSTR, lw=fs.LW_THIN, ls=(0, (3, 2)),
            zorder=1)
fs.points(ax, pts[1:], s=14)
ax.scatter([0], [0], s=30, color=fs.ORANGE, edgecolors="white", linewidths=0.8,
           zorder=6)
ax.text(0.10, -0.20, r"$p$", ha="left", va="center")

# radius labels beside the dashed radii, between consecutive circles
for i, (rad, side) in enumerate(zip([0.55, 1.30, 1.95], [+1, -1, +1])):
    a = np.radians(angles_deg[i][0])
    u = np.array([np.cos(a), np.sin(a)])
    nrm = side * np.array([-u[1], u[0]])
    m = rad * u + 0.15 * nrm
    ax.text(m[0], m[1], rf"$r_{i + 1}$", ha="center", va="center",
            fontsize=8, color=fs.MUTED)

# annotation on the outer circle
a = np.radians(50)
tip = 2.3 * np.array([np.cos(a), np.sin(a)])
ax.annotate(r"at most $3$ points", xy=tip, xytext=(1.30, 2.45),
            ha="left", va="bottom", fontsize=8,
            arrowprops=dict(arrowstyle="-", color=fs.MUTED, lw=fs.LW_THIN,
                            shrinkA=0, shrinkB=1))

fs.geometry_axes(ax)
ax.set_xlim(-2.45, 2.45)
ax.set_ylim(-2.45, 2.75)
fs.save(fig, "fig_circles_around_p")

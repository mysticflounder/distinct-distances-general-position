"""fig_examples: the 4x4 grid (D = 9), the regular 16-gon (D = 8) and 16
points in general position (D = 120)."""
import itertools
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

import figstyle as fs  # noqa: E402

fs.setup()

# (a) 4x4 grid: exact integer squared distances.
grid = [(i, j) for i in range(4) for j in range(4)]
sq = sorted({(p[0] - q[0]) ** 2 + (p[1] - q[1]) ** 2
             for p, q in itertools.combinations(grid, 2)})
D_grid = len(sq)
assert sq == [1, 2, 4, 5, 8, 9, 10, 13, 18], sq
triple = [(0, 1), (1, 2), (2, 3)]   # collinear: on the line y = x + 1
assert fs.orient(*triple) == 0
# exactly three grid points on that line
assert sum(1 for p in grid if p[1] == p[0] + 1) == 3

# (b) regular 16-gon.
m = 16
ang = 2 * np.pi * np.arange(m) / m + np.pi / m
gon = np.c_[np.cos(ang), np.sin(ang)]
D_gon, gap_gon = fs.distinct_sq_distances(gon)

# (c) 16 points in general position.  Points are added one at a time from a
# seeded stream; a candidate is rejected if it is too close to a previous
# point or makes a nearly collinear triple or a nearly concyclic quadruple.
rng = np.random.default_rng(20261026)
T_ORIENT, T_INCIRC, T_SEP = 0.03, 0.005, 0.30
gp = []
tries = 0
while len(gp) < 16:
    tries += 1
    c = rng.uniform(-1, 1, 2)
    if c @ c > 1:
        continue
    if any(np.hypot(*(c - q)) < T_SEP for q in gp):
        continue
    if any(abs(fs.orient(a, b, c)) < T_ORIENT
           for a, b in itertools.combinations(gp, 2)):
        continue
    if any(abs(fs.incircle(a, b, d, c)) < T_INCIRC
           for a, b, d in itertools.combinations(gp, 3)):
        continue
    gp.append(c)
gp = np.array(gp)
mo, mc = fs.gp_margins(gp)
D_gp, gap_gp = fs.distinct_sq_distances(gp)
assert D_gp == 120

print(f"(a) grid: D = {D_grid}, squared distances {sq}")
print(f"(b) 16-gon: D = {D_gon}, min gap between distinct squared distances {gap_gon:.4f}")
print(f"(c) GP set (in unit disk, seed 20261026, {tries} candidates): D = {D_gp}, "
      f"min gap between distinct squared distances {gap_gp:.2e}")
print(f"    min |orientation det| = {mo:.4f}, min |incircle det| = {mc:.4f}")

# ---------------------------------------------------------------------
fig, axes = plt.subplots(1, 3, figsize=(fs.FULL_WIDTH, 2.25))
fig.subplots_adjust(left=0.01, right=0.99, bottom=0.12, top=0.97, wspace=0.12)
R = 1.18   # common half-width of every panel's window

# (a): grid centred at 0, scaled to the unit-disk size of the others.
ax = axes[0]
g = (np.array(grid, float) - 1.5) / 1.5 * 0.80
t3 = (np.array(triple, float) - 1.5) / 1.5 * 0.80
d = t3[2] - t3[0]
ext = np.array([t3[0] - 0.30 * d, t3[2] + 0.30 * d])
ax.plot(ext[:, 0], ext[:, 1], color=fs.ORANGE, lw=fs.LW, zorder=2)
fs.points(ax, g, s=14)
fs.points(ax, t3, s=22, color=fs.ORANGE, zorder=6)

# (b)
ax = axes[1]
tt = np.linspace(0, 2 * np.pi, 400)
ax.plot(np.cos(tt), np.sin(tt), color=fs.ORANGE, lw=fs.LW, zorder=2)
fs.points(ax, gon, s=14)

# (c)
ax = axes[2]
fs.points(ax, gp, s=14)

for ax, lab, D in zip(axes, "abc", (D_grid, D_gon, D_gp)):
    fs.geometry_axes(ax)
    ax.set_xlim(-R, R)
    ax.set_ylim(-R, R)
    ax.text(0.0, -R - 0.02, rf"({lab})\quad $D = {D}$", ha="center", va="top")

fs.save(fig, "fig_examples")

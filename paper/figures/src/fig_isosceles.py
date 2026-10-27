"""fig_isosceles: (a) at most two apexes on bis(a,b); (b) the proof that
kappa is injective (Lemma 2): a, sigma a, b, sigma b on one circle."""
import os
import sys
from fractions import Fraction as Fr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

import figstyle as fs  # noqa: E402

fs.setup()


def ticks(ax, p, q, k, size=0.10, gap=0.055, color=fs.INK):
    """k equal-length tick marks across the midpoint of segment pq."""
    p, q = np.asarray(p, float), np.asarray(q, float)
    u = (q - p) / np.hypot(*(q - p))
    nrm = np.array([-u[1], u[0]])
    mid = (p + q) / 2
    for j in range(k):
        c = mid + (j - (k - 1) / 2) * gap * u
        e = np.array([c - size / 2 * nrm, c + size / 2 * nrm])
        ax.plot(e[:, 0], e[:, 1], color=color, lw=0.8, zorder=4,
                solid_capstyle="butt")


fig, (axa, axb) = plt.subplots(1, 2, figsize=(fs.FULL_WIDTH, 2.95))
fig.subplots_adjust(left=0.0, right=1.0, bottom=0.0, top=1.0, wspace=0.08)

# ------------------------------------------------------------- panel (a)
ax = axa
a, b = np.array([-1.0, 0.0]), np.array([1.0, 0.0])
z, zp = np.array([0.0, 1.45]), np.array([0.0, -1.05])
assert np.isclose(np.hypot(*(z - a)), np.hypot(*(z - b)))
assert np.isclose(np.hypot(*(zp - a)), np.hypot(*(zp - b)))
ax.plot([0, 0], [-1.33, 2.05], color=fs.BLUE, lw=fs.LW, zorder=1)
ax.text(0.07, 2.0, r"$\mathrm{bis}(a,b)$", ha="left", va="top")
ax.plot([a[0], b[0]], [a[1], b[1]], color=fs.CONSTR, lw=fs.LW_THIN, zorder=1)
for apex in (z, zp):
    for e in (a, b):
        ax.plot([apex[0], e[0]], [apex[1], e[1]], color=fs.INK, lw=0.7, zorder=2)
ticks(ax, z, a, 1)
ticks(ax, z, b, 1)
ticks(ax, zp, a, 2)
ticks(ax, zp, b, 2)
fs.points(ax, [a, b], s=16)
fs.points(ax, [z, zp], s=22, color=fs.ORANGE, zorder=6)
ax.text(a[0] - 0.08, a[1], r"$a$", ha="right", va="center")
ax.text(b[0] + 0.08, b[1], r"$b$", ha="left", va="center")
ax.text(z[0] + 0.1, z[1] + 0.04, r"$z$", ha="left", va="bottom")
ax.text(zp[0] + 0.1, zp[1] - 0.04, r"$z'$", ha="left", va="top")
ax.text(0.0, -1.47, "a third apex would be a third\n"
        r"point of $P$ on the line $\mathrm{bis}(a,b)$",
        ha="center", va="top", fontsize=8, color=fs.MUTED, linespacing=1.3)
fs.geometry_axes(ax)
ax.set_xlim(-1.6, 1.6)
ax.set_ylim(-1.75, 2.12)
fs.panel_label(ax, "(a)")

# ------------------------------------------------------------- panel (b)
ax = axb
# exact rational coordinates; ell is the line x = 0, sigma(x, y) = (-x, y)
A = (Fr(-21, 20), Fr(11, 20))
B = (Fr(-9, 20), Fr(-15, 20))
sA, sB = (-A[0], A[1]), (-B[0], B[1])
# c = bis(a, b) meets ell: |c - a|^2 = |c - b|^2 with c = (0, y)
cy = (A[0] ** 2 + A[1] ** 2 - B[0] ** 2 - B[1] ** 2) / (2 * (A[1] - B[1]))
C = (Fr(0), cy)
r2 = {tuple(P): (P[0] - C[0]) ** 2 + (P[1] - C[1]) ** 2 for P in (A, sA, B, sB)}
assert len(set(r2.values())) == 1, r2   # exact: the four points are concyclic
rad2 = next(iter(r2.values()))
print(f"(b) c = (0, {cy}) = (0, {float(cy):.5f}); |ca|^2 = {rad2} "
      f"(radius {float(rad2) ** 0.5:.5f}); equal for a, sigma a, b, sigma b")

af, bf, saf, sbf, cf = (np.array([float(P[0]), float(P[1])])
                        for P in (A, B, sA, sB, C))
z, zp = np.array([0.0, 1.78]), np.array([0.0, -1.32])
assert abs(np.hypot(*(z - cf)) ** 2 - float(rad2)) > 0.1
assert abs(np.hypot(*(zp - cf)) ** 2 - float(rad2)) > 0.1

ax.plot([0, 0], [-1.6, 2.05], color=fs.BLUE, lw=fs.LW, zorder=1)
ax.text(0.07, 2.0, r"$\ell$", ha="left", va="top")
# the pairs f and g (perpendicular to ell)
for P, Q in ((af, saf), (bf, sbf)):
    ax.plot([P[0], Q[0]], [P[1], Q[1]], color=fs.CONSTR, lw=fs.LW_THIN, zorder=1)
# bis(a, b): through the midpoint of ab, perpendicular to ab, through c
mid = (af + bf) / 2
d = bf - af
u = np.array([-d[1], d[0]]) / np.hypot(*d)
assert abs((cf - mid)[0] * u[1] - (cf - mid)[1] * u[0]) < 1e-12      # c lies on bis(a, b)
seg = np.array([mid - 0.95 * u, mid + 1.25 * u])
ax.plot(seg[:, 0], seg[:, 1], color=fs.CONSTR, lw=fs.LW_THIN,
        ls=(0, (4, 2)), zorder=1)
ax.plot([af[0], bf[0]], [af[1], bf[1]], color=fs.CONSTR, lw=fs.LW_THIN, zorder=1)
ax.text(seg[0][0] - 0.05, seg[0][1] - 0.02, r"$\mathrm{bis}(a,b)$",
        ha="right", va="top", fontsize=8, color=fs.MUTED)
# the circle through the four points, centred at c
tt = np.linspace(0, 2 * np.pi, 600)
rad = float(rad2) ** 0.5
ax.plot(cf[0] + rad * np.cos(tt), cf[1] + rad * np.sin(tt), color=fs.ORANGE,
        lw=1.1, ls=(0, (5, 2.5)), zorder=2)
fs.points(ax, [af, saf, bf, sbf], s=16)
fs.points(ax, [z, zp], s=22, color=fs.ORANGE, zorder=6)
ax.scatter([cf[0]], [cf[1]], s=14, facecolors="white", edgecolors=fs.INK,
           linewidths=0.8, zorder=6)
ax.text(af[0] - 0.08, af[1] + 0.03, r"$a$", ha="right", va="bottom")
ax.text(saf[0] + 0.08, saf[1] + 0.03, r"$\sigma a$", ha="left", va="bottom")
ax.text(bf[0] - 0.08, bf[1] - 0.03, r"$b$", ha="right", va="top")
ax.text(sbf[0] + 0.08, sbf[1] - 0.03, r"$\sigma b$", ha="left", va="top")
ax.text(cf[0] + 0.09, cf[1] + 0.02, r"$c$", ha="left", va="bottom")
ax.text(z[0] + 0.1, z[1], r"$z$", ha="left", va="center")
ax.text(zp[0] + 0.1, zp[1], r"$z'$", ha="left", va="center")
ax.text(0.5 * (af[0] + saf[0]) + 0.45, af[1] + 0.03, r"$f$", ha="center",
        va="bottom", fontsize=8, color=fs.MUTED)
ax.text(0.5 * (bf[0] + sbf[0]) + 0.25, bf[1] + 0.03, r"$g$", ha="center",
        va="bottom", fontsize=8, color=fs.MUTED)
fs.geometry_axes(ax)
ax.set_xlim(-1.6, 1.6)
ax.set_ylim(-1.75, 2.12)
fs.panel_label(ax, "(b)")

fs.save(fig, "fig_isosceles")

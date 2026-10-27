"""fig_chord: Lemma C (the same-side chord).  The circle C(p, r) meets the
annulus in two arcs, mirror images in the line through 0 and p; two of any
three points of P on C(p, r) lie on one arc, so they are close."""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Rectangle, ConnectionPatch, Arc  # noqa: E402

import figstyle as fs  # noqa: E402

fs.setup()

a, r = 0.92, 0.8
p = np.array([a, 0.0])
L = (fs.SBAR2 - a * a - r * r) / (2 * a * r)
w = (1 - fs.SBAR2) / (2 * a * r)


def cl(x):
    return min(1.0, max(-1.0, x))


phi_lo, phi_hi = np.arccos(cl(L + w)), np.arccos(cl(L))
f = phi_hi - phi_lo


def q(phi):
    return p + r * np.array([np.cos(phi), np.sin(phi)])


# check the arc ends lie on the two boundary circles
assert np.isclose(np.hypot(*q(phi_lo)), 1.0) and np.isclose(np.hypot(*q(phi_hi)), fs.SBAR)
q1, q2 = q(phi_lo + 0.18 * f), q(phi_lo + 0.86 * f)
q3 = q(-(phi_lo + 0.45 * f))
for x in (q1, q2, q3):
    assert fs.SBAR <= np.hypot(*x) <= 1.0
chord = np.hypot(*(q1 - q2))
bound = 2 * r * np.sin(f / 2)
assert chord <= bound
print(f"a = {a}, r = {r}: L = {L:.5f}, w = {w:.5f}, L + w = {L + w:.5f}")
print(f"arc: phi in [{np.degrees(phi_lo):.3f}, {np.degrees(phi_hi):.3f}] deg, "
      f"f = {f:.5f} rad = {np.degrees(f):.3f} deg")
print(f"|q1 q2| = {chord:.5f} <= 2 r sin(f/2) = {bound:.5f}; "
      f"squared bound 2 r^2 (1 - cos f) = {2 * r * r * (1 - np.cos(f)):.5f}")

fig = plt.figure(figsize=(5.94, 2.43))
axa = fig.add_axes([0.0, 0.0, 0.53, 1.0])
axb = fig.add_axes([0.585, 0.02, 0.415, 0.96])

tt = np.linspace(0, 2 * np.pi, 1200)
circ = p[None, :] + r * np.c_[np.cos(tt), np.sin(tt)]


def draw_common(ax, label_pts=True):
    fs.draw_annulus(ax)
    ax.plot([-1.1, 1.85], [0, 0], color=fs.CONSTR, lw=fs.LW_THIN,
            ls=(0, (4, 2)), zorder=1)
    ax.plot(circ[:, 0], circ[:, 1], color=fs.CONSTR, lw=fs.LW, zorder=2)
    for sgn in (+1, -1):
        ph = sgn * np.linspace(phi_lo, phi_hi, 200)
        arc = p[None, :] + r * np.c_[np.cos(ph), np.sin(ph)]
        ax.plot(arc[:, 0], arc[:, 1], color=fs.ORANGE, lw=fs.LW_THICK,
                solid_capstyle="butt", zorder=3)
    ax.plot([q1[0], q2[0]], [q1[1], q2[1]], color=fs.BLUE, lw=1.1, zorder=4,
            solid_capstyle="butt")
    fs.points(ax, [q1, q2, q3], s=12, zorder=6)


# ------------------------------------------------------------- panel (a)
ax = axa
draw_common(ax)
# rays from p to the ends of the upper arc and the angle f
for ph in (phi_lo, phi_hi):
    e = q(ph)
    ax.plot([p[0], e[0]], [p[1], e[1]], color=fs.INK, lw=fs.LW_THIN,
            zorder=2)
ax.add_patch(Arc(p, 0.9, 0.9, theta1=np.degrees(phi_lo),
                 theta2=np.degrees(phi_hi), color=fs.INK, lw=0.6, zorder=3))
mid = (phi_lo + phi_hi) / 2
lf = p + 0.47 * np.array([np.cos(phi_lo - 0.13), np.sin(phi_lo - 0.13)])
ax.text(lf[0], lf[1], r"$f$", ha="center", va="center", fontsize=9)
fs.points(ax, [p], s=14, zorder=6)
ax.text(p[0] + 0.05, p[1] - 0.05, r"$p$", ha="left", va="top")
c = 0.03
ax.plot([-c, c], [-c, c], color=fs.INK, lw=0.7, zorder=6)
ax.plot([-c, c], [c, -c], color=fs.INK, lw=0.7, zorder=6)
ax.text(0.0, -0.06, r"$0$", ha="center", va="top")
ax.text(q3[0] - 0.07, q3[1] - 0.02, r"$q_3$", ha="right", va="top")
cp = q(np.radians(-40))
ax.text(cp[0] + 0.05, cp[1] - 0.05, r"$C(p,r)$", ha="left", va="top")
# radius label on the lower half
rm = p + 0.5 * r * np.array([np.cos(np.radians(-150)), np.sin(np.radians(-150))])
ax.plot([p[0], q(np.radians(-150))[0]], [p[1], q(np.radians(-150))[1]],
        color=fs.CONSTR, lw=fs.LW_THIN, ls=(0, (2, 1.5)), zorder=2)
ax.text(rm[0] + 0.02, rm[1] - 0.05, r"$r$", ha="center", va="top", fontsize=9)
# zoom box
zc = q(mid)
hw, hh = 0.15, 0.135
zx0, zx1, zy0, zy1 = zc[0] - hw, zc[0] + hw, zc[1] - hh, zc[1] + hh
ax.add_patch(Rectangle((zx0, zy0), zx1 - zx0, zy1 - zy0, fill=False,
                       ec=fs.MUTED, lw=0.5, zorder=7))
fs.geometry_axes(ax)
ax.set_xlim(-1.05, 1.76)
ax.set_ylim(-1.05, 1.05)
fs.panel_label(ax, "(a)", x=0.0, y=1.0)

# ------------------------------------------------------------- panel (b)
ax = axb
draw_common(ax)
for ph in (phi_lo, phi_hi):
    e = q(ph)
    ax.plot([p[0], e[0]], [p[1], e[1]], color=fs.INK, lw=fs.LW_THIN, zorder=2)
ax.text(q1[0] + 0.012, q1[1] + 0.006, r"$q_1$", ha="left", va="bottom")
ax.text(q2[0] - 0.012, q2[1] + 0.006, r"$q_2$", ha="right", va="bottom")
cm = (q1 + q2) / 2
ax.annotate(r"short chord", xy=cm, xytext=(zc[0] - 0.035, zc[1] + 0.105),
            ha="center", va="center", fontsize=8,
            arrowprops=dict(arrowstyle="-", color=fs.MUTED, lw=fs.LW_THIN,
                            shrinkA=2, shrinkB=1))
# Boundary labels: |q| = 1 in the white region outside the annulus (top
# right of the window), |q| = sbar in the hole, left of the ray from p.
lab_out = np.array([zx1 - 0.055, zy1 - 0.022])
lab_in = np.array([zx0 + 0.046, zy0 + 0.030])
assert np.hypot(*lab_out) > 1.04 and np.hypot(*lab_in) < fs.SBAR - 0.04
ax.text(*lab_out, r"$|q| = 1$", ha="center", va="center", fontsize=8,
        color=fs.MUTED)
ax.text(*lab_in, r"$|q| = \bar s$", ha="center", va="center", fontsize=8,
        color=fs.MUTED)
fs.geometry_axes(ax)
ax.set_xlim(zx0, zx1)
ax.set_ylim(zy0, zy1)
ax.add_patch(Rectangle((zx0, zy0), zx1 - zx0, zy1 - zy0, fill=False,
                       ec=fs.MUTED, lw=0.5, zorder=7, clip_on=False))
ax.text(zx0 + 0.008, zy1 - 0.008, "(b)", ha="left", va="top", fontsize=9,
        zorder=8)

for (xa, ya), (xb, yb) in (((zx1, zy1), (zx0, zy1)), ((zx1, zy0), (zx0, zy0))):
    con = ConnectionPatch(xyA=(xa, ya), coordsA=axa.transData,
                          xyB=(xb, yb), coordsB=axb.transData,
                          color=fs.MUTED, lw=0.4, ls=(0, (2, 2)))
    fig.add_artist(con)

fs.save(fig, "fig_chord")

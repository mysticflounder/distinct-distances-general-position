"""fig_bisector_gap: Lemma J.  For a long pair {a, b} in the annulus, bis(a,b)
passes within h0 of 0 and meets the annulus in two short segments; the pair
{z, z'} on it is either short (same segment) or long (opposite segments)."""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

import figstyle as fs  # noqa: E402

fs.setup()

lam = 1.2
h0 = fs.OMEGA / (2 * lam)
sigma0 = np.sqrt(1 - h0 ** 2) - np.sqrt(fs.SBAR2 - h0 ** 2)
Lambda0 = 2 * np.sqrt(fs.SBAR2 - h0 ** 2)


def polar(r, deg):
    t = np.radians(deg)
    return r * np.array([np.cos(t), np.sin(t)])


a, b = polar(0.95, 200), polar(0.86, 320)
ab = np.hypot(*(b - a))
assert ab >= lam
nrm = (b - a) / ab
off = (b @ b - a @ a) / (2 * ab)     # bis(a,b) = {x : <x, nrm> = off}
h = abs(off)
c = off * nrm                        # foot of the perpendicular from 0
u = np.array([-nrm[1], nrm[0]])
if u[1] < 0:
    u = -u                           # u points "up"
t_in, t_out = np.sqrt(fs.SBAR2 - h * h), np.sqrt(1 - h * h)
sig_h, Lam_h = t_out - t_in, 2 * t_in
assert h <= h0 and sig_h <= sigma0 and Lam_h >= Lambda0


def on(t):
    return c + t * u


# panel (a): same segment; panel (b): opposite segments
za, zpa = on(t_in + 0.14 * sig_h), on(t_in + 0.86 * sig_h)
zb, zpb = on(t_in + 0.55 * sig_h), on(-(t_in + 0.40 * sig_h))
for z in (za, zpa, zb, zpb):
    assert fs.SBAR - 1e-12 <= np.hypot(*z) <= 1 + 1e-12
    assert np.isclose(np.hypot(*(z - a)), np.hypot(*(z - b)))
d_same, d_opp = np.hypot(*(za - zpa)), np.hypot(*(zb - zpb))
assert d_same <= sigma0 and d_opp >= Lambda0

print(f"lambda = {lam}, omega = {fs.OMEGA:.6f}")
print(f"h0 = omega/(2 lambda) = {h0:.5f}, sigma0 = {sigma0:.5f}, Lambda0 = {Lambda0:.5f}")
print(f"|a| = {np.hypot(*a):.3f}, |b| = {np.hypot(*b):.3f}, |ab| = {ab:.5f} >= {lam}")
print(f"h = ||b|^2 - |a|^2|/(2|ab|) = {h:.5f} <= h0; sigma(h) = {sig_h:.5f} <= sigma0; "
      f"Lambda(h) = {Lam_h:.5f} >= Lambda0")
print(f"(a) |zz'| = {d_same:.5f} <= sigma0 = {sigma0:.5f}")
print(f"(b) |zz'| = {d_opp:.5f} >= Lambda0 = {Lambda0:.5f}")

fig, axes = plt.subplots(1, 2, figsize=(fs.FULL_WIDTH, 2.95))
fig.subplots_adjust(left=0.0, right=1.0, bottom=0.0, top=1.0, wspace=0.06)
R = 1.10
for ax, (z, zp), lab in zip(axes, ((za, zpa), (zb, zpb)), "ab"):
    fs.draw_annulus(ax)
    # bisector line across the window, and its two segments in the annulus
    ext = np.array([on(-1.12), on(1.12)])
    ax.plot(ext[:, 0], ext[:, 1], color=fs.BLUE, lw=fs.LW_THIN, zorder=2)
    for sgn in (+1, -1):
        s = np.array([on(sgn * t_in), on(sgn * t_out)])
        ax.plot(s[:, 0], s[:, 1], color=fs.ORANGE, lw=fs.LW_THICK,
                solid_capstyle="butt", zorder=3)
    # perpendicular from 0 to the foot c
    ax.plot([0, c[0]], [0, c[1]], color=fs.INK, lw=0.6, zorder=4)
    cr = 0.022
    ax.plot([-cr, cr], [-cr, cr], color=fs.INK, lw=0.7, zorder=6)
    ax.plot([-cr, cr], [cr, -cr], color=fs.INK, lw=0.7, zorder=6)
    ax.scatter([c[0]], [c[1]], s=12, facecolors="white", edgecolors=fs.INK,
               linewidths=0.7, zorder=6)
    # the apex segments za, zb
    for e in (a, b):
        ax.plot([z[0], e[0]], [z[1], e[1]], color=fs.INK, lw=0.5,
                ls=(0, (3, 2)), zorder=2)
    fs.points(ax, [a, b], s=14, zorder=6)
    ax.scatter([z[0], zp[0]], [z[1], zp[1]], s=11, color=fs.INK,
               edgecolors="white", linewidths=0.6, zorder=7)
    ax.text(a[0] - 0.05, a[1] - 0.03, r"$a$", ha="right", va="top")
    ax.text(b[0] + 0.05, b[1] - 0.03, r"$b$", ha="left", va="top")
    ax.text(0.04, -0.035, r"$0$", ha="left", va="top", fontsize=8)
    ax.text(c[0] - 0.04, c[1] + 0.01, r"$c$", ha="right", va="bottom", fontsize=8)
    bot = on(-1.12)
    ax.text(bot[0] - 0.04, bot[1] + 0.02, r"$\mathrm{bis}(a,b)$", ha="right",
            va="bottom", fontsize=8)
    fs.geometry_axes(ax)
    ax.set_xlim(-R, R)
    ax.set_ylim(-R - 0.03, R + 0.03)

# panel labels carry the two cases; labels of z, z'
fs.panel_label(axes[0], r"(a) \ $|zz'| \le \sigma_0$", x=0.0, y=1.0)
fs.panel_label(axes[1], r"(b) \ $|zz'| \ge \Lambda_0$", x=0.0, y=1.0)
ax = axes[0]
ax.text(za[0] - 0.05, za[1] + 0.03, r"$z$", ha="right", va="center")
ax.text(zpa[0] + 0.045, zpa[1] + 0.01, r"$z'$", ha="left", va="center")
ax = axes[1]
ax.text(zb[0] + 0.05, zb[1] + 0.01, r"$z$", ha="left", va="center")
ax.text(zpb[0] + 0.05, zpb[1], r"$z'$", ha="left", va="center")

fs.save(fig, "fig_bisector_gap")

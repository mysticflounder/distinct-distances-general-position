"""fig_profile: the floating-point LP point at delta = 1/39, ybar = 0.154.
(a) A(rho), 3R(rho) and X(rho) = 3R(rho) - A(rho);
(b) the fraction of pairs per distance bin of width 0.05."""
import bisect
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

import figstyle as fs  # noqa: E402

fs.setup()

with open(os.path.join(fs.DATADIR, "profile_1_39.json")) as fh:
    d = json.load(fh)
G, A, R = (np.array(d[k]) for k in ("G", "A", "R"))
A = np.maximum(A, 0.0)            # the LP point stores -0.0
X = 3 * R - A
qq = 3 * d["dd"] - 1 + d["ee"]
Gl = list(G)


def A_at(x):
    """A at the largest grid point <= x (A = 0 below the grid)."""
    i = bisect.bisect_right(Gl, x + 1e-12) - 1
    return A[i] if i >= 0 else 0.0


def frac(lo, hi):
    return (A_at(hi) - A_at(lo)) / A[-1]


print(f"delta = {d['delta']:.6f}, ybar = {d['ybar0']}, A(2) = {A[-1]:.6f}, "
      f"3R(2) = {3 * R[-1]:.6f}, max X = {X.max():.6f}, q = {qq:.6f}")
first = G[np.argmax(A > 0)]
print(f"A = 0 for rho <= {G[np.argmax(A > 0) - 1]:.6f}; first positive at {first:.6f}")
for lo, hi in ((0, 0.17), (0.17, 0.2), (0.7, 0.8), (1.5, 1.8)):
    print(f"fraction of pairs in ({lo}, {hi}] = {frac(lo, hi):.4f}")

edges = np.round(np.arange(0, 2.0001, 0.05), 10)
bins = np.array([frac(lo, hi) for lo, hi in zip(edges[:-1], edges[1:])])
assert abs(bins.sum() - 1) < 1e-9

fig = plt.figure(figsize=(fs.FULL_WIDTH, 2.75))
ax1 = fig.add_axes([0.085, 0.42, 0.37, 0.52])
ax2 = fig.add_axes([0.085, 0.14, 0.37, 0.20], sharex=ax1)
ax3 = fig.add_axes([0.60, 0.14, 0.39, 0.80])

# ---- (a) top: A and 3R
ax1.plot(G, A, color=fs.ORANGE, lw=1.3, drawstyle="steps-post",
         label=r"$\mathcal A(\rho)$")
ax1.plot(G, 3 * R, color=fs.BLUE, lw=1.3, ls=(0, (4, 1.6)),
         drawstyle="steps-post", label=r"$3\mathcal R(\rho)$")
ax1.set_ylim(0, 1.12)
ax1.set_yticks([0, 0.5, 1])
ax1.set_yticklabels([r"$0$", r"$0.5$", r"$1$"])
ax1.legend(loc="upper left", borderaxespad=0.1)
ax1.tick_params(labelbottom=False)
# ---- (a) bottom: X = 3R - A <= q
ax2.plot(G, X, color=fs.GREEN, lw=1.3, drawstyle="steps-post")
ax2.axhline(qq, color=fs.CONSTR, lw=fs.LW_THIN, ls=(0, (3, 2)))
ax2.text(0.03, qq + 0.004, r"$\mathfrak q$", ha="left", va="bottom", fontsize=8,
         color=fs.MUTED)
ax2.set_ylim(0, 0.105)
ax2.set_yticks([0, 0.05])
ax2.set_yticklabels([r"$0$", r"$0.05$"])
ax2.set_ylabel(r"$\mathcal X = 3\mathcal R - \mathcal A$", fontsize=8, labelpad=4)
ax2.set_xlim(0, 2)
ax2.set_xticks([0, 0.5, 1, 1.5, 2])
ax2.set_xticklabels([r"$0$", r"$0.5$", r"$1$", r"$1.5$", r"$2$"])
ax2.set_xlabel(r"$\rho$")
for ax in (ax1, ax2):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
ax1.text(-0.22, 1.0, "(a)", transform=ax1.transAxes, ha="left", va="top")

# ---- (b) bar chart of the distance distribution
ax = ax3
ax.bar(edges[:-1], bins, width=0.05, align="edge", color=fs.BLUE,
       edgecolor="white", linewidth=0.8)
ax.set_xlim(0, 2)
ax.set_xticks([0, 0.5, 1, 1.5, 2])
ax.set_xticklabels([r"$0$", r"$0.5$", r"$1$", r"$1.5$", r"$2$"])
print(f"largest bin: ({edges[bins.argmax()]:.2f}, {edges[bins.argmax() + 1]:.2f}] "
      f"with {bins.max():.4f}")
ax.set_ylim(0, bins.max() + 0.035)
ax.set_yticks([0, 0.05, 0.10, 0.15, 0.20])
ax.set_yticklabels([r"$0$", r"$0.05$", r"$0.10$", r"$0.15$", r"$0.20$"])
ax.set_xlabel(r"distance $|pq|$ (bins of width $0.05$)")
ax.set_ylabel("fraction of pairs")
for side in ("top", "right"):
    ax.spines[side].set_visible(False)
# selective direct labels: the three groups that carry most of the mass
groups = ((0.15, 0.20), (0.70, 0.80), (1.50, 1.80))
for lo, hi in groups:
    sel = (edges[:-1] >= lo - 1e-9) & (edges[1:] <= hi + 1e-9)
    top = bins[sel].max() + 0.008
    f = bins[sel].sum()
    ax.plot([lo + 0.005, lo + 0.005, hi - 0.005, hi - 0.005],
            [top - 0.004, top, top, top - 0.004], color=fs.MUTED, lw=0.5)
    ax.text((lo + hi) / 2, top + 0.003, rf"${100 * f:.0f}\%$", ha="center",
            va="bottom", fontsize=8)
ax.text(-0.17, 1.0, "(b)", transform=ax.transAxes, ha="left", va="top")

fs.save(fig, "fig_profile")

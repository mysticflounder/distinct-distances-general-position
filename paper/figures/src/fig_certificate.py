"""fig_certificate: the 60 case certificates.  Top: the Farkas gap of each
case (log scale) over its interval of ybar; bottom: the rows each case uses.
Data: data/certificate_cases.json, parsed from the output of
verifier/verify_constant_gain_1_40.py (see extract_data.py)."""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.ticker import LogLocator, NullFormatter  # noqa: E402

import figstyle as fs  # noqa: E402

fs.setup()

with open(os.path.join(fs.DATADIR, "certificate_cases.json")) as fh:
    d = json.load(fh)
cases = d["cases"]
assert len(cases) == 60 and d["verified"]
lo = np.array([c["ylo"] for c in cases])
hi = np.array([c["yhi"] for c in cases])
gap = np.array([c["gap"] for c in cases])
rows = np.array([c["rows"] for c in cases])
kmin = int(gap.argmin())
assert kmin == 30, kmin
print(f"60 cases on [0, omega], omega = {d['omega']:.6f}; total rows {rows.sum()}")
print(f"smallest gap: case {kmin}, ybar in [{lo[kmin]:.4f}, {hi[kmin]:.4f}], "
      f"gap = {gap[kmin]:.4e}")
print(f"largest gap: case {gap.argmax()}, gap = {gap.max():.4e}")
print(f"rows per case: min {rows.min()} (case {rows.argmin()}), "
      f"max {rows.max()} (case {rows.argmax()})")

fig = plt.figure(figsize=(fs.FULL_WIDTH, 3.5))
ax1 = fig.add_axes([0.10, 0.44, 0.88, 0.53])
ax2 = fig.add_axes([0.10, 0.12, 0.88, 0.24], sharex=ax1)

# ---- top: gaps as horizontal segments over each interval
pad = 0.0004
for k in range(60):
    col = fs.ORANGE if k == kmin else fs.BLUE
    ax1.plot([lo[k] + pad, hi[k] - pad], [gap[k], gap[k]], color=col,
             lw=2.0, solid_capstyle="butt")
# light connectors between consecutive cases so the trend is readable
mids = (lo + hi) / 2
ax1.plot(mids, gap, color=fs.LIGHT, lw=0.5, zorder=0)
ax1.set_yscale("log")
ax1.set_ylim(1e-6, 3e-2)
ax1.yaxis.set_major_locator(LogLocator(base=10, numticks=10))
ax1.yaxis.set_minor_formatter(NullFormatter())
ax1.set_ylabel(r"gap (log scale)")
ax1.tick_params(labelbottom=False)
ax1.annotate(rf"case {kmin}: $\bar y \in [{lo[kmin]:.4f}, {hi[kmin]:.4f}]$," "\n"
             rf"gap $\approx {gap[kmin] * 1e6:.2f}\cdot 10^{{-6}}$",
             xy=(mids[kmin], gap[kmin]), xytext=(mids[kmin] + 0.035, 4.5e-6),
             ha="left", va="center", fontsize=8,
             arrowprops=dict(arrowstyle="-", color=fs.MUTED, lw=fs.LW_THIN,
                             shrinkA=2, shrinkB=3))
ax1.text(0.005, 1.6e-2, r"each case: one interval of $\bar y$ of width $\omega/60$",
         ha="left", va="center", fontsize=8, color=fs.MUTED)

# ---- bottom: rows per case
w = hi - lo
cols = [fs.ORANGE if k == kmin else fs.BLUE for k in range(60)]
ax2.bar(lo, rows, width=w, align="edge", color=cols, edgecolor="white",
        linewidth=0.6)
ax2.set_ylabel("rows used")
ax2.set_ylim(0, rows.max() * 1.08)
ax2.set_yticks([0, 1000, 2000, 3000, 4000])
ax2.set_yticklabels([r"$0$", r"$1000$", r"$2000$", r"$3000$", r"$4000$"])
ax2.set_xlim(0, d["omega"])
ticks = [0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30]
ax2.set_xticks(ticks)
ax2.set_xticklabels([rf"${t:g}$" for t in ticks])
ax2.set_xlabel(r"$\bar y$ \ (the axis ends at $\omega \approx 0.3004$)")
for ax in (ax1, ax2):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
fig.align_ylabels([ax1, ax2])

fs.save(fig, "fig_certificate")

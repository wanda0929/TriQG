"""
K = 1 vs K = 0.95 comparison plot.

Reads brute_force_results.csv (K = 0.95) and brute_force_results_K1.csv
(K = 1) and produces:
   brute_force_K_comparison.png   -- two-panel figure:
     (a) side-by-side per-axis pass-rate bars
     (b) Pareto cloud of both, overlaid
"""

import csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

matplotlib.rcParams.update({
    "font.family": "serif",
    "mathtext.fontset": "cm",
    "axes.linewidth": 0.8,
    "axes.labelsize": 10,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 9,
})


def load(path):
    with open(path) as f:
        rows = list(csv.DictReader(f))
    for r in rows:
        for k, v in r.items():
            try:
                r[k] = float(v)
            except ValueError:
                pass
    return rows


K095 = load("brute_force_results.csv")
K100 = load("brute_force_results_K1.csv")

AXES = ["a_um", "omega_p_MHz", "ratio", "delta_MHz", "omega_c_MHz"]
LABELS = {
    "a_um":         r"$a$ (μm)",
    "omega_p_MHz":  r"$\Omega_p / (2\pi)$ (MHz)",
    "ratio":        r"$\Omega_R / \Omega_p$",
    "delta_MHz":    r"$\Delta / (2\pi)$ (MHz)",
    "omega_c_MHz":  r"$\Omega_c / (2\pi)$ (MHz)",
}
VALUES = {ax: sorted({r[ax] for r in K095}) for ax in AXES}
THRESH = 0.99

# =====================================================================
# Top row: per-axis pass-rate, K=0.95 vs K=1 side-by-side bars
# =====================================================================
fig, (axs1, axs2) = plt.subplots(2, 5, figsize=(13.5, 5.6),
                                 gridspec_kw={"height_ratios": [1, 1]})
# Make second row a single subplot via span
gs = fig.add_gridspec(2, 5, height_ratios=[1, 1.05])
for r in (axs1, axs2):
    for ax in r:
        ax.remove()
axs_top = [fig.add_subplot(gs[0, i]) for i in range(5)]
ax_bot  = fig.add_subplot(gs[1, :])

for i, ax in enumerate(AXES):
    vals = VALUES[ax]
    rate_K95 = []
    rate_K1  = []
    for v in vals:
        s1 = [r for r in K095 if abs(r[ax] - v) < 1e-9]
        s2 = [r for r in K100 if abs(r[ax] - v) < 1e-9]
        rate_K95.append(sum(1 for r in s1 if r["F_bar"] > THRESH) / len(s1))
        rate_K1.append(sum(1 for r in s2 if r["F_bar"] > THRESH) / len(s2))
    x = np.arange(len(vals))
    w = 0.38
    axs_top[i].bar(x - w / 2, np.array(rate_K95) * 100, w,
                   color="#4477aa", label=r"$K = 0.95$")
    axs_top[i].bar(x + w / 2, np.array(rate_K1)  * 100, w,
                   color="#cc6677", label=r"$K = 1.00$")
    axs_top[i].set_xticks(x)
    axs_top[i].set_xticklabels([f"{v:g}" for v in vals])
    axs_top[i].set_xlabel(LABELS[ax])
    axs_top[i].set_ylim(0, 100)
    axs_top[i].axhline(50, color="grey", lw=0.4, ls=":")
    if i == 0:
        axs_top[i].set_ylabel(r"% cells with $\overline{F}_{OR}>0.99$")
        axs_top[i].legend(loc="upper right", fontsize=8, frameon=True,
                          framealpha=0.9)

fig.suptitle("Marginal pass-rate comparison:  "
             "$K = 0.95$ (sub-$\\pi/4$) vs $K = 1.00$ (canonical $\\pi/4$)",
             y=0.99, fontsize=11)

# =====================================================================
# Bottom panel: Pareto cloud of both
# =====================================================================
for rows, color, lbl in [
    (K095, "#4477aa", r"$K = 0.95$"),
    (K100, "#cc6677", r"$K = 1.00$"),
]:
    xs = [r["t_total_ns"] for r in rows]
    ys = [r["F_bar"]      for r in rows]
    # Faded all
    ax_bot.scatter(
        [x for x, y in zip(xs, ys) if y <= THRESH],
        [y for x, y in zip(xs, ys) if y <= THRESH],
        c=color, alpha=0.18, s=12, edgecolors="none",
    )
    # Bold pass
    ax_bot.scatter(
        [x for x, y in zip(xs, ys) if y > THRESH],
        [y for x, y in zip(xs, ys) if y > THRESH],
        c=color, alpha=0.85, s=22, edgecolor="k", linewidth=0.3,
        label=lbl + f"  ({sum(1 for y in ys if y>THRESH)} cells)",
    )
    # Pareto
    pareto = []
    bestF = 0
    for x, y in sorted(zip(xs, ys), key=lambda p: p[0]):
        if y > bestF:
            pareto.append((x, y)); bestF = y
    if pareto:
        px, py = zip(*pareto)
        ax_bot.plot(px, py, "--", color=color, lw=1.2)
ax_bot.axhline(THRESH, color="red", lw=0.7, ls=":")
ax_bot.set_xlabel(r"Total gate time $T_{\mathrm{tot}}$ (ns)")
ax_bot.set_ylabel(r"$\overline{F}_{OR}$")
ax_bot.set_xlim(150, 1000)
ax_bot.set_ylim(0.85, 1.0)
ax_bot.legend(loc="lower right")
ax_bot.set_title("Pareto cloud:  every cell as a point;  "
                 "Pareto frontier dashed",
                 fontsize=10)

fig.tight_layout(rect=[0, 0, 1, 0.965])
fig.savefig("brute_force_K_comparison.png", dpi=170, bbox_inches="tight")
plt.close(fig)
print("Wrote brute_force_K_comparison.png")

"""
Visualization of the 5-D brute-force OR-gate sweep.

Produces three figures:
  brute_force_axis_marginals.png  -- pass-rate vs each axis value
  brute_force_2D_a_Op.png         -- 2-D marginal heatmap (a vs Omega_p),
                                      best F_bar over (ratio, Delta, Omega_c)
  brute_force_pareto.png          -- gate-time vs fidelity Pareto cloud
                                      colored by robust_score
"""

import csv

import matplotlib.pyplot as plt
import matplotlib
import numpy as np

matplotlib.rcParams.update({
    "font.family": "serif",
    "mathtext.fontset": "cm",
    "axes.linewidth": 0.8,
    "axes.labelsize": 10,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 9,
})

with open("brute_force_robust.csv") as f:
    robust = list(csv.DictReader(f))
with open("brute_force_results.csv") as f:
    all_rows = list(csv.DictReader(f))
for tbl in (robust, all_rows):
    for r in tbl:
        for k, v in r.items():
            try:
                r[k] = float(v)
            except ValueError:
                pass

# =====================================================================
# Helpers
# =====================================================================
AXES = ["a_um", "omega_p_MHz", "ratio", "delta_MHz", "omega_c_MHz"]
LABELS = {
    "a_um":         r"$a$ (μm)",
    "omega_p_MHz":  r"$\Omega_p / (2\pi)$ (MHz)",
    "ratio":        r"$\Omega_R / \Omega_p$",
    "delta_MHz":    r"$\Delta / (2\pi)$ (MHz)",
    "omega_c_MHz":  r"$\Omega_c / (2\pi)$ (MHz)",
}
VALUES = {ax: sorted({r[ax] for r in all_rows}) for ax in AXES}

THRESH = 0.99

# =====================================================================
# 1) Axis marginals: pass-rate
# =====================================================================
fig, axs = plt.subplots(1, 5, figsize=(13.5, 2.8))
for i, ax in enumerate(AXES):
    vals = VALUES[ax]
    rates = []
    medians = []
    maxs = []
    for v in vals:
        sel = [r for r in all_rows if abs(r[ax] - v) < 1e-9]
        Fs = [r["F_bar"] for r in sel]
        rates.append(sum(1 for f in Fs if f > THRESH) / len(Fs))
        medians.append(np.median(Fs))
        maxs.append(max(Fs))
    bars = axs[i].bar(range(len(vals)), [r * 100 for r in rates],
                      color="#4477aa", alpha=0.85)
    axs[i].set_xticks(range(len(vals)))
    axs[i].set_xticklabels([f"{v:g}" for v in vals], rotation=0)
    axs[i].set_xlabel(LABELS[ax])
    axs[i].set_ylim(0, 100)
    axs[i].axhline(50, color="grey", lw=0.5, ls=":")
    for j, b in enumerate(bars):
        axs[i].text(j, rates[j] * 100 + 1.5, f"{rates[j]*100:.0f}",
                    ha="center", va="bottom", fontsize=8)
        axs[i].text(j, 3, f"max\n{maxs[j]:.4f}",
                    ha="center", va="bottom", fontsize=7,
                    color="white", weight="bold")
    if i == 0:
        axs[i].set_ylabel(r"% cells with $\overline{F}_{OR} > 0.99$")
fig.suptitle("Marginal pass-rate (each bar averages over the other four axes)",
             fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.92])
fig.savefig("brute_force_axis_marginals.png", dpi=170,
            bbox_inches="tight")
plt.close(fig)
print("Wrote brute_force_axis_marginals.png")

# =====================================================================
# 2) 2D heatmaps: best F_bar maximized over the other axes
# =====================================================================
pairs = [
    ("a_um", "omega_p_MHz"),
    ("a_um", "delta_MHz"),
    ("omega_p_MHz", "ratio"),
    ("delta_MHz", "ratio"),
]

fig, axs = plt.subplots(1, 4, figsize=(14, 3.2))
for ax_obj, (xn, yn) in zip(axs, pairs):
    xs = VALUES[xn]
    ys = VALUES[yn]
    Z = np.full((len(ys), len(xs)), np.nan)
    for i, yv in enumerate(ys):
        for j, xv in enumerate(xs):
            sel = [r for r in all_rows if abs(r[xn] - xv) < 1e-9
                                       and abs(r[yn] - yv) < 1e-9]
            if sel:
                Z[i, j] = max(r["F_bar"] for r in sel)
    im = ax_obj.imshow(Z, origin="lower", aspect="auto", cmap="viridis",
                       vmin=0.97, vmax=0.995,
                       extent=[-0.5, len(xs) - 0.5, -0.5, len(ys) - 0.5])
    ax_obj.set_xticks(range(len(xs)))
    ax_obj.set_xticklabels([f"{v:g}" for v in xs])
    ax_obj.set_yticks(range(len(ys)))
    ax_obj.set_yticklabels([f"{v:g}" for v in ys])
    ax_obj.set_xlabel(LABELS[xn])
    ax_obj.set_ylabel(LABELS[yn])
    # overlay 0.99 contour and text values
    for i in range(len(ys)):
        for j in range(len(xs)):
            if not np.isnan(Z[i, j]):
                col = "white" if Z[i, j] < 0.984 else "black"
                ax_obj.text(j, i, f"{Z[i,j]:.4f}",
                            ha="center", va="center",
                            fontsize=7, color=col)
    cbar = fig.colorbar(im, ax=ax_obj, pad=0.02)
    cbar.set_label(r"$\max\overline{F}_{OR}$", fontsize=9)
    cbar.ax.tick_params(labelsize=8)
fig.suptitle("Best $\\overline{F}_{OR}$ over each (axis, axis) plane (max over the other 3 axes)",
             fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.93])
fig.savefig("brute_force_2D_marginals.png", dpi=170, bbox_inches="tight")
plt.close(fig)
print("Wrote brute_force_2D_marginals.png")

# =====================================================================
# 3) Pareto:  gate time vs fidelity, colored by robust_score
# =====================================================================
# Use full results table; robustness comes from robust CSV.  Map fidelities
# to robust scores by matching (a, Op, r, D, Oc).
key_to_row = {}
for r in robust:
    key = tuple(r[k] for k in AXES)
    key_to_row[key] = r

fig, ax = plt.subplots(figsize=(7.2, 4.4))
xs = []
ys = []
cs = []
for r in all_rows:
    key = tuple(r[k] for k in AXES)
    rb = key_to_row.get(key)
    if rb:
        cs.append(rb["robust_score"])
    else:
        cs.append(-0.05)  # below scale; we'll fade these
    xs.append(r["t_total_ns"])
    ys.append(r["F_bar"])

cs = np.array(cs)
xs = np.array(xs)
ys = np.array(ys)
# Plot failing cells faded
mask_fail = cs < 0
ax.scatter(xs[mask_fail], ys[mask_fail], c="lightgrey", s=14,
           alpha=0.5, label=r"$\overline{F}_{OR} \leq 0.99$")
# Pass cells colored
sc = ax.scatter(xs[~mask_fail], ys[~mask_fail], c=cs[~mask_fail],
                cmap="plasma", vmin=0.0, vmax=1.0, s=22, edgecolor="k",
                linewidth=0.2, label=r"$\overline{F}_{OR} > 0.99$")
# Pareto frontier
sorted_pass = sorted(
    [(x, y) for x, y, m in zip(xs, ys, mask_fail) if not m],
    key=lambda p: p[0],
)
pareto = []
bestF = 0
for x, y in sorted_pass:
    if y > bestF:
        pareto.append((x, y))
        bestF = y
if pareto:
    px, py = zip(*pareto)
    ax.plot(px, py, "k--", lw=1.0, label="Pareto front")
ax.axhline(0.99, color="red", lw=0.7, ls=":")
ax.set_xlabel("Total gate time $T_\\mathrm{total}$ (ns)")
ax.set_ylabel(r"$\overline{F}_{OR}$")
ax.set_ylim(0.83, 1.0)
ax.set_xlim(150, 1000)
cbar = fig.colorbar(sc, ax=ax)
cbar.set_label(r"robust score (fraction of 1-step neighbors with $\overline{F}>0.99$)")
ax.legend(loc="lower right")
fig.tight_layout()
fig.savefig("brute_force_pareto.png", dpi=170, bbox_inches="tight")
plt.close(fig)
print("Wrote brute_force_pareto.png")

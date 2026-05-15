"""
Plot the error budget produced by a4_K1_error_budget.py.

Generates a stacked bar comparing OR vs CCX, where each segment of
the bar is the basis-infidelity contributed by one source of error.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).parent
data = json.loads((HERE / "a4_K1_error_budget.json").read_text())


def channel_shares(gate_block: dict) -> dict:
    """Return ordered (channel_name -> infidelity) including the
    coherent floor as a fourth bucket."""
    b = gate_block["budget"]
    return {
        "Coherent floor\n(non-removable)":  b["coherent_floor"],
        "Same-species $V_{cc}$ (Cs-Cs)":     b["delta_F_no_Vcc"],
        "Intermediate $|P\\rangle$ decay":   b["delta_F_no_P"],
        "Rydberg decay $|R\\rangle,|r\\rangle$": b["delta_F_no_Ryd"],
    }


or_shares = channel_shares(data["OR_gate"])
ccx_shares = channel_shares(data["CCX_gate"])

# Use a consistent color per channel
COLORS = {
    "Coherent floor\n(non-removable)":   "#7a7a7a",
    "Same-species $V_{cc}$ (Cs-Cs)":     "#d95f02",
    "Intermediate $|P\\rangle$ decay":   "#7570b3",
    "Rydberg decay $|R\\rangle,|r\\rangle$": "#1b9e77",
}
ORDER = list(COLORS.keys())

fig, axes = plt.subplots(
    1, 2, figsize=(11, 5.0), gridspec_kw={"width_ratios": [1.05, 1.0]},
)

# ------ Panel (a): stacked bars on absolute scale ----------------------
ax = axes[0]
gates = ["OR\n(244 ns)", "CCX\n(95 ns)"]
shares_per_gate = [or_shares, ccx_shares]
x = np.arange(len(gates))
bottom = np.zeros(len(gates))
for ch in ORDER:
    vals = np.array([sh[ch] for sh in shares_per_gate])
    ax.bar(x, vals * 1e3, bottom=bottom * 1e3, label=ch,
           color=COLORS[ch], edgecolor="white", linewidth=0.7)
    bottom += vals

# Total infidelity annotation
totals = [data["OR_gate"]["budget"]["baseline_infidelity"],
          data["CCX_gate"]["budget"]["baseline_infidelity"]]
for xi, tot in zip(x, totals):
    ax.text(xi, tot * 1e3 + 0.10, f"$1-\\overline{{F}} = {tot*1e3:.2f}\\times 10^{{-3}}$",
            ha="center", va="bottom", fontsize=10, weight="bold")

ax.set_xticks(x)
ax.set_xticklabels(gates, fontsize=11)
ax.set_ylabel(r"Basis infidelity contribution  $\Delta(1-\overline{F})\ [10^{-3}]$",
              fontsize=11)
ax.set_title("(a)  Absolute error budget", loc="left", fontsize=12, weight="bold")
ax.set_ylim(0, max(totals) * 1.18 * 1e3)
ax.grid(True, axis="y", alpha=0.3, linestyle="--")
ax.set_axisbelow(True)

# ------ Panel (b): percent share of removable error --------------------
ax = axes[1]
# Drop coherent floor for the percentage view: it isn't *removable* by
# turning off a single error channel.
PERCENT_ORDER = [
    "Same-species $V_{cc}$ (Cs-Cs)",
    "Intermediate $|P\\rangle$ decay",
    "Rydberg decay $|R\\rangle,|r\\rangle$",
]
shares_pct = []
for sh in shares_per_gate:
    total = sum(sh[c] for c in PERCENT_ORDER)
    shares_pct.append({c: 100.0 * sh[c] / total if total > 0 else 0.0
                       for c in PERCENT_ORDER})

bottom = np.zeros(len(gates))
for ch in PERCENT_ORDER:
    vals = np.array([sp[ch] for sp in shares_pct])
    bars = ax.bar(x, vals, bottom=bottom, label=ch, color=COLORS[ch],
                  edgecolor="white", linewidth=0.7)
    for xi, v, b in zip(x, vals, bottom):
        if v > 4:
            ax.text(xi, b + v / 2, f"{v:.1f}%", ha="center", va="center",
                    fontsize=10, color="white", weight="bold")
    bottom += vals

ax.set_xticks(x)
ax.set_xticklabels(gates, fontsize=11)
ax.set_ylabel("Share of removable (decoherence + $V_{cc}$) error  [%]",
              fontsize=11)
ax.set_title("(b)  Relative share of removable error",
             loc="left", fontsize=12, weight="bold")
ax.set_ylim(0, 102)
ax.grid(True, axis="y", alpha=0.3, linestyle="--")
ax.set_axisbelow(True)

# Single shared legend
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc="lower center", ncol=4, frameon=False,
           fontsize=9.5, bbox_to_anchor=(0.5, -0.02))

fig.suptitle(
    r"Error budget at $a = 4.0\ \mu$m, $K = 1$  "
    r"(OR + CCX gates, ARC 0 K lifetimes)",

    fontsize=12.5, weight="bold", y=1.00,
)
fig.tight_layout(rect=[0, 0.05, 1, 0.97])
out = HERE / "a4_K1_error_budget.png"
fig.savefig(out, dpi=180, bbox_inches="tight")
print(f"Wrote {out}")

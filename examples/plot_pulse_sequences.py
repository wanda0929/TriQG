"""
Pulse Sequence Plots for OR Gate and CCX Gate
==============================================

Creates publication-quality pulse sequence diagrams following
NeurIPS-2025 plotting conventions (PaperBanana style reference):

1. OR gate  -- broken-axis plot (3 segments) so the short control
   pi-pulses (~10 ns) and the long Hanning probe pulse (~653 ns)
   are all visible at the same visual scale.

2. CCX gate -- standard plot (all pulses are on similar timescales).

Run from the repository root:
    PYTHONPATH=. python3 examples/plot_pulse_sequences.py
"""

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

from triqg.pulses import omega_c, omega_p, omega_R, omega_cc, omega_t1, omega_t2

# =====================================================================
# Publication style -- scaled for 4:3 figures
# =====================================================================
plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.size": 11,
        "axes.titlesize": 14,
        "axes.labelsize": 13,
        "xtick.labelsize": 11,
        "ytick.labelsize": 11,
    }
)

# Color palette -- NeurIPS soft pastels
COLOR_CTRL = "#5B8DBE"  # Steel blue  (control pulses)
COLOR_PROBE = "#E76F51"  # Burnt sienna (probe / target-1)
COLOR_COUPLE = "#7FB685"  # Sage green  (coupling / target-2)

LW = 2.8  # line width for pulse traces
AXIS_LW = 1.4  # visible axis/spine thickness
TICK_W = 1.2  # tick mark thickness
BREAK_LW = 1.2  # broken-axis slash thickness


# =====================================================================
# OR Gate Parameters
# =====================================================================
omega_c_amp = 2 * np.pi * 50  # Control pulse amplitude [MHz, angular]
omega_p_amp = 2 * np.pi * 70  # Target probe pulse amplitude
omega_R_amp = 2.5 * omega_p_amp  # Target Rydberg coupling amplitude

T_c = np.pi / omega_c_amp  # Control pi-pulse duration: ~0.01 us (10 ns)
T_f = 0.32653  # Target pulse half-window [us]

or_args = {
    "omega_c_amp": omega_c_amp,
    "omega_p_amp": omega_p_amp,
    "omega_R_amp": omega_R_amp,
    "T_c": T_c,
    "T_f": T_f,
}

t_total_or = 2 * T_c + 2 * T_f

print("OR gate timings:")
print(f"  T_c      = {T_c * 1e3:.2f} ns  (control pi-pulse)")
print(f"  2*T_f    = {2 * T_f * 1e3:.2f} ns  (Hanning probe window)")
print(f"  t_total  = {t_total_or * 1e3:.2f} ns")
print(f"  Ratio    = 2*T_f / T_c = {2 * T_f / T_c:.1f}x\n")

# =====================================================================
# CCX Gate Parameters
# =====================================================================
omega_cc_amp = 2 * np.pi * 100  # Control pulse amplitude [MHz, angular]
omega_t_amp = 2 * np.pi * 50  # Target pulse amplitude

T_cc = np.pi / omega_cc_amp  # Control pi-pulse duration: ~0.005 us (5 ns)
T_t = np.pi / omega_t_amp  # Target sub-pulse duration: ~0.01 us (10 ns)

ccx_args = {
    "omega_cc_amp": omega_cc_amp,
    "omega_t_amp": omega_t_amp,
    "T_cc": T_cc,
    "T_t": T_t,
}

t_total_ccx = 2 * T_cc + 3 * T_t

print("CCX gate timings:")
print(f"  T_cc     = {T_cc * 1e3:.2f} ns  (control pi-pulse)")
print(f"  T_t      = {T_t * 1e3:.2f} ns  (target sub-pulse)")
print(f"  t_total  = {t_total_ccx * 1e3:.2f} ns\n")


# =====================================================================
# Helpers
# =====================================================================
def eval_pulse(func, tlist, args):
    """Evaluate a pulse function at every point in tlist."""
    return np.array([func(t, args) for t in tlist])


def mask_zero(vals, threshold=1e-10):
    """Replace zero-valued regions with NaN, but keep boundary zeros.

    Interior zeros (far from any non-zero sample) become NaN so
    matplotlib draws nothing there.  Boundary zeros -- those directly
    adjacent to a non-zero sample -- are kept so that square-pulse
    vertical edges are preserved.
    """
    out = vals.astype(float).copy()
    is_zero = np.abs(vals) < threshold
    n = len(vals)
    for i in range(n):
        if is_zero[i]:
            left_nz = (i > 0) and not is_zero[i - 1]
            right_nz = (i < n - 1) and not is_zero[i + 1]
            if not (left_nz or right_nz):
                out[i] = np.nan
    return out


# =====================================================================
# FIGURE 1: OR Gate -- Broken-Axis Pulse Sequence
# =====================================================================
# Strategy: evaluate every pulse over the FULL time range once, then
# plot the same data on all three axes.  Each panel clips via xlim,
# so pulses (especially the constant Omega_R) run unbroken to the
# panel edges -- no vertical artefacts, no gaps.

# -- Single dense time array spanning the entire gate --
N_full = 4000
t_full = np.linspace(0, t_total_or, N_full)
t_ns = t_full * 1e3  # nanoseconds for plotting

or_pulses = [
    (omega_c * 2, COLOR_CTRL, r"$\Omega_c$ (control)"),
    (omega_p * 2, COLOR_PROBE, r"$\Omega_p$ (probe)"),
    (omega_R * 2, COLOR_COUPLE, r"$\Omega_R$ (coupling)"),
]

# Pre-evaluate all pulses (in MHz, divided by 2pi), then mask zeros
or_vals = {}
for fn, _, lbl in or_pulses:
    raw = eval_pulse(fn, t_full, or_args) / (2 * np.pi)
    or_vals[lbl] = mask_zero(raw)

# -- Panel x-limits (ns) --
eps = 0.0005 * 1e3  # 0.5 ns -- smaller than any pulse feature
pad_c = 0.35 * T_c * 1e3  # ~3.5 ns visual padding after ctrl pulse

xlim_L = (0, T_c * 1e3 + pad_c)
xlim_M = (T_c * 1e3 + eps, (T_c + 2 * T_f) * 1e3 - eps)
xlim_R = ((T_c + 2 * T_f) * 1e3 - pad_c, t_total_or * 1e3)

# Width ratios -- inflate control segments so their ~10 ns pulses
# are visually comparable to the ~653 ns probe window.
width_ratios = [0.18, 0.64, 0.18]

fig1, (ax_L, ax_M, ax_R) = plt.subplots(
    1,
    3,
    sharey=True,
    figsize=(7, 5.25),
    gridspec_kw={"width_ratios": width_ratios, "wspace": 0.06},
)

# -- Plot the SAME full-range arrays on every panel --
for ax, xlim in [(ax_L, xlim_L), (ax_M, xlim_M), (ax_R, xlim_R)]:
    for _, color, label in or_pulses:
        ax.plot(t_ns, or_vals[label], color=color, linewidth=LW, label=label, zorder=3)
    ax.set_xlim(xlim)

# -- Tighten y-axis around data: pulses range from -25 to 87.5 MHz --
ax_L.set_ylim(-35, 100)

# -- Move x-axis (bottom spine) to y=0 on all panels --
for ax in (ax_L, ax_M, ax_R):
    ax.spines["bottom"].set_position(("data", 0))
    ax.spines["bottom"].set_linewidth(AXIS_LW)

# -- Axis labels --
ax_L.set_ylabel(r"Amplitude / 2$\pi$ (MHz)", fontsize=13)
fig1.text(0.53, 0.08, "Time (ns)", ha="center", fontsize=13)

# -- Spine & tick styling for broken axis --
ax_L.spines["top"].set_visible(False)
ax_L.spines["right"].set_visible(False)
ax_L.spines["left"].set_linewidth(AXIS_LW)
ax_L.tick_params(right=False, width=TICK_W)

for sp in ("top", "left", "right"):
    ax_M.spines[sp].set_visible(False)
ax_M.tick_params(left=False, right=False, width=TICK_W)

ax_R.spines["top"].set_visible(False)
ax_R.spines["left"].set_visible(False)
ax_R.tick_params(left=False, labelleft=False, width=TICK_W)

# -- Y-axis grid on all panels --
for ax in (ax_L, ax_M, ax_R):
    ax.grid(axis="y", linestyle="--", linewidth=0.5, color="#cccccc", zorder=0)
    ax.set_axisbelow(True)

# -- Broken-axis diagonal marks (centred on the x-axis at y=0) --
ylim = ax_L.get_ylim()
y0_frac = (0.0 - ylim[0]) / (ylim[1] - ylim[0])

d = 0.015  # mark size in axes fraction
bk = dict(color="k", clip_on=False, linewidth=BREAK_LW)

# Right edge of left panel (on x-axis only, no top marks)
ax_L.plot((1 - d, 1 + d), (y0_frac - d, y0_frac + d), transform=ax_L.transAxes, **bk)

# Both edges of middle panel
ax_M.plot((-d, +d), (y0_frac - d, y0_frac + d), transform=ax_M.transAxes, **bk)
ax_M.plot((1 - d, 1 + d), (y0_frac - d, y0_frac + d), transform=ax_M.transAxes, **bk)

# Left edge of right panel
ax_R.plot((-d, +d), (y0_frac - d, y0_frac + d), transform=ax_R.transAxes, **bk)

# -- Title --
fig1.suptitle(
    "OR Gate Pulse Sequence",
    fontweight="bold",
    fontsize=14,
    y=0.95,
)

# -- Legend (figure-level, below title) --
handles, labels = ax_L.get_legend_handles_labels()
fig1.legend(
    handles,
    labels,
    loc="upper center",
    ncol=3,
    bbox_to_anchor=(0.5, 0.92),
    frameon=False,
    fontsize=10,
    handlelength=1.5,
    columnspacing=1.0,
)

fig1.subplots_adjust(top=0.88, bottom=0.12, left=0.10, right=0.97)
fig1.savefig("results/OR/or_gate_pulse_sequence.pdf", bbox_inches="tight")
fig1.savefig("results/OR/or_gate_pulse_sequence.png", dpi=200, bbox_inches="tight")
print("Saved: results/OR/or_gate_pulse_sequence.{pdf,png}")


# =====================================================================
# FIGURE 2: CCX Gate -- Standard Pulse Sequence
# =====================================================================
tlist_ccx = np.linspace(0, t_total_ccx, 1000)

fig2, ax2 = plt.subplots(figsize=(7, 5.25), constrained_layout=False)

ccx_pulses = [
    (omega_cc * 2, r"$\Omega_{cc}$ (control)", COLOR_CTRL),
    (omega_t1 * 2, r"$\Omega_{t1}$ ($|B\rangle\!\leftrightarrow\!|R\rangle$)", COLOR_PROBE),
    (
        omega_t2 * 2,
        r"$\Omega_{t2}$ ($|A\rangle\!\leftrightarrow\!|R\rangle$)",
        COLOR_COUPLE,
    ),
]

for func, label, color in ccx_pulses:
    raw = eval_pulse(func, tlist_ccx, ccx_args) / (2 * np.pi)
    ax2.plot(
        tlist_ccx * 1e3,
        mask_zero(raw),
        linewidth=LW,
        label=label,
        color=color,
        zorder=3,
    )

ax2.set_ylabel(r"Amplitude / 2$\pi$ (MHz)", fontsize=13)
ax2.set_xlim(0, t_total_ccx * 1e3)

# Tighten y-axis: pulses range from -50 to +50 MHz
ax2.set_ylim(-60, 60)

# Move x-axis to y=0
ax2.spines["bottom"].set_position(("data", 0))
ax2.spines["bottom"].set_linewidth(AXIS_LW)
ax2.spines["left"].set_linewidth(AXIS_LW)
for spine in ("top", "right"):
    ax2.spines[spine].set_visible(False)
ax2.tick_params(axis="both", which="major", length=3, pad=3, width=TICK_W)
ax2.grid(axis="y", linestyle="--", linewidth=0.5, color="#cccccc", zorder=0)
ax2.set_axisbelow(True)

# xlabel close to the figure bottom
fig2.text(0.53, 0.08, "Time (ns)", ha="center", fontsize=13)

# Title
fig2.suptitle(
    "CCX Gate Pulse Sequence",
    fontweight="bold",
    fontsize=14,
    y=0.95,
)

# Figure-level legend, below title
handles2, labels2 = ax2.get_legend_handles_labels()
fig2.legend(
    handles2,
    labels2,
    loc="upper center",
    ncol=3,
    bbox_to_anchor=(0.53, 0.92),
    frameon=False,
    fontsize=10,
    handlelength=1.5,
    columnspacing=1.0,
)

fig2.subplots_adjust(top=0.88, bottom=0.12, left=0.12, right=0.97)
fig2.savefig("results/CCX/ccx_gate_pulse_sequence.pdf", bbox_inches="tight")
fig2.savefig("results/CCX/ccx_gate_pulse_sequence.png", dpi=200, bbox_inches="tight")
print("Saved: results/CCX/ccx_gate_pulse_sequence.{pdf,png}")

plt.close("all")
print("\nDone.")

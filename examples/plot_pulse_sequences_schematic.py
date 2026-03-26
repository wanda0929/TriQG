"""
Schematic Pulse Sequence Diagrams (not to scale)

Shows pulse shapes and types only — no real amplitudes or durations.
Saves:
  results/OR/or_gate_pulse_sequence_schematic.pdf/.png
  results/CCX/ccx_gate_pulse_sequence_schematic.pdf/.png

Run from repository root:
    PYTHONPATH=. python3 examples/plot_pulse_sequences_schematic.py
"""

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({"font.family": "sans-serif", "font.size": 11})

COLOR_CTRL   = "#5B8DBE"
COLOR_PROBE  = "#E76F51"
COLOR_COUPLE = "#7FB685"
LW    = 2.5
ALPHA = 0.18
N     = 3000


# ── helpers ──────────────────────────────────────────────────────────

def step_func(t, segments):
    """Piecewise-constant waveform.  segments = [(t0, t1, amplitude), ...]"""
    out = np.zeros_like(t, dtype=float)
    for t0, t1, a in segments:
        out[(t >= t0) & (t < t1)] = a
    return out


def super_gaussian(t, tc, hw, order=6):
    """Super-Gaussian of given order centred at tc, half-width hw."""
    sigma = hw / (np.log(100) ** (1.0 / order))
    out = np.zeros_like(t, dtype=float)
    mask = (t >= tc - hw) & (t <= tc + hw)
    out[mask] = np.exp(-((t[mask] - tc) / sigma) ** order)
    return out


def draw_channel(ax, t, vals, color, ylabel, ylim=(-1.5, 1.5)):
    ax.plot(t, vals, color=color, lw=LW)
    ax.fill_between(t, 0, vals, where=(vals > 0), alpha=ALPHA, color=color)
    ax.fill_between(t, 0, vals, where=(vals < 0), alpha=ALPHA, color=color)
    ax.axhline(0, color="k", lw=0.8, zorder=0)
    ax.set_xlim(t[0], t[-1])
    ax.set_ylim(*ylim)
    ax.set_yticks([])
    ax.set_ylabel(ylabel, rotation=0, ha="right", va="center",
                  fontsize=11, labelpad=6)
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)
    ax.tick_params(left=False, bottom=False)


def note(ax, x, y, text, color, va="bottom"):
    ax.text(x, y, text, ha="center", va=va, fontsize=9,
            color=color, style="italic")


def finalize_rows(axes):
    """Hide intermediate x-axis spines; keep only bottom one."""
    for ax in axes[:-1]:
        ax.spines["bottom"].set_visible(False)
    axes[-1].spines["bottom"].set_linewidth(0.8)
    axes[-1].set_xticks([])
    axes[-1].set_xlabel("Time (arb. u.)", fontsize=12)


GRID = dict(hspace=0.08, top=0.87, bottom=0.12, left=0.28, right=0.97)


# =====================================================================
# Figure 1 — OR Gate Schematic
# =====================================================================
#
# Schematic timeline (arbitrary units, not to scale):
#
#  [0, 1]   Ω_c  : +π square pulse      ─── control atoms |1⟩↔|r⟩
#  [1, 6]   Ω_g  : super-Gaussian probe  ─── target |A,B⟩↔|P⟩
#  [0, 7]   Ω_R  : constant coupling     ─── target |P⟩↔|R⟩
#  [6, 7]   Ω_c  : −π square pulse
#
t_or = np.linspace(0, 7, N)

oc_or = step_func(t_or, [(0, 1, 1.0), (6, 7, -1.0)])
og_or = super_gaussian(t_or, tc=3.5, hw=2.5)
oR_or = np.ones_like(t_or)

fig, axes = plt.subplots(3, 1, sharex=True, figsize=(6.5, 4.2),
                          gridspec_kw=GRID)

# row 0 — Ω_c
draw_channel(axes[0], t_or, oc_or, COLOR_CTRL,
             r"$\Omega_c$" + "\n" + r"$|1\rangle\!\leftrightarrow\!|r\rangle$")
note(axes[0], 0.5,  1.1,  r"$\pi$-pulse",  COLOR_CTRL)
note(axes[0], 6.5, -1.05, r"$-\pi$-pulse", COLOR_CTRL, va="top")

# row 1 — Ω_g
draw_channel(axes[1], t_or, og_or, COLOR_PROBE,
             r"$\Omega_g$" + "\n" + r"$|A,B\rangle\!\leftrightarrow\!|P\rangle$",
             ylim=(-0.3, 1.5))
note(axes[1], 3.5, 1.1, "super-Gaussian", COLOR_PROBE)

# row 2 — Ω_R
draw_channel(axes[2], t_or, oR_or, COLOR_COUPLE,
             r"$\Omega_R$" + "\n" + r"$|P\rangle\!\leftrightarrow\!|R\rangle$",
             ylim=(-0.3, 1.5))
note(axes[2], 3.5, 1.1, "constant", COLOR_COUPLE)

finalize_rows(axes)
fig.suptitle("OR Gate Pulse Sequence (Schematic)", fontweight="bold", fontsize=13)
fig.savefig("results/OR/or_gate_pulse_sequence_schematic.pdf", bbox_inches="tight")
fig.savefig("results/OR/or_gate_pulse_sequence_schematic.png", dpi=200, bbox_inches="tight")
print("Saved: results/OR/or_gate_pulse_sequence_schematic.{pdf,png}")
plt.close(fig)


# =====================================================================
# Figure 2 — CCX Gate Schematic
# =====================================================================
#
# Schematic timeline (arbitrary units, not to scale):
#
#  [0,   1.5]  Ω_cc : +π square pulse   ─── control atoms |0⟩↔|r⟩
#  [1.5, 3.5]  Ω_t1 : sub-pulse 1       ─── target |B⟩↔|R⟩
#  [3.5, 5.5]  Ω_t2 : sub-pulse 2       ─── target |A⟩↔|R⟩
#  [5.5, 7.5]  Ω_t1 : sub-pulse 3       ─── target |B⟩↔|R⟩
#  [7.5, 9.0]  Ω_cc : −π square pulse
#
t_ccx = np.linspace(0, 9, N)

occ = step_func(t_ccx, [(0, 1.5, 1.0), (7.5, 9.0, -1.0)])
ot1 = step_func(t_ccx, [(1.5, 3.5, 1.0), (5.5, 7.5, 1.0)])
ot2 = step_func(t_ccx, [(3.5, 5.5, 1.0)])

fig, axes = plt.subplots(3, 1, sharex=True, figsize=(6.5, 4.2),
                          gridspec_kw=GRID)

# row 0 — Ω_cc
draw_channel(axes[0], t_ccx, occ, COLOR_CTRL,
             r"$\Omega_{cc}$" + "\n" + r"$|0\rangle\!\leftrightarrow\!|r\rangle$")
note(axes[0], 0.75,  1.1,  r"$\pi$-pulse",  COLOR_CTRL)
note(axes[0], 8.25, -1.05, r"$-\pi$-pulse", COLOR_CTRL, va="top")

# row 1 — Ω_t1
draw_channel(axes[1], t_ccx, ot1, COLOR_PROBE,
             r"$\Omega_{t1}$" + "\n" + r"$|B\rangle\!\leftrightarrow\!|R\rangle$",
             ylim=(-0.3, 1.5))
note(axes[1], 2.5, 1.1, "sub-pulse 1", COLOR_PROBE)
note(axes[1], 6.5, 1.1, "sub-pulse 3", COLOR_PROBE)

# row 2 — Ω_t2
draw_channel(axes[2], t_ccx, ot2, COLOR_COUPLE,
             r"$\Omega_{t2}$" + "\n" + r"$|A\rangle\!\leftrightarrow\!|R\rangle$",
             ylim=(-0.3, 1.5))
note(axes[2], 4.5, 1.1, "sub-pulse 2", COLOR_COUPLE)

finalize_rows(axes)
fig.suptitle("CCX Gate Pulse Sequence (Schematic)", fontweight="bold", fontsize=13)
fig.savefig("results/CCX/ccx_gate_pulse_sequence_schematic.pdf", bbox_inches="tight")
fig.savefig("results/CCX/ccx_gate_pulse_sequence_schematic.png", dpi=200, bbox_inches="tight")
print("Saved: results/CCX/ccx_gate_pulse_sequence_schematic.{pdf,png}")
plt.close(fig)

print("\nDone.")

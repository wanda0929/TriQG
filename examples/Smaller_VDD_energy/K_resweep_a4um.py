"""
K-coefficient re-sweep at the brute-force champion (a = 4 um) and at
the rev. 4 FINAL cell (a = 5 um), with the explicit goal of testing
the *finite-blockade amplitude-rescaling prediction* from

    Theis, Motzoi, Wilhelm, Saffman, PRA 94, 032306 (2016)
    arXiv:1605.08891

namely that the optimal two-photon pulse area at finite V_ct sits
*below* the perfect-blockade nominal value pi/4 by a relative amount
~ 1 / M_2, where M_2 = V_ct (2 Delta) / (Omega_p Omega_R) is the
strong-blockade margin (see coefficients_trial.typ Sec. 1).

In our notation:

    K = area / (pi/4),
    K_opt(predicted) ~ 1 - 1 / M_2.

What this script does
---------------------
1. Sweeps K in a fine grid around K = 1 at five representative cells:
      A. rev. 4 FINAL              (a = 5.0,  Op = 50, r = 3.5, D = 500, Oc = 50)
      B. champion (BF rank 1)      (a = 4.0,  Op = 60, r = 4.0, D = 500, Oc = 70)
      C. fastest robust            (a = 4.0,  Op = 60, r = 4.0, D = 400, Oc = 70)
      D. a = 4 um, low-D variant   (a = 4.0,  Op = 50, r = 3.5, D = 500, Oc = 50)
      E. a = 4.5 um, ratio = 3.5   (a = 4.5,  Op = 60, r = 3.5, D = 500, Oc = 70)

   For each cell we compute F_bar at K in [0.88, ..., 1.02] (step 0.01).

2. Estimates K_opt for each cell by parabolic interpolation around the
   argmax of F_bar.

3. Compares K_opt(empirical) against the Theis-style prediction
      K_opt(predicted) = 1 - c / M_2,
   fitted across the five cells (linear regression through the origin
   in (1 - K_opt) vs 1/M_2 space).

Outputs
-------
   K_resweep_a4um.csv          -- every (cell, K) row, fidelities, M_2, etc.
   K_resweep_a4um_optimum.csv  -- one row per cell: K_opt (empirical),
                                   K_opt (predicted), residual, M_2
   K_resweep_a4um.png          -- two-panel plot:
                                   (left)  F_bar vs K for each cell
                                   (right) (1 - K_opt) vs 1/M_2 fit

Wall time
---------
   5 cells * 15 K-points * 8 mesolve runs ~ 600 mesolve calls.
   At ~0.15 s each on Apple Silicon -> ~ 90 s.
"""

from __future__ import annotations

import csv
import time
from dataclasses import dataclass

import numpy as np
import matplotlib.pyplot as plt

from triqg.atoms import CsAtom, RbAtom, composite_basis_state
from triqg.pulses import omega_gaussian, compute_pulse_area
from triqg.hamiltonian import build_hamiltonian
from triqg.decoherence import build_collapse_operators
from triqg.solver import simulate
from triqg.analysis import state_fidelity, average_gate_fidelity

from level_parameters import C3_tilde, C6_RbRb, C6_CsCs, gamma_r, gamma_R, gamma_P


# =====================================================================
# Configuration
# =====================================================================
@dataclass
class Cell:
    name: str
    a_um: float
    omega_p_MHz: float
    ratio: float
    delta_MHz: float
    omega_c_MHz: float


CELLS = [
    Cell("REV4_FINAL",    a_um=5.00, omega_p_MHz=50.0, ratio=3.5,
         delta_MHz=500.0, omega_c_MHz=50.0),
    Cell("CHAMPION",      a_um=4.00, omega_p_MHz=60.0, ratio=4.0,
         delta_MHz=500.0, omega_c_MHz=70.0),
    Cell("FASTEST",       a_um=4.00, omega_p_MHz=60.0, ratio=4.0,
         delta_MHz=400.0, omega_c_MHz=70.0),
    Cell("a4_lowOp",      a_um=4.00, omega_p_MHz=50.0, ratio=3.5,
         delta_MHz=500.0, omega_c_MHz=50.0),
    Cell("a45_r35",       a_um=4.50, omega_p_MHz=60.0, ratio=3.5,
         delta_MHz=500.0, omega_c_MHz=70.0),
]

# K-sweep grid (fine spacing near 1).  We extend a touch below the
# rev. 3 K = 0.95 to make sure we bracket the optimum at a = 4 um.
K_GRID = np.round(np.arange(0.88, 1.021, 0.01), 4)

ALPHA = 4.0
I_INF = 2 * 0.92770 * 2 ** (-1.0 / 6.0)   # ~ 1.6534

# =====================================================================
# Build basis once
# =====================================================================
cs = CsAtom()
rb = RbAtom()
ctrl_levels = [cs.level_index["0"], cs.level_index["1"]]
tgt_levels  = [rb.level_index["A"], rb.level_index["B"]]

basis_inputs = []
for c1 in range(2):
    for c2 in range(2):
        for t in range(2):
            psi_in = composite_basis_state(
                ctrl_levels[c1], ctrl_levels[c2], tgt_levels[t]
            )
            psi_id = composite_basis_state(
                ctrl_levels[c1], ctrl_levels[c2],
                tgt_levels[1 - t if (c1 or c2) else t],
            )
            label = f"|{c1}{c2}{['A','B'][t]}>"
            basis_inputs.append((label, c1, c2, t, psi_in, psi_id))

c_ops = build_collapse_operators(gamma_r, gamma_R, gamma_P)
solver_opts = {"store_final_state": True, "nsteps": 200000}


# =====================================================================
# Helpers
# =====================================================================
def cell_constants(cell: Cell):
    """Return (V_ct_MHz, V_cc_MHz, V_DD_MHz, omega_*, delta, T_c, M_2)."""
    r_DA = cell.a_um / np.sqrt(2)
    r_DD = cell.a_um
    r_AA = cell.a_um * np.sqrt(2)
    V_ct_MHz = 1000.0 * C3_tilde / r_DA ** 3
    V_DD_MHz = 1000.0 * C6_RbRb  / r_DD ** 6
    V_cc_MHz = 1000.0 * C6_CsCs  / r_AA ** 6
    V_ct = 2 * np.pi * V_ct_MHz
    V_cc = 2 * np.pi * V_cc_MHz

    omega_p_amp = 2 * np.pi * cell.omega_p_MHz
    omega_R_MHz = cell.ratio * cell.omega_p_MHz
    omega_R_amp = 2 * np.pi * omega_R_MHz
    omega_c_amp = 2 * np.pi * cell.omega_c_MHz
    delta       = 2 * np.pi * cell.delta_MHz
    T_c         = np.pi / omega_c_amp

    two_photon_R = cell.omega_p_MHz * omega_R_MHz / (2 * cell.delta_MHz)
    M_2 = V_ct_MHz / two_photon_R

    return dict(
        V_ct_MHz=V_ct_MHz, V_DD_MHz=V_DD_MHz, V_cc_MHz=V_cc_MHz,
        V_ct=V_ct, V_cc=V_cc,
        omega_p_amp=omega_p_amp, omega_R_amp=omega_R_amp,
        omega_R_MHz=omega_R_MHz, omega_c_amp=omega_c_amp,
        delta=delta, T_c=T_c, M_2=M_2,
    )


def run_K(cell: Cell, K: float, const: dict):
    """Return (F_bar, F_00, F_01, F_10, F_11, sigma, T_f, t_total, area)."""
    sigma_pi4 = ((2 * np.pi * const["delta"])
                 / (const["omega_p_amp"] ** 2 * I_INF)) ** 3
    sigma = (K ** 3) * sigma_pi4
    T_f   = (ALPHA * sigma) ** (1.0 / 3.0)
    T_c   = const["T_c"]
    t_total = 2 * T_c + 2 * T_f

    H = build_hamiltonian(
        const["delta"], const["V_ct"],
        pulse_p=omega_gaussian, V_cc=const["V_cc"],
    )
    args = {
        "omega_c_amp": const["omega_c_amp"],
        "omega_p_amp": const["omega_p_amp"],
        "omega_R_amp": const["omega_R_amp"],
        "T_c": T_c, "T_f": T_f, "sigma": sigma,
    }
    area = compute_pulse_area(
        omega_gaussian, const["delta"], T_c, T_c + 2 * T_f, args,
    )

    tlist = np.linspace(0, t_total, 400)
    pairs = []
    per_F = []
    for lbl, c1, c2, t, psi_in, psi_id in basis_inputs:
        res = simulate(
            method="mesolve", H=H, psi0=psi_in, tlist=tlist,
            c_ops=c_ops, e_ops=[], options=solver_opts, args=args,
        )
        pairs.append((res.final_state, psi_id))
        per_F.append(state_fidelity(res.final_state, psi_id))
    F_bar = average_gate_fidelity(pairs)
    F_00 = 0.5 * (per_F[0] + per_F[1])
    F_01 = 0.5 * (per_F[2] + per_F[3])
    F_10 = 0.5 * (per_F[4] + per_F[5])
    F_11 = 0.5 * (per_F[6] + per_F[7])
    return F_bar, F_00, F_01, F_10, F_11, sigma, T_f, t_total, area


def parabolic_argmax(K_vals: np.ndarray, F_vals: np.ndarray) -> float:
    """Three-point parabolic interpolation around the discrete argmax."""
    i = int(np.argmax(F_vals))
    if i == 0 or i == len(K_vals) - 1:
        return float(K_vals[i])
    x0, x1, x2 = K_vals[i - 1], K_vals[i], K_vals[i + 1]
    y0, y1, y2 = F_vals[i - 1], F_vals[i], F_vals[i + 1]
    denom = (y0 - 2 * y1 + y2)
    if abs(denom) < 1e-14:
        return float(x1)
    dx = 0.5 * (y0 - y2) / denom
    # dx in units of grid step; convert
    h = x1 - x0
    return float(x1 + dx * h)


# =====================================================================
# Sweep
# =====================================================================
print("=" * 88)
print("K-resweep at five (a, Omega_p, ratio, Delta, Omega_c) cells")
print("Testing the finite-blockade amplitude-rescaling prediction:")
print("                   K_opt ~ 1 - 1 / M_2")
print("=" * 88)
print(f"K grid: {K_GRID[0]} -> {K_GRID[-1]}  ({len(K_GRID)} points, step "
      f"{K_GRID[1]-K_GRID[0]:.3f})")
print(f"alpha = {ALPHA}\n")

rows = []
cell_summary = []
t_start = time.time()

for cell in CELLS:
    const = cell_constants(cell)
    print("-" * 88)
    print(f"Cell {cell.name}:  a = {cell.a_um:.2f} um,  "
          f"Omega_p = {cell.omega_p_MHz:.0f},  r = {cell.ratio},  "
          f"Delta = {cell.delta_MHz:.0f},  Omega_c = {cell.omega_c_MHz:.0f}")
    print(f"          V_ct = {const['V_ct_MHz']:.2f} MHz   "
          f"M_2 = {const['M_2']:.2f}   "
          f"K_pred (1 - 1/M_2) = {1.0 - 1.0/const['M_2']:.4f}")
    print(f"{'K':>6s} {'area/(pi/4)':>11s} {'sigma(ns)':>10s} {'T_f(ns)':>9s} "
          f"{'2Tf(ns)':>9s} {'Ttot(ns)':>9s} {'F_bar':>9s} "
          f"{'F00':>8s} {'F01':>8s} {'F11':>8s}")

    F_vals = []
    for K in K_GRID:
        F_bar, F00, F01, F10, F11, sigma, T_f, t_total, area = run_K(
            cell, float(K), const
        )
        F_vals.append(F_bar)
        rows.append({
            "cell": cell.name,
            "a_um": cell.a_um,
            "omega_p_MHz": cell.omega_p_MHz,
            "ratio": cell.ratio,
            "delta_MHz": cell.delta_MHz,
            "omega_c_MHz": cell.omega_c_MHz,
            "V_ct_MHz": const["V_ct_MHz"],
            "M_2": const["M_2"],
            "K": float(K),
            "area": area,
            "area_over_pi4": area / (np.pi / 4),
            "sigma_ns": sigma * 1e3,
            "T_f_ns":   T_f   * 1e3,
            "t_total_ns": t_total * 1e3,
            "F_bar": F_bar,
            "F_00": F00, "F_01": F01, "F_10": F10, "F_11": F11,
        })
        print(f"{K:6.3f} {area/(np.pi/4):11.5f} {sigma*1e3:10.4f} "
              f"{T_f*1e3:9.3f} {2*T_f*1e3:9.3f} {t_total*1e3:9.3f} "
              f"{F_bar:9.6f} {F00:8.5f} {F01:8.5f} {F11:8.5f}")

    F_vals = np.array(F_vals)
    K_opt = parabolic_argmax(K_GRID, F_vals)
    K_pred = 1.0 - 1.0 / const["M_2"]
    F_at_opt = float(F_vals.max())
    print(f"   -> K_opt(empirical) ~ {K_opt:.4f}   "
          f"K_opt(pred 1 - 1/M2) = {K_pred:.4f}   "
          f"diff = {K_opt - K_pred:+.4f}")
    cell_summary.append({
        "cell": cell.name,
        "a_um": cell.a_um,
        "omega_p_MHz": cell.omega_p_MHz,
        "ratio": cell.ratio,
        "delta_MHz": cell.delta_MHz,
        "omega_c_MHz": cell.omega_c_MHz,
        "V_ct_MHz": const["V_ct_MHz"],
        "M_2": const["M_2"],
        "inv_M2": 1.0 / const["M_2"],
        "K_opt_emp": K_opt,
        "K_opt_pred_1minus1overM2": K_pred,
        "K_opt_residual": K_opt - K_pred,
        "F_bar_at_K_opt": F_at_opt,
    })

# =====================================================================
# Linear regression of (1 - K_opt) on (1 / M_2), through the origin
# =====================================================================
inv_M2 = np.array([c["inv_M2"]      for c in cell_summary])
delta_K = np.array([1.0 - c["K_opt_emp"] for c in cell_summary])
# Slope via least-squares through origin:  c = sum(x*y) / sum(x^2)
slope = float(np.dot(inv_M2, delta_K) / np.dot(inv_M2, inv_M2))
# R^2 against the through-origin fit
y_fit = slope * inv_M2
ss_res = float(np.sum((delta_K - y_fit) ** 2))
ss_tot = float(np.sum(delta_K ** 2))   # variance about 0 (through-origin)
R2 = 1.0 - (ss_res / ss_tot if ss_tot > 0 else 0.0)

print("\n" + "=" * 88)
print("Theis-style scaling test:    1 - K_opt  =  c / M_2 + residual")
print("=" * 88)
print(f"{'cell':<14s} {'M_2':>6s} {'1/M_2':>8s} "
      f"{'1-Kopt':>9s} {'c/M_2 (fit)':>13s} {'resid':>9s}")
for c, dKi, fit in zip(cell_summary, delta_K, y_fit):
    print(f"{c['cell']:<14s} {c['M_2']:>6.2f} {c['inv_M2']:>8.5f} "
          f"{dKi:>9.5f} {fit:>13.5f} {dKi-fit:>9.5f}")
print(f"\nFitted slope c = {slope:.4f}   (Theis prediction: c ~ 1)")
print(f"R^2 of through-origin fit (variance about 0): {R2:.4f}")

# =====================================================================
# Persist
# =====================================================================
with open("K_resweep_a4um.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)
print("\nWrote K_resweep_a4um.csv "
      f"({len(rows)} rows = {len(CELLS)} cells x {len(K_GRID)} K values).")

with open("K_resweep_a4um_optimum.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(cell_summary[0].keys()))
    w.writeheader()
    w.writerows(cell_summary)
print("Wrote K_resweep_a4um_optimum.csv (one row per cell).")

# =====================================================================
# Plot
# =====================================================================
fig, axs = plt.subplots(1, 2, figsize=(13, 5.0))

# --- (a) F_bar vs K for each cell --------------------------------------
ax = axs[0]
colors = plt.get_cmap("tab10").colors
for i, cell in enumerate(CELLS):
    cs_rows = [r for r in rows if r["cell"] == cell.name]
    Ks = np.array([r["K"] for r in cs_rows])
    Fs = np.array([r["F_bar"] for r in cs_rows])
    c = colors[i % len(colors)]
    ax.plot(Ks, Fs, marker="o", lw=1.6, color=c,
            label=f"{cell.name} (M$_2$={cell_summary[i]['M_2']:.1f})")
    K_opt = cell_summary[i]["K_opt_emp"]
    K_pred = cell_summary[i]["K_opt_pred_1minus1overM2"]
    ax.axvline(K_opt, color=c, ls=":", alpha=0.5)
    ax.axvline(K_pred, color=c, ls="--", alpha=0.25)
ax.axhline(0.99, color="grey", lw=0.7, ls=":")
ax.axvline(1.0, color="black", lw=0.6)
ax.set_xlabel(r"$K = \mathrm{area} / (\pi/4)$")
ax.set_ylabel(r"$\overline{F}_{\mathrm{OR}}$")
ax.set_title(
    "K-resweep: F_bar vs pulse-area scale\n"
    "(dotted = K_opt empirical, dashed faint = 1 - 1/M_2 prediction)"
)
ax.legend(loc="lower center", fontsize=8, frameon=True)

# --- (b) (1 - K_opt) vs 1/M_2 ------------------------------------------
ax = axs[1]
for i, c in enumerate(cell_summary):
    ax.scatter(c["inv_M2"], 1.0 - c["K_opt_emp"],
               s=80, color=colors[i % len(colors)],
               label=c["cell"], zorder=3)
xs = np.linspace(0.0, max(inv_M2) * 1.15, 100)
ax.plot(xs, slope * xs, color="black", lw=1.2,
        label=f"fit: $1-K_{{\\mathrm{{opt}}}} = {slope:.2f} / M_2$  "
              f"($R^2 = {R2:.2f}$)")
ax.plot(xs, 1.0 * xs, color="grey", lw=1.0, ls="--",
        label=r"Theis prediction: $1 - K_{\mathrm{opt}} = 1 / M_2$")
ax.set_xlabel(r"$1 / M_2 \;=\; \Omega_p \Omega_R / [V_{ct} \cdot 2\Delta]$")
ax.set_ylabel(r"$1 - K_{\mathrm{opt}}$")
ax.set_title("Finite-blockade amplitude rescaling vs 1 / M$_2$")
ax.set_xlim(0.0, max(inv_M2) * 1.15)
ax.grid(True, alpha=0.3)
ax.legend(loc="upper left", fontsize=8, frameon=True)

plt.tight_layout()
plt.savefig("K_resweep_a4um.png", dpi=160)
print("Wrote K_resweep_a4um.png")

print(f"\nTotal wall time: {time.time() - t_start:.1f} s")

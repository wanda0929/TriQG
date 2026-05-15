"""
CCX (Omega_cc, Omega_t) sweep at a = 4 um, K = 1, ARC lifetimes
================================================================

Why this exists
---------------
`a4_K1_full_fidelity_analysis.py` evaluated the CCX at the rev. 4
amplitudes (Omega_cc = 50 MHz, Omega_t = 20 MHz) -- but those were
chosen back when the lattice was a = 5 um, V_ct/(2 pi) = 226 MHz.

At a = 4 um the interspecies blockade is V_ct/(2 pi) = 442 MHz --
roughly 2x larger -- so the old amplitudes sit deeper in the blockade
window than necessary (V_ct/Omega_t = 22) and pay a 95 ns gate-time
penalty.  At the new ARC lifetimes (tau_R = 83.4 us, tau_r = 83.9 us)
each ns of gate time costs ~1.2e-5 in survival; shaving 50 ns can
recover ~6e-4 of fidelity if leakage does not blow up.

This script sweeps over a coarse 4 x 4 grid of (Omega_cc, Omega_t),
finds the best F_bar_basis, then runs the full Pedersen + phase-
corrected audit on the optimum so we can compare F_PC directly against
the 0.99790 baseline.

Grid
----
    Omega_cc / (2 pi) MHz  in  {30, 50, 70, 100}
    Omega_t  / (2 pi) MHz  in  {20, 30, 40, 50}

Cost
----
    16 basis-only sweeps   (8 mesolve runs each, gate < 100 ns)  ~  3 min
    1 full Pedersen audit  (64 mesolve runs + propagator)        ~  1 min

Output
------
    ccx_omega_sweep_a4.csv     all 16 cells
    ccx_omega_sweep_a4.log     full stdout
    ccx_omega_sweep_a4.json    best-cell Pedersen + phase-correction
"""

from __future__ import annotations

import csv
import json
import time
from pathlib import Path

import numpy as np
import qutip
from scipy.optimize import minimize

from triqg.atoms import CsAtom, RbAtom, composite_basis_state
from triqg.hamiltonian import build_ccx_hamiltonian
from triqg.decoherence import build_collapse_operators
from triqg.pedersen import (
    build_permutation_unitary,
    choi_on_subspace_via_mesolve,
    pedersen_average_fidelity_from_choi,
)
from triqg.analysis import state_fidelity
from triqg.solver import simulate

from level_parameters import C3_tilde, C6_RbRb, C6_CsCs

# =====================================================================
# Physical setup at a = 4 um (same as a4_K1_full_fidelity_analysis.py)
# =====================================================================
a_um = 4.00
r_DA = a_um / np.sqrt(2)
r_DD = a_um
r_AA = a_um * np.sqrt(2)
V_ct_MHz = 1000.0 * C3_tilde / r_DA ** 3   # = 441.94 MHz
V_DD_MHz = 1000.0 * C6_RbRb  / r_DD ** 6   # = 2.214 MHz
V_cc_MHz = 1000.0 * C6_CsCs  / r_AA ** 6   # = -2.778 MHz
V_ct = 2 * np.pi * V_ct_MHz
V_cc = 2 * np.pi * V_cc_MHz

# ARC-verified lifetimes
tau_R_us = 83.41
tau_r_us = 83.90
tau_P_us = 0.131
gamma_r = 1.0 / tau_r_us
gamma_R = 1.0 / tau_R_us
gamma_P = 1.0 / tau_P_us
c_ops = build_collapse_operators(gamma_r, gamma_R, gamma_P)

# 8 computational basis kets
cs = CsAtom()
rb = RbAtom()
ctrl = [cs.level_index["0"], cs.level_index["1"]]
tgt  = [rb.level_index["A"], rb.level_index["B"]]
comp_kets = []
labels    = []
for c1 in range(2):
    for c2 in range(2):
        for t in range(2):
            comp_kets.append(composite_basis_state(ctrl[c1], ctrl[c2], tgt[t]))
            labels.append(f"|{c1},{c2},{['A','B'][t]}>")
V_iso = np.array([k.full().flatten() for k in comp_kets]).T


def ccx_perm(idx: int) -> int:
    c1 = (idx >> 2) & 1
    c2 = (idx >> 1) & 1
    t  = idx & 1
    return (c1 << 2) | (c2 << 1) | (t ^ (c1 & c2))


# =====================================================================
# Per-cell evaluator: 8 mesolve runs, return F_bar_basis + T_tot
# =====================================================================
H_ccx = build_ccx_hamiltonian(V_ct, V_cc=V_cc)


def eval_basis(omega_cc_MHz: float, omega_t_MHz: float) -> dict:
    omega_cc_amp = 2 * np.pi * omega_cc_MHz
    omega_t_amp  = 2 * np.pi * omega_t_MHz
    T_cc = np.pi / omega_cc_amp
    T_t  = np.pi / omega_t_amp
    t_total = 2 * T_cc + 3 * T_t
    args = {
        "omega_cc_amp": omega_cc_amp,
        "omega_t_amp":  omega_t_amp,
        "T_cc": T_cc, "T_t": T_t,
    }
    max_freq = max(omega_cc_amp, omega_t_amp) / (2 * np.pi)
    solver_opts = {
        "store_final_state": True, "nsteps": 200000,
        "max_step": 1.0 / (20 * max_freq),
        "atol": 1e-10, "rtol": 1e-8,
    }

    F_k = np.zeros(8)
    for k in range(8):
        psi_in = comp_kets[k]
        psi_id = comp_kets[ccx_perm(k)]
        tlist  = np.linspace(0, t_total, 350)
        res = simulate(
            method="mesolve", H=H_ccx, psi0=psi_in, tlist=tlist,
            c_ops=c_ops, e_ops=[], options=solver_opts, args=args,
        )
        F_k[k] = state_fidelity(res.final_state, psi_id)
    return {
        "omega_cc_MHz": omega_cc_MHz,
        "omega_t_MHz":  omega_t_MHz,
        "T_tot_ns":     1e3 * t_total,
        "F_bar_basis":  float(F_k.mean()),
        "F_k":          F_k.tolist(),
        "V_ct_over_Omega_t":  V_ct_MHz / omega_t_MHz,
        "V_ct_over_Omega_cc": V_ct_MHz / omega_cc_MHz,
    }


# =====================================================================
# Pedersen + phase correction (on top cell only)
# =====================================================================
def _F_pro_with_phase(choi_tensor, U0, phis):
    d = U0.shape[0]
    D = np.diag(np.exp(1j * phis))
    U_target = D @ U0
    F = np.einsum("ki,ijkl,lj->", U_target.conj(), choi_tensor, U_target) / (d * d)
    return float(np.real(F))


def _optimize_phase(choi_tensor, U0, n_restarts=20, seed=0):
    d = U0.shape[0]
    rng = np.random.default_rng(seed)
    best_phis = np.zeros(d)
    best_F = _F_pro_with_phase(choi_tensor, U0, best_phis)
    for _ in range(n_restarts):
        x0 = rng.uniform(-np.pi, np.pi, size=d)
        res = minimize(
            lambda phis: -_F_pro_with_phase(choi_tensor, U0, phis),
            x0, method="L-BFGS-B",
        )
        if -res.fun > best_F:
            best_F = -res.fun
            best_phis = res.x
    F_bar_pc = (d * best_F + 1) / (d + 1)
    return best_F, F_bar_pc, best_phis


def full_audit(omega_cc_MHz: float, omega_t_MHz: float) -> dict:
    omega_cc_amp = 2 * np.pi * omega_cc_MHz
    omega_t_amp  = 2 * np.pi * omega_t_MHz
    T_cc = np.pi / omega_cc_amp
    T_t  = np.pi / omega_t_amp
    t_total = 2 * T_cc + 3 * T_t
    args = {
        "omega_cc_amp": omega_cc_amp,
        "omega_t_amp":  omega_t_amp,
        "T_cc": T_cc, "T_t": T_t,
    }
    max_freq = max(omega_cc_amp, omega_t_amp) / (2 * np.pi)
    solver_opts = {
        "store_final_state": True, "nsteps": 200000,
        "max_step": 1.0 / (20 * max_freq),
        "atol": 1e-10, "rtol": 1e-8,
    }
    U0 = build_permutation_unitary(comp_kets, ccx_perm)

    # Basis fidelities
    F_k = np.zeros(8)
    for k in range(8):
        psi_in = comp_kets[k]
        psi_id = comp_kets[ccx_perm(k)]
        tlist  = np.linspace(0, t_total, 350)
        res = simulate(
            method="mesolve", H=H_ccx, psi0=psi_in, tlist=tlist,
            c_ops=c_ops, e_ops=[], options=solver_opts, args=args,
        )
        F_k[k] = state_fidelity(res.final_state, psi_id)

    # Choi tensor, Pedersen
    choi = choi_on_subspace_via_mesolve(
        H_ccx, c_ops, t_total, comp_kets,
        args=args, options=solver_opts, verbose=False,
    )
    F_bar_ped, F_pro, survival = pedersen_average_fidelity_from_choi(choi, U0)

    # Noiseless unitary
    U_full = qutip.propagator(
        H_ccx, t_total, c_ops=[], args=args,
        options={"atol": 1e-12, "rtol": 1e-10, "nsteps": 200000},
    )
    U_eff = V_iso.conj().T @ U_full.full() @ V_iso
    tr = np.trace(U0.conj().T @ U_eff)
    F_pro_u = abs(tr) ** 2 / 64.0
    F_bar_u = (8 * F_pro_u + 1) / 9.0

    # Phase correction
    F_pro_pc, F_bar_pc, best_phis = _optimize_phase(choi, U0)

    return {
        "omega_cc_MHz": omega_cc_MHz,
        "omega_t_MHz":  omega_t_MHz,
        "T_tot_ns":     1e3 * t_total,
        "F_k_basis_fidelities": F_k.tolist(),
        "F_bar_basis":   float(F_k.mean()),
        "F_pro_Pedersen_raw":            F_pro,
        "F_bar_Pedersen_raw":            F_bar_ped,
        "survival":                      survival,
        "F_pro_unitary":                 F_pro_u,
        "F_bar_unitary":                 F_bar_u,
        "F_pro_Pedersen_phase_corrected": F_pro_pc,
        "F_bar_Pedersen_phase_corrected": F_bar_pc,
        "optimal_phase_correction_over_pi": (best_phis / np.pi).tolist(),
    }


# =====================================================================
# Run sweep
# =====================================================================
GRID_OMEGA_CC = [30.0, 50.0, 70.0, 100.0]
GRID_OMEGA_T  = [20.0, 30.0, 40.0, 50.0]

print("=" * 78)
print("  CCX (Omega_cc, Omega_t) sweep  at  a = 4 um  with ARC lifetimes")
print("=" * 78)
print(f"V_ct / (2 pi) = +{V_ct_MHz:7.3f} MHz")
print(f"V_cc / (2 pi) = {V_cc_MHz:+8.4f} MHz")
print(f"Lifetimes (ARC):  tau_R = {tau_R_us:.2f} us,  "
      f"tau_r = {tau_r_us:.2f} us,  tau_P = {tau_P_us:.3f} us")
print()
print(f"Grid:  Omega_cc in {GRID_OMEGA_CC}  MHz")
print(f"       Omega_t  in {GRID_OMEGA_T}  MHz")
print()
print(f"{'Omega_cc':>9} {'Omega_t':>9} {'V_ct/Ot':>9} {'T_tot ns':>10} "
      f"{'F_bar_basis':>14}")
print("-" * 78)

t_start = time.time()
cells = []
for oc in GRID_OMEGA_CC:
    for ot in GRID_OMEGA_T:
        t0 = time.time()
        r = eval_basis(oc, ot)
        elapsed = time.time() - t0
        cells.append(r)
        print(
            f"{oc:>9.1f} {ot:>9.1f} {r['V_ct_over_Omega_t']:>9.2f} "
            f"{r['T_tot_ns']:>10.2f} {r['F_bar_basis']:>14.6f}   "
            f"({elapsed:.1f} s)"
        )

# Save full grid as CSV
here = Path(__file__).parent
csv_path = here / "ccx_omega_sweep_a4.csv"
with open(csv_path, "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["omega_cc_MHz", "omega_t_MHz", "T_tot_ns",
                "V_ct_over_Omega_t", "V_ct_over_Omega_cc", "F_bar_basis"])
    for c in cells:
        w.writerow([c["omega_cc_MHz"], c["omega_t_MHz"], c["T_tot_ns"],
                    c["V_ct_over_Omega_t"], c["V_ct_over_Omega_cc"],
                    c["F_bar_basis"]])
print()
print(f"Wrote {csv_path}")

# Top 3 by F_bar_basis
ranked = sorted(cells, key=lambda c: -c["F_bar_basis"])
print()
print("Top 3 cells by F_bar_basis:")
print(f"  {'rank':>4} {'Omega_cc':>9} {'Omega_t':>9} {'T_tot ns':>10} "
      f"{'F_bar_basis':>14}")
for i, c in enumerate(ranked[:3], 1):
    print(f"  {i:>4} {c['omega_cc_MHz']:>9.1f} {c['omega_t_MHz']:>9.1f} "
          f"{c['T_tot_ns']:>10.2f} {c['F_bar_basis']:>14.6f}")

# Full Pedersen audit on the top cell
print()
print("=" * 78)
print("Full Pedersen + phase-correction audit on top cell:")
print("=" * 78)
best = ranked[0]
t0 = time.time()
audit = full_audit(best["omega_cc_MHz"], best["omega_t_MHz"])
print(f"  audit completed in {time.time()-t0:.1f} s")
print()
print(f"  Omega_cc = {audit['omega_cc_MHz']:.1f} MHz, "
      f"Omega_t = {audit['omega_t_MHz']:.1f} MHz, "
      f"T_tot = {audit['T_tot_ns']:.2f} ns")
print()
print(f"  F_bar_basis (mean F_k)              = {audit['F_bar_basis']:.6f}")
print(f"  F_bar_Pedersen_raw                  = {audit['F_bar_Pedersen_raw']:.6f}")
print(f"  F_bar_Pedersen_phase_corrected      = "
      f"{audit['F_bar_Pedersen_phase_corrected']:.6f}")
print(f"  F_bar_unitary  (no decoherence)     = {audit['F_bar_unitary']:.6f}")
print(f"  Computational survival              = {audit['survival']:.6f}")

# Compare against rev. 4 a4 baseline (50, 20)
print()
print("Baseline comparison (Omega_cc=50, Omega_t=20 -- existing JSON):")
print(f"  Baseline F_bar_Pedersen_phase_corrected = 0.997897")
gain = audit["F_bar_Pedersen_phase_corrected"] - 0.997897
print(f"  Top-cell - baseline                     = {gain:+.6f}  "
      f"({1e3*gain:+.3f} x 10^-3)")

# Persist JSON for the best cell
json_path = here / "ccx_omega_sweep_a4.json"
out = {
    "physical_params": {
        "a_um": a_um,
        "V_ct_over_2pi_MHz": V_ct_MHz,
        "V_cc_over_2pi_MHz": V_cc_MHz,
        "tau_r_us_ARC": tau_r_us,
        "tau_R_us_ARC": tau_R_us,
        "tau_P_us":     tau_P_us,
    },
    "grid": {
        "omega_cc_MHz": GRID_OMEGA_CC,
        "omega_t_MHz":  GRID_OMEGA_T,
    },
    "all_cells_basis": cells,
    "best_cell_full_audit": audit,
}
json_path.write_text(json.dumps(out, indent=2))
print()
print(f"Wrote {json_path}")
print()
print(f"Total wall time: {time.time()-t_start:.1f} s")

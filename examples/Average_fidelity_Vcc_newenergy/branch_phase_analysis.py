"""
Branch-phase analysis for the OR and CCX gates.
================================================

Extracts the residual branch-dependent phases left after the nominal
+/- amplitude pulse design, and computes three fidelity metrics for each
gate:

    F_bar_basis         = (1/d) sum_k F_k                (Yu et al. Eq. 7)
    F_bar_Pedersen      = (d F_pro + 1)/(d+1)            (Pedersen Eq. 5, raw)
    F_bar_Pedersen_PC   = F_bar after optimal diagonal   (virtual-Z absorbed)
                          phase absorption into U_0

Also dumps the 8 branch phases phi_k = arg<k|U_eff U_0^dag|k> computed
from the noiseless unitary propagator, so we can see exactly which
V_cc / Raman / 3-pi-pulse contributions the +/- amplitude design fails
to cancel.

Writes everything to ``pedersen_fidelity_analysis.json`` in this
directory for consumption by the reference note.
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
import qutip
from scipy.optimize import minimize

from triqg.atoms import CsAtom, RbAtom, DIMS, composite_basis_state
from triqg.pulses import omega_gaussian
from triqg.hamiltonian import build_hamiltonian, build_ccx_hamiltonian
from triqg.decoherence import build_collapse_operators
from triqg.pedersen import (
    build_permutation_unitary,
    choi_on_subspace_via_mesolve,
    pedersen_average_fidelity_from_choi,
)

# =============================================================================
# Shared physical parameters
# =============================================================================
a_um = 5.00
r_DA = a_um / np.sqrt(2)
r_AA = a_um * np.sqrt(2)
C3_tilde = 22.84
C6_CsCs = -692.9
V_ct_MHz = 1000.0 * C3_tilde / r_DA**3
V_cc_MHz = 1000.0 * C6_CsCs / r_AA**6
V_ct = 2 * np.pi * V_ct_MHz
V_cc = 2 * np.pi * V_cc_MHz

gamma_r = 1.0 / 142.73
gamma_R = 1.0 / 134.87
gamma_P = 1.0 / 0.131
c_ops = build_collapse_operators(gamma_r, gamma_R, gamma_P)

# Computational-subspace basis (index: c1*4 + c2*2 + t)
cs = CsAtom(); rb = RbAtom()
ctrl = [cs.level_index["0"], cs.level_index["1"]]
tgt = [rb.level_index["A"], rb.level_index["B"]]
comp_kets = []
labels = []
for c1 in range(2):
    for c2 in range(2):
        for t in range(2):
            comp_kets.append(composite_basis_state(ctrl[c1], ctrl[c2], tgt[t]))
            labels.append(f"|{c1},{c2},{['A','B'][t]}>")
V = np.array([k.full().flatten() for k in comp_kets]).T  # 36 x 8

def or_perm(idx):
    c1 = (idx >> 2) & 1; c2 = (idx >> 1) & 1; t = idx & 1
    return (c1 << 2) | (c2 << 1) | (t ^ (c1 | c2))

def ccx_perm(idx):
    c1 = (idx >> 2) & 1; c2 = (idx >> 1) & 1; t = idx & 1
    return (c1 << 2) | (c2 << 1) | (t ^ (c1 & c2))


# =============================================================================
# Helper: process fidelity with a diagonal phase correction
# =============================================================================
def _process_fidelity_phase_corrected(choi_tensor: np.ndarray, U0: np.ndarray,
                                      phis: np.ndarray) -> float:
    """F_pro between E (given by its subspace Choi tensor) and (D U_0)
    where D = diag(exp(i phi_0), ..., exp(i phi_{d-1}))."""
    d = U0.shape[0]
    D = np.diag(np.exp(1j * phis))
    U_target = D @ U0
    F = np.einsum("ki,ijkl,lj->", U_target.conj(), choi_tensor, U_target) / (d * d)
    return float(np.real(F))


def _optimize_phase_correction(choi_tensor: np.ndarray, U0: np.ndarray,
                               n_restarts: int = 20, seed: int = 0):
    """Find the best diagonal-D phase correction maximizing Pedersen F_bar."""
    d = U0.shape[0]
    rng = np.random.default_rng(seed)
    best_phis = np.zeros(d)
    best_F_pro = _process_fidelity_phase_corrected(choi_tensor, U0, best_phis)

    for r in range(n_restarts):
        x0 = rng.uniform(-np.pi, np.pi, size=d)
        res = minimize(
            lambda phis: -_process_fidelity_phase_corrected(choi_tensor, U0, phis),
            x0, method="L-BFGS-B",
        )
        F_pro = -res.fun
        if F_pro > best_F_pro:
            best_F_pro = F_pro
            best_phis = res.x
    best_F_bar = (d * best_F_pro + 1) / (d + 1)
    return best_F_pro, best_F_bar, best_phis


# =============================================================================
# Per-gate analysis
# =============================================================================
def analyze_gate(name: str, H, args, t_total, gate_perm, solver_opts):
    print(f"\n======  {name}  ======")
    U0 = build_permutation_unitary(comp_kets, gate_perm)

    # ---------- Noiseless unitary propagator ----------
    print(" (1) noiseless unitary propagator ...", flush=True)
    t0 = time.time()
    U_full = qutip.propagator(
        H, t_total, c_ops=[], args=args,
        options={"atol": 1e-12, "rtol": 1e-10, "nsteps": 200000},
    )
    print(f"     done in {time.time()-t0:.1f} s")
    U_eff = V.conj().T @ U_full.full() @ V           # 8 x 8 projected

    tr = np.trace(U0.conj().T @ U_eff)
    F_pro_u = abs(tr) ** 2 / 64.0
    F_bar_u = (8 * F_pro_u + 1) / 9.0
    print(f"     |Tr(U0^dag U_eff)|^2 / d^2 = {F_pro_u:.6f}")
    print(f"     noiseless F_bar_Pedersen  = {F_bar_u:.6f}")

    # Branch phases and magnitudes on the diagonal of U_eff U_0^dag
    D = U_eff @ U0.conj().T
    diag = np.diag(D)
    branch_magnitudes = np.abs(diag).tolist()
    branch_phases_pi = (np.angle(diag) / np.pi).tolist()

    print("     branch magnitudes and phases (phases in pi-units):")
    for lbl, m, p in zip(labels, branch_magnitudes, branch_phases_pi):
        print(f"        {lbl:<12}  |d|={m:.4f}, phi={p:+.4f} pi")

    # ---------- Decoherent Choi via 64 mesolve runs ----------
    print(" (2) decoherent 64-run Choi tensor ...", flush=True)
    t0 = time.time()
    choi = choi_on_subspace_via_mesolve(
        H, c_ops, t_total, comp_kets, args=args, options=solver_opts,
        verbose=False,
    )
    print(f"     done in {time.time()-t0:.1f} s")

    F_bar_ped, F_pro, survival = pedersen_average_fidelity_from_choi(choi, U0)
    print(f"     survival T_P          = {survival:.6f}")
    print(f"     F_pro (raw)           = {F_pro:.6f}")
    print(f"     F_bar_Pedersen (raw)  = {F_bar_ped:.6f}")

    # ---------- Per-input basis fidelity from Choi diagonal ----------
    F_k = np.array([
        np.real(choi[i, i, gate_perm(i), gate_perm(i)]) for i in range(8)
    ])
    F_bar_basis = float(F_k.mean())
    print(f"     F_bar_basis (Yu eq.7) = {F_bar_basis:.6f}")

    # ---------- Phase-corrected Pedersen ----------
    print(" (3) phase-corrected (virtual-Z) Pedersen ...", flush=True)
    t0 = time.time()
    F_pro_pc, F_bar_pc, best_phis = _optimize_phase_correction(choi, U0)
    print(f"     done in {time.time()-t0:.2f} s")
    print(f"     F_pro (phase-corr.)   = {F_pro_pc:.6f}")
    print(f"     F_bar_PC              = {F_bar_pc:.6f}")
    print(f"     optimal phases / pi   = {[f'{p/np.pi:+.4f}' for p in best_phis]}")

    return {
        "name": name,
        "t_total_ns": 1e3 * t_total,
        "labels": labels,
        "branch_magnitudes": branch_magnitudes,
        "branch_phases_pi": branch_phases_pi,
        "F_bar_noiseless_Pedersen": F_bar_u,
        "F_pro_noiseless": F_pro_u,
        "F_bar_basis": F_bar_basis,
        "F_k_basis": F_k.tolist(),
        "F_bar_Pedersen_raw": F_bar_ped,
        "F_pro_raw": F_pro,
        "survival": survival,
        "F_bar_Pedersen_phase_corrected": F_bar_pc,
        "F_pro_phase_corrected": F_pro_pc,
        "optimal_phase_correction_over_pi": (np.array(best_phis) / np.pi).tolist(),
    }


# =============================================================================
# Run the two gates
# =============================================================================
# OR
omega_c_amp = 2 * np.pi * 50
omega_p_amp = 2 * np.pi * 50
omega_R_amp = 3.5 * omega_p_amp
delta = 2 * np.pi * 500
T_c = np.pi / omega_c_amp
T_f = 0.15
sigma = 0.001771
args_or = dict(omega_c_amp=omega_c_amp, omega_p_amp=omega_p_amp,
               omega_R_amp=omega_R_amp, T_c=T_c, T_f=T_f, sigma=sigma)
H_or = build_hamiltonian(delta, V_ct, pulse_p=omega_gaussian, V_cc=V_cc)
t_total_or = 2 * T_c + 2 * T_f
solver_opts_or = {"nsteps": 200000, "atol": 1e-10, "rtol": 1e-8}

# CCX
omega_cc_amp = 2 * np.pi * 100
omega_t_amp = 2 * np.pi * 50
T_cc = np.pi / omega_cc_amp
T_t = np.pi / omega_t_amp
args_ccx = dict(omega_cc_amp=omega_cc_amp, omega_t_amp=omega_t_amp,
                T_cc=T_cc, T_t=T_t)
H_ccx = build_ccx_hamiltonian(V_ct, V_cc=V_cc)
t_total_ccx = 2 * T_cc + 3 * T_t
max_freq = max(omega_cc_amp, omega_t_amp) / (2 * np.pi)
solver_opts_ccx = {"nsteps": 200000, "max_step": 1.0 / (20 * max_freq),
                   "atol": 1e-10, "rtol": 1e-8}

results = {
    "physical_params": {
        "a_um": a_um,
        "r_DA_um": r_DA,
        "r_AA_um": r_AA,
        "C3_tilde_GHz_um3": C3_tilde,
        "C6_CsCs_GHz_um6": C6_CsCs,
        "V_ct_over_2pi_MHz": V_ct_MHz,
        "V_cc_over_2pi_MHz": V_cc_MHz,
        "tau_r_us": 1.0 / gamma_r,
        "tau_R_us": 1.0 / gamma_R,
        "tau_P_us": 1.0 / gamma_P,
        "Omega_c_over_2pi_MHz": omega_c_amp / (2 * np.pi),
        "Omega_p_over_2pi_MHz": omega_p_amp / (2 * np.pi),
        "Omega_R_over_2pi_MHz": omega_R_amp / (2 * np.pi),
        "Omega_cc_over_2pi_MHz": omega_cc_amp / (2 * np.pi),
        "Omega_t_over_2pi_MHz": omega_t_amp / (2 * np.pi),
        "delta_over_2pi_MHz": delta / (2 * np.pi),
    },
    "OR": analyze_gate("OR gate (Option A, a=5 um, super-Gaussian)", H_or,
                      args_or, t_total_or, or_perm, solver_opts_or),
    "CCX": analyze_gate("CCX gate (Option A, a=5 um)", H_ccx,
                       args_ccx, t_total_ccx, ccx_perm, solver_opts_ccx),
}

out = Path(__file__).parent / "pedersen_fidelity_analysis.json"
out.write_text(json.dumps(results, indent=2))
print(f"\nWrote: {out}")

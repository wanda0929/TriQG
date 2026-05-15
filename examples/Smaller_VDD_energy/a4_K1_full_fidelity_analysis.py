"""
Full fidelity analysis at a = 4 um, K = 1 -- OR and CCX gates.

For each gate we report:
    F_k                      8 per-input computational-basis state fidelities
    F_bar_basis              (1/d) sum_k F_k                        (Yu et al. average)
    F_bar_Pedersen_raw       Pedersen Eq. (5) Haar-averaged fidelity     (Nielsen / Horodecki form)
    F_bar_Pedersen_PC        Pedersen after optimal diagonal phase    (virtual-Z absorbed)
                             correction  D = diag(exp(i phi_k))
    F_bar_unitary            Pedersen on the noiseless unitary       (separates coherent err)

Lifetimes (ARC 3.10.2, T = 0 K, spontaneous emission only;
           includeLevelsUpTo = n + 30):
    tau_R (Rb 54 D_3/2) = 164.55 us      (was 83.41 us at 300 K)
    tau_r (Cs 62 D_5/2) = 138.87 us      (was 83.90 us at 300 K)
    tau_P (Rb  7 P_3/2) =   0.270 us     (ARC; BBR negligible at n = 7)

OR-gate config (this report's robust champion):
    Omega_p = 2 pi * 65 MHz,  ratio = 3.0  (-> Omega_R = 2 pi * 195 MHz),
    Omega_c = 2 pi * 60 MHz,  Delta = 2 pi * 500 MHz,  K = 1.0,  alpha = 4.0

CCX-gate config (kept at the rev. 4 amplitudes; sub-pulse times fixed by pi-pulse):
    Omega_cc = 2 pi * 50 MHz,  Omega_t = 2 pi * 20 MHz
    T_cc = pi/Omega_cc = 10 ns,  T_t = pi/Omega_t = 25 ns,  T_tot = 95 ns

Output:  a4_K1_full_fidelity_analysis.json (machine-readable)
         a4_K1_full_fidelity_analysis.log  (human-readable; this run's stdout)
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
import qutip
from scipy.optimize import minimize

from triqg.atoms import CsAtom, RbAtom, composite_basis_state
from triqg.pulses import omega_gaussian
from triqg.hamiltonian import build_hamiltonian, build_ccx_hamiltonian
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
# Physical setup at a = 4.0 um with the redesigned levels
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
R_DD = V_ct_MHz / V_DD_MHz
R_AA = V_ct_MHz / abs(V_cc_MHz)

# =====================================================================
# ARC-verified lifetimes (T = 0 K, spontaneous emission only,
#                         re-checked by check_lifetimes_arc.py)
# =====================================================================
tau_R_us = 164.55   # Rb 54 D_3/2  (ARC, 0 K)
tau_r_us = 138.87   # Cs 62 D_5/2  (ARC, 0 K)
tau_P_us = 0.270    # Rb  7 P_3/2  (ARC, 0 K; BBR negligible at n = 7)
gamma_r = 1.0 / tau_r_us
gamma_R = 1.0 / tau_R_us
gamma_P = 1.0 / tau_P_us
c_ops = build_collapse_operators(gamma_r, gamma_R, gamma_P)

# =====================================================================
# Computational subspace: 8 ket states |c1 c2 t>
# =====================================================================
cs = CsAtom()
rb = RbAtom()
ctrl = [cs.level_index["0"], cs.level_index["1"]]
tgt  = [rb.level_index["A"], rb.level_index["B"]]

comp_kets = []
labels = []
for c1 in range(2):
    for c2 in range(2):
        for t in range(2):
            comp_kets.append(composite_basis_state(ctrl[c1], ctrl[c2], tgt[t]))
            labels.append(f"|{c1},{c2},{['A','B'][t]}>")

# 36 x 8 isometry
V_iso = np.array([k.full().flatten() for k in comp_kets]).T


# Index packing: idx = c1*4 + c2*2 + t   (matches comp_kets order)
def or_perm(idx: int) -> int:
    """OR truth table:  target flips iff (c1 OR c2)."""
    c1 = (idx >> 2) & 1
    c2 = (idx >> 1) & 1
    t  = idx & 1
    t_new = t ^ (c1 | c2)
    return (c1 << 2) | (c2 << 1) | t_new


def ccx_perm(idx: int) -> int:
    """CCX truth table:  target flips iff (c1 AND c2)."""
    c1 = (idx >> 2) & 1
    c2 = (idx >> 1) & 1
    t  = idx & 1
    t_new = t ^ (c1 & c2)
    return (c1 << 2) | (c2 << 1) | t_new


# =====================================================================
# Phase-correction helpers
# =====================================================================
def _F_pro_with_phase(choi_tensor: np.ndarray, U0: np.ndarray,
                      phis: np.ndarray) -> float:
    """Pedersen process fidelity between channel E (via Choi) and (D U0)
    where D = diag(exp(i phis))."""
    d = U0.shape[0]
    D = np.diag(np.exp(1j * phis))
    U_target = D @ U0
    F = np.einsum("ki,ijkl,lj->", U_target.conj(), choi_tensor, U_target) / (d * d)
    return float(np.real(F))


def _optimize_phase(choi_tensor: np.ndarray, U0: np.ndarray,
                    n_restarts: int = 20, seed: int = 0):
    """Find the diagonal phase D = diag(exp(i phi)) maximizing Pedersen
    F_bar. Returns (F_pro_pc, F_bar_pc, phis_best_over_pi)."""
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


# =====================================================================
# Per-gate analyser
# =====================================================================
def analyze_gate(name: str, H, args, t_total, gate_perm,
                 solver_opts, basis_perm_label: str):
    print()
    print("=" * 78)
    print(f"  {name}")
    print("=" * 78)
    print(f"Total gate time:  T_tot = {t_total*1e3:.3f} ns")

    U0 = build_permutation_unitary(comp_kets, gate_perm)

    # ----- (a) Per-input basis fidelities via 8 mesolve runs --------------
    print()
    print("  [a] 8 per-input basis fidelities  (mesolve, with decoherence)")
    F_k = np.zeros(8)
    populations = []   # for diagnostic, not exported
    for k in range(8):
        psi_in = comp_kets[k]
        psi_id = comp_kets[gate_perm(k)]
        tlist  = np.linspace(0, t_total, 350)
        res = simulate(
            method="mesolve", H=H, psi0=psi_in, tlist=tlist,
            c_ops=c_ops, e_ops=[], options=solver_opts, args=args,
        )
        F_k[k] = state_fidelity(res.final_state, psi_id)
        print(f"      F( {labels[k]:<10s} -> {labels[gate_perm(k)]:<10s} ) "
              f"= {F_k[k]:.6f}")

    F_bar_basis = float(F_k.mean())
    print(f"      F_bar_basis  (1/d) sum_k F_k         = {F_bar_basis:.6f}")

    # ----- (b) Pedersen Choi tensor via 64 mesolve runs -------------------
    print()
    print("  [b] Pedersen Choi tensor  (64 mesolve runs, with decoherence)")
    t0 = time.time()
    choi = choi_on_subspace_via_mesolve(
        H, c_ops, t_total, comp_kets,
        args=args, options=solver_opts, verbose=False,
    )
    print(f"      Choi tensor built in {time.time()-t0:.1f} s")

    F_bar_ped, F_pro, survival = pedersen_average_fidelity_from_choi(choi, U0)
    print(f"      Survival T_P                          = {survival:.6f}")
    print(f"      F_pro  (process fidelity, raw)        = {F_pro:.6f}")
    print(f"      F_bar  (Pedersen Haar-averaged, raw)  = {F_bar_ped:.6f}")

    # ----- (c) Noiseless unitary Pedersen ---------------------------------
    print()
    print("  [c] Noiseless unitary Pedersen  (decoherence OFF)")
    t0 = time.time()
    U_full = qutip.propagator(
        H, t_total, c_ops=[], args=args,
        options={"atol": 1e-12, "rtol": 1e-10, "nsteps": 200000},
    )
    U_eff = V_iso.conj().T @ U_full.full() @ V_iso        # 8 x 8 projected
    tr = np.trace(U0.conj().T @ U_eff)
    F_pro_u = abs(tr) ** 2 / 64.0
    F_bar_u = (8 * F_pro_u + 1) / 9.0
    print(f"      noiseless propagator in {time.time()-t0:.1f} s")
    print(f"      F_pro_u   = |Tr(U0^dag U_eff)|^2 / 64 = {F_pro_u:.6f}")
    print(f"      F_bar_u   (Pedersen, no decoherence)  = {F_bar_u:.6f}")

    # Branch phases / magnitudes (D = U_eff U_0^dag, take its diagonal)
    D = U_eff @ U0.conj().T
    diag_D = np.diag(D)
    branch_magnitudes = np.abs(diag_D).tolist()
    branch_phases_pi = (np.angle(diag_D) / np.pi).tolist()

    print("      branch | mag.   | phi/pi   | label")
    for lbl, m, p in zip(labels, branch_magnitudes, branch_phases_pi):
        print(f"             | {m:.4f} | {p:+.4f}  | {lbl}")

    # ----- (d) Phase-corrected Pedersen -----------------------------------
    print()
    print("  [d] Phase-corrected Pedersen  (optimal diagonal D absorbed)")
    t0 = time.time()
    F_pro_pc, F_bar_pc, best_phis = _optimize_phase(choi, U0)
    print(f"      optimization in {time.time()-t0:.2f} s")
    print(f"      F_pro_pc   (process, virtual-Z-absorbed)  = {F_pro_pc:.6f}")
    print(f"      F_bar_pc   (Pedersen + phase correction)  = {F_bar_pc:.6f}")
    print(f"      optimal phases / pi:")
    for lbl, p in zip(labels, best_phis):
        print(f"          {lbl:<10s}  phi/pi = {p/np.pi:+.4f}")

    # ----- collected outputs ---------------------------------------------
    return {
        "name": name,
        "permutation": basis_perm_label,
        "t_total_ns": 1e3 * t_total,
        "labels": labels,
        "F_k_basis_fidelities": F_k.tolist(),
        "F_bar_basis": F_bar_basis,
        "F_pro_Pedersen_raw": F_pro,
        "F_bar_Pedersen_raw": F_bar_ped,
        "survival": survival,
        "F_pro_unitary": F_pro_u,
        "F_bar_unitary": F_bar_u,
        "branch_magnitudes": branch_magnitudes,
        "branch_phases_pi": branch_phases_pi,
        "F_pro_Pedersen_phase_corrected": F_pro_pc,
        "F_bar_Pedersen_phase_corrected": F_bar_pc,
        "optimal_phase_correction_over_pi": (best_phis / np.pi).tolist(),
    }


# =====================================================================
# Header
# =====================================================================
print("=" * 78)
print("  Full fidelity analysis at a = 4.0 um, K = 1")
print("  Energy levels: Rb 54 D_3/2 + Cs 62 D_5/2 (smaller-V_DD redesign)")
print("=" * 78)
print(f"V_ct / (2 pi) = +{V_ct_MHz:7.3f} MHz")
print(f"V_DD / (2 pi) = {V_DD_MHz:+8.4f} MHz  ->  R_DD = {R_DD:.1f}")
print(f"V_cc / (2 pi) = {V_cc_MHz:+8.4f} MHz  ->  R_AA = {R_AA:.1f}")
print(f"Lifetimes (ARC):  tau_R = {tau_R_us:.2f} us,  "
      f"tau_r = {tau_r_us:.2f} us,  tau_P = {tau_P_us:.3f} us")


# =====================================================================
# OR-gate parameters: this report's K = 1 champion
# =====================================================================
omega_p_MHz = 65.0
omega_R_MHz = 3.0 * omega_p_MHz     # = 195.0
omega_c_MHz = 60.0
delta_MHz   = 500.0
K_or        = 1.00
ALPHA       = 4.0
I_inf       = 2 * 0.92770 * 2 ** (-1.0 / 6.0)

omega_p_amp = 2 * np.pi * omega_p_MHz
omega_R_amp = 2 * np.pi * omega_R_MHz
omega_c_amp = 2 * np.pi * omega_c_MHz
delta_or    = 2 * np.pi * delta_MHz
T_c = np.pi / omega_c_amp
sigma_pi4 = ((2 * np.pi * delta_or) / (omega_p_amp ** 2 * I_inf)) ** 3
sigma_or  = (K_or ** 3) * sigma_pi4
T_f       = (ALPHA * sigma_or) ** (1.0 / 3.0)
t_total_or = 2 * T_c + 2 * T_f

H_or = build_hamiltonian(delta_or, V_ct, pulse_p=omega_gaussian, V_cc=V_cc)
args_or = {
    "omega_c_amp": omega_c_amp,
    "omega_p_amp": omega_p_amp,
    "omega_R_amp": omega_R_amp,
    "T_c": T_c, "T_f": T_f, "sigma": sigma_or,
}
solver_opts_or = {"store_final_state": True, "nsteps": 200000,
                  "atol": 1e-10, "rtol": 1e-8}

print()
print(f"OR drive params (K = 1 champion):")
print(f"  Omega_p / (2 pi) = {omega_p_MHz:.1f} MHz")
print(f"  Omega_R / Omega_p = {omega_R_MHz/omega_p_MHz:.2f}")
print(f"  Omega_c / (2 pi) = {omega_c_MHz:.1f} MHz")
print(f"  Delta   / (2 pi) = {delta_MHz:.1f} MHz")
print(f"  K = {K_or},  alpha = {ALPHA}")
print(f"  sigma = {sigma_or*1e3:.3f} ns, T_f = {T_f*1e3:.2f} ns, "
      f"T_c = {T_c*1e3:.2f} ns")
print(f"  T_tot_OR = {t_total_or*1e3:.2f} ns")

results_OR = analyze_gate(
    name="OR gate  (K = 1, a = 4 um champion)",
    H=H_or, args=args_or, t_total=t_total_or,
    gate_perm=or_perm, solver_opts=solver_opts_or,
    basis_perm_label="t XOR (c1 OR c2)",
)


# =====================================================================
# CCX-gate parameters: rev. 4 amplitudes, same a = 4 um lattice
# =====================================================================
omega_cc_MHz = 50.0
omega_t_MHz  = 20.0
omega_cc_amp = 2 * np.pi * omega_cc_MHz
omega_t_amp  = 2 * np.pi * omega_t_MHz
T_cc = np.pi / omega_cc_amp
T_t  = np.pi / omega_t_amp
t_total_ccx = 2 * T_cc + 3 * T_t

H_ccx = build_ccx_hamiltonian(V_ct, V_cc=V_cc)
args_ccx = {
    "omega_cc_amp": omega_cc_amp,
    "omega_t_amp":  omega_t_amp,
    "T_cc": T_cc, "T_t": T_t,
}
max_freq = max(omega_cc_amp, omega_t_amp) / (2 * np.pi)
solver_opts_ccx = {"store_final_state": True, "nsteps": 200000,
                   "max_step": 1.0 / (20 * max_freq),
                   "atol": 1e-10, "rtol": 1e-8}

print()
print(f"CCX drive params (rev. 4 amplitudes, same lattice):")
print(f"  Omega_cc / (2 pi) = {omega_cc_MHz:.1f} MHz  -> T_cc = {T_cc*1e3:.2f} ns")
print(f"  Omega_t  / (2 pi) = {omega_t_MHz:.1f} MHz  -> T_t  = {T_t *1e3:.2f} ns")
print(f"  T_tot_CCX = 2 T_cc + 3 T_t = {t_total_ccx*1e3:.2f} ns")

results_CCX = analyze_gate(
    name="CCX gate  (rev. 4 amplitudes, a = 4 um lattice)",
    H=H_ccx, args=args_ccx, t_total=t_total_ccx,
    gate_perm=ccx_perm, solver_opts=solver_opts_ccx,
    basis_perm_label="t XOR (c1 AND c2)",
)


# =====================================================================
# Persist
# =====================================================================
out_path = Path(__file__).parent / "a4_K1_full_fidelity_analysis.json"
all_results = {
    "physical_params": {
        "a_um": a_um,
        "r_DA_um": r_DA, "r_DD_um": r_DD, "r_AA_um": r_AA,
        "C3_tilde_GHz_um3": C3_tilde,
        "C6_RbRb_GHz_um6":  C6_RbRb,
        "C6_CsCs_GHz_um6":  C6_CsCs,
        "V_ct_over_2pi_MHz": V_ct_MHz,
        "V_DD_over_2pi_MHz": V_DD_MHz,
        "V_cc_over_2pi_MHz": V_cc_MHz,
        "R_DD": R_DD, "R_AA": R_AA,
        "tau_r_us_ARC": tau_r_us,
        "tau_R_us_ARC": tau_R_us,
        "tau_P_us":     tau_P_us,
    },
    "OR_gate": results_OR,
    "CCX_gate": results_CCX,
}
out_path.write_text(json.dumps(all_results, indent=2))
print()
print("=" * 78)
print(f"Wrote: {out_path}")
print("=" * 78)


# =====================================================================
# Final summary table
# =====================================================================
def fmt(x): return f"{x:.6f}" if isinstance(x, float) else str(x)

print()
print("Summary  (a = 4.0 um, K = 1, ARC lifetimes; all values dimensionless):")
print("=" * 78)
print(f"  {'Metric':<42s} {'OR gate':>14s} {'CCX gate':>14s}")
print("-" * 78)
print(f"  {'F_bar_basis (1/d) sum F_k':<42s} "
      f"{fmt(results_OR['F_bar_basis']):>14s} "
      f"{fmt(results_CCX['F_bar_basis']):>14s}")
print(f"  {'F_bar_Pedersen (raw)':<42s} "
      f"{fmt(results_OR['F_bar_Pedersen_raw']):>14s} "
      f"{fmt(results_CCX['F_bar_Pedersen_raw']):>14s}")
print(f"  {'F_bar_Pedersen_phase_corrected':<42s} "
      f"{fmt(results_OR['F_bar_Pedersen_phase_corrected']):>14s} "
      f"{fmt(results_CCX['F_bar_Pedersen_phase_corrected']):>14s}")
print(f"  {'F_bar_unitary (no decoherence)':<42s} "
      f"{fmt(results_OR['F_bar_unitary']):>14s} "
      f"{fmt(results_CCX['F_bar_unitary']):>14s}")
print(f"  {'Computational survival T_P':<42s} "
      f"{fmt(results_OR['survival']):>14s} "
      f"{fmt(results_CCX['survival']):>14s}")
print(f"  {'T_tot (ns)':<42s} "
      f"{results_OR['t_total_ns']:>14.3f} "
      f"{results_CCX['t_total_ns']:>14.3f}")
print("=" * 78)

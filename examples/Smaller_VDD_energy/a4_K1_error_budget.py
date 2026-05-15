"""
Error-budget decomposition of the a = 4.0 um, K = 1 OR-gate champion
(and the rev. 4 CCX gate at the same lattice).

Three sources of error are claimed in the protocol:
    (P)  intermediate Rb 7 P_3/2 state spontaneous decay   [gamma_P]
    (Y)  Rydberg-state spontaneous decay
            |R> = Rb 54 D_3/2    target                    [gamma_R]
            |r> = Cs 62 D_5/2    ancillas                  [gamma_r]
    (V)  same-species van der Waals interaction
            V_cc = -2.78 MHz  on |r r, *>    Cs-Cs (ancilla-ancilla)
            V_DD = +2.21 MHz  on |R R> :  no Rb-Rb pair exists in the 3-atom
                  box (only one target), so it does NOT enter this simulation.
                  We re-check this by *not* changing anything but V_cc.

For each gate we run *five* simulations on the champion drive parameters:

    "baseline"          full physics (gamma_P, gamma_R, gamma_r, V_cc) on
    "no_P_decay"        gamma_P  := 0,   everything else unchanged
    "no_Rydberg_decay"  gamma_R  := 0, gamma_r := 0
    "no_same_species"   V_cc     := 0   (and V_DD, but V_DD has no effect)
    "clean"             everything off  -> pure coherent floor

For each scenario we report
    F_bar_basis    (1/d) sum_k F_k       (classical truth-table average)
    1 - F_bar      basis-infidelity
    delta_F        = F_bar - F_bar_baseline   (how much the error costs)
    share          = delta_F / (F_bar_clean - F_bar_baseline)  in %

The "share" column is the fraction of total infidelity removed by that
single channel.  Shares sum to ~100% (modulo small cross-terms between
channels, which we print explicitly).

Run time: ~1.5 min on Apple Silicon (10 x 8 mesolve calls).

Outputs:
    a4_K1_error_budget.json   (machine-readable)
    a4_K1_error_budget.log    (this run's stdout)
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np

from triqg.atoms import CsAtom, RbAtom, composite_basis_state
from triqg.pulses import omega_gaussian
from triqg.hamiltonian import build_hamiltonian, build_ccx_hamiltonian
from triqg.decoherence import build_collapse_operators
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
V_DD_MHz = 1000.0 * C6_RbRb  / r_DD ** 6   # = 2.214 MHz   (would-be Rb-Rb)
V_cc_MHz = 1000.0 * C6_CsCs  / r_AA ** 6   # = -2.778 MHz  Cs-Cs
V_ct = 2 * np.pi * V_ct_MHz
V_cc_full = 2 * np.pi * V_cc_MHz

# =====================================================================
# ARC-verified lifetimes (T = 0 K, spontaneous emission only)
# =====================================================================
tau_R_us = 164.55
tau_r_us = 138.87
tau_P_us = 0.270
gamma_r_full = 1.0 / tau_r_us
gamma_R_full = 1.0 / tau_R_us
gamma_P_full = 1.0 / tau_P_us

# =====================================================================
# Computational subspace (8 kets) and permutations
# =====================================================================
cs = CsAtom(); rb = RbAtom()
ctrl = [cs.level_index["0"], cs.level_index["1"]]
tgt  = [rb.level_index["A"], rb.level_index["B"]]

comp_kets = []
labels = []
for c1 in range(2):
    for c2 in range(2):
        for t in range(2):
            comp_kets.append(composite_basis_state(ctrl[c1], ctrl[c2], tgt[t]))
            labels.append(f"|{c1},{c2},{['A','B'][t]}>")


def or_perm(idx: int) -> int:
    c1 = (idx >> 2) & 1; c2 = (idx >> 1) & 1; t = idx & 1
    return (c1 << 2) | (c2 << 1) | (t ^ (c1 | c2))


def ccx_perm(idx: int) -> int:
    c1 = (idx >> 2) & 1; c2 = (idx >> 1) & 1; t = idx & 1
    return (c1 << 2) | (c2 << 1) | (t ^ (c1 & c2))


# =====================================================================
# OR-gate champion drive
# =====================================================================
omega_p_MHz = 65.0
omega_R_MHz = 3.0 * omega_p_MHz
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

args_or = {
    "omega_c_amp": omega_c_amp, "omega_p_amp": omega_p_amp,
    "omega_R_amp": omega_R_amp,
    "T_c": T_c, "T_f": T_f, "sigma": sigma_or,
}
solver_opts_or = {"store_final_state": True, "nsteps": 200000,
                  "atol": 1e-10, "rtol": 1e-8}

# =====================================================================
# CCX-gate drive (rev. 4 amplitudes)
# =====================================================================
omega_cc_MHz = 50.0
omega_t_MHz  = 20.0
omega_cc_amp = 2 * np.pi * omega_cc_MHz
omega_t_amp  = 2 * np.pi * omega_t_MHz
T_cc = np.pi / omega_cc_amp
T_t  = np.pi / omega_t_amp
t_total_ccx = 2 * T_cc + 3 * T_t

args_ccx = {
    "omega_cc_amp": omega_cc_amp, "omega_t_amp": omega_t_amp,
    "T_cc": T_cc, "T_t": T_t,
}
max_freq = max(omega_cc_amp, omega_t_amp) / (2 * np.pi)
solver_opts_ccx = {"store_final_state": True, "nsteps": 200000,
                   "max_step": 1.0 / (20 * max_freq),
                   "atol": 1e-10, "rtol": 1e-8}


# =====================================================================
# Run one scenario
# =====================================================================
def run_scenario(name: str, build_H, gate_perm, args, solver_opts,
                 t_total: float,
                 gamma_r: float, gamma_R: float, gamma_P: float,
                 V_cc: float) -> dict:
    """Run all 8 basis inputs and return the average basis fidelity."""
    H = build_H(V_cc)
    c_ops = build_collapse_operators(gamma_r, gamma_R, gamma_P)

    F_k = np.zeros(8)
    tlist = np.linspace(0, t_total, 350)
    t0 = time.time()
    for k in range(8):
        psi_in = comp_kets[k]
        psi_id = comp_kets[gate_perm(k)]
        res = simulate(
            method="mesolve", H=H, psi0=psi_in, tlist=tlist,
            c_ops=c_ops, e_ops=[], options=solver_opts, args=args,
        )
        F_k[k] = state_fidelity(res.final_state, psi_id)
    elapsed = time.time() - t0
    F_bar = float(F_k.mean())

    return {
        "name": name,
        "gamma_r": gamma_r, "gamma_R": gamma_R, "gamma_P": gamma_P,
        "V_cc_over_2pi_MHz": V_cc / (2 * np.pi),
        "F_k": F_k.tolist(),
        "F_bar_basis": F_bar,
        "infidelity": 1.0 - F_bar,
        "elapsed_s": elapsed,
    }


# =====================================================================
# Header
# =====================================================================
print("=" * 78)
print("  Error budget at a = 4.0 um, K = 1  (OR + CCX gates)")
print("=" * 78)
print(f"V_ct / (2 pi) = +{V_ct_MHz:7.3f} MHz")
print(f"V_DD / (2 pi) = {V_DD_MHz:+8.4f} MHz   (would-be Rb-Rb; not in 3-atom box)")
print(f"V_cc / (2 pi) = {V_cc_MHz:+8.4f} MHz   (Cs-Cs ancilla-ancilla, enters H)")
print(f"gamma_P = {gamma_P_full:.3f} / us   (tau_P = {tau_P_us:.3f} us, Rb 7P_3/2)")
print(f"gamma_R = {gamma_R_full:.5f} / us   (tau_R = {tau_R_us:.2f} us, Rb 54D_3/2)")
print(f"gamma_r = {gamma_r_full:.5f} / us   (tau_r = {tau_r_us:.2f} us, Cs 62D_5/2)")


# =====================================================================
# Scenario factory
# =====================================================================
SCENARIOS = [
    # (label, gamma_r, gamma_R, gamma_P, V_cc)
    ("baseline",           gamma_r_full, gamma_R_full, gamma_P_full, V_cc_full),
    ("no_P_decay",         gamma_r_full, gamma_R_full, 0.0,           V_cc_full),
    ("no_Rydberg_decay",   0.0,           0.0,           gamma_P_full, V_cc_full),
    ("no_same_species",    gamma_r_full, gamma_R_full, gamma_P_full, 0.0),
    ("clean",              0.0,           0.0,           0.0,           0.0),
]


def run_gate(name: str, build_H, gate_perm, args, solver_opts, t_total):
    print()
    print("=" * 78)
    print(f"  {name}")
    print("=" * 78)
    print(f"T_tot = {t_total*1e3:.2f} ns")
    print()

    results = []
    for (sc_name, gr, gR, gP, Vcc) in SCENARIOS:
        print(f"  running scenario  {sc_name:<22s}  "
              f"(gP={gP:.3f}, gR={gR:.4f}, gr={gr:.4f}, "
              f"Vcc/2pi={Vcc/(2*np.pi):+7.3f} MHz)")
        r = run_scenario(sc_name, build_H, gate_perm, args, solver_opts,
                         t_total, gr, gR, gP, Vcc)
        print(f"    F_bar = {r['F_bar_basis']:.6f}, "
              f"1-F = {r['infidelity']:.3e}, "
              f"({r['elapsed_s']:.1f} s)")
        results.append(r)
    return results


# =====================================================================
# Build-H closures
# =====================================================================
def build_or_H(V_cc_val):
    return build_hamiltonian(delta_or, V_ct, pulse_p=omega_gaussian, V_cc=V_cc_val)


def build_ccx_H(V_cc_val):
    return build_ccx_hamiltonian(V_ct, V_cc=V_cc_val)


# =====================================================================
# Run both gates
# =====================================================================
OR_results = run_gate(
    "OR gate (K = 1, a = 4 um champion)",
    build_or_H, or_perm, args_or, solver_opts_or, t_total_or,
)

CCX_results = run_gate(
    "CCX gate (rev. 4 amplitudes, a = 4 um lattice)",
    build_ccx_H, ccx_perm, args_ccx, solver_opts_ccx, t_total_ccx,
)


# =====================================================================
# Decomposition table
# =====================================================================
def print_budget(gate_name: str, results: list):
    by_name = {r["name"]: r for r in results}
    F_base  = by_name["baseline"]["F_bar_basis"]
    F_clean = by_name["clean"]["F_bar_basis"]
    inf_base = 1.0 - F_base
    inf_clean = 1.0 - F_clean
    total_decoh = inf_base - inf_clean

    print()
    print("=" * 78)
    print(f"  ERROR BUDGET  ({gate_name})")
    print("=" * 78)
    print(f"  Baseline infidelity                       1 - F = {inf_base:.4e}")
    print(f"  Coherent floor (clean run)                1 - F = {inf_clean:.4e}")
    print(f"  Total decoherence + V_cc contribution     {total_decoh:.4e}")
    print()
    print(f"  {'Scenario':<22s} {'F_bar':>9s} {'1-F':>11s} "
          f"{'delta_F':>11s} {'share':>9s}")
    print(f"  {'-'*22} {'-'*9} {'-'*11} {'-'*11} {'-'*9}")
    for r in results:
        F = r["F_bar_basis"]
        inf = 1.0 - F
        delta_F = F - F_base               # >0 = error removed
        if r["name"] == "baseline":
            share_str = "-"
        elif r["name"] == "clean":
            share_str = "100.00%"
        else:
            share = 100.0 * delta_F / total_decoh if total_decoh > 0 else 0.0
            share_str = f"{share:6.2f}%"
        print(f"  {r['name']:<22s} {F:>9.6f} {inf:>11.3e} "
              f"{delta_F:>+11.3e} {share_str:>9s}")

    # Cross-term sanity check
    F_noP   = by_name["no_P_decay"]["F_bar_basis"] - F_base
    F_noYR  = by_name["no_Rydberg_decay"]["F_bar_basis"] - F_base
    F_noVcc = by_name["no_same_species"]["F_bar_basis"] - F_base
    F_clean_gain = F_clean - F_base
    sum_singles = F_noP + F_noYR + F_noVcc
    cross = F_clean_gain - sum_singles
    print()
    print(f"  Linearity check:")
    print(f"    sum of single-channel gains  = {sum_singles:.4e}")
    print(f"    clean - baseline             = {F_clean_gain:.4e}")
    print(f"    cross-term (clean-sum)       = {cross:+.4e}   "
          f"({100*cross/F_clean_gain:+.2f}% of total)")
    print()
    # Headline verdict
    rank = sorted(
        [("P_decay",        F_noP),
         ("Rydberg_decay",  F_noYR),
         ("same_species",   F_noVcc)],
        key=lambda kv: -kv[1],
    )
    print(f"  Dominant error channel  -->  {rank[0][0]}  "
          f"(removes delta_F = {rank[0][1]:.3e})")
    print(f"  Ranking (largest -> smallest):")
    for tag, dF in rank:
        share = 100.0 * dF / total_decoh if total_decoh > 0 else 0.0
        print(f"    {tag:<14s}  delta_F = {dF:+.3e}   share = {share:5.2f}%")
    return {"baseline_infidelity": inf_base,
            "coherent_floor": inf_clean,
            "total_decoh": total_decoh,
            "delta_F_no_P":      F_noP,
            "delta_F_no_Ryd":    F_noYR,
            "delta_F_no_Vcc":    F_noVcc,
            "cross_term":        cross,
            "ranking": rank}


OR_budget  = print_budget("OR gate",  OR_results)
CCX_budget = print_budget("CCX gate", CCX_results)


# =====================================================================
# Persist
# =====================================================================
out = {
    "physical_params": {
        "a_um": a_um, "r_DA_um": r_DA, "r_DD_um": r_DD, "r_AA_um": r_AA,
        "V_ct_over_2pi_MHz": V_ct_MHz,
        "V_DD_over_2pi_MHz": V_DD_MHz,
        "V_cc_over_2pi_MHz": V_cc_MHz,
        "tau_R_us": tau_R_us, "tau_r_us": tau_r_us, "tau_P_us": tau_P_us,
        "gamma_R": gamma_R_full, "gamma_r": gamma_r_full,
        "gamma_P": gamma_P_full,
    },
    "OR_gate": {
        "drive": {
            "omega_p_MHz": omega_p_MHz, "omega_R_MHz": omega_R_MHz,
            "omega_c_MHz": omega_c_MHz, "delta_MHz": delta_MHz,
            "K": K_or, "alpha": ALPHA, "T_tot_ns": 1e3 * t_total_or,
        },
        "scenarios": OR_results,
        "budget": OR_budget,
    },
    "CCX_gate": {
        "drive": {
            "omega_cc_MHz": omega_cc_MHz, "omega_t_MHz": omega_t_MHz,
            "T_tot_ns": 1e3 * t_total_ccx,
        },
        "scenarios": CCX_results,
        "budget": CCX_budget,
    },
}

out_path = Path(__file__).parent / "a4_K1_error_budget.json"
out_path.write_text(json.dumps(out, indent=2, default=float))
print()
print("=" * 78)
print(f"Wrote: {out_path}")
print("=" * 78)

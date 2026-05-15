"""
Candidate "Rev. 4" OR-gate configuration at a = 4 um.

Goal:  beat the FINAL F_bar = 0.99282 by suppressing the V_cc-during-pi-pulse
       error that bit Case B (a = 4, Delta = 500, K = 0.98) -- that
       configuration matched FINAL on M_2 and gate time but lost 0.6% on
       |11,*> from imperfect Cs pi-rotation when both controls move
       |1,1> -> |r,r> simultaneously while V_cc detunes the |r,r>
       component by ~5.6% of Omega_c.

Two changes from Case B / Case C:
  - Delta retuned from 500 -> 250 MHz so M_2 ~ 25 (same as FINAL).
  - Omega_c raised from 50 -> {75, 100, 150, 200} MHz so T_c = pi/Omega_c
    shrinks by the same factor. The V_cc-induced pi-pulse error
    scales as (V_cc / Omega_c)^2.

Predicted scaling at a = 4 um (V_cc = -2.78 MHz):
    Omega_c = 50  MHz:  (V_cc/Omega_c)^2 = (2.78/50)^2  = 3.1e-3
                        x 4 (2 pi-pulses x 2 controls) = 1.24%
    Omega_c = 100 MHz:  (V_cc/Omega_c)^2 = (2.78/100)^2 = 7.7e-4
                        x 4                              = 0.31%
    Omega_c = 200 MHz:  (V_cc/Omega_c)^2 = (2.78/200)^2 = 1.9e-4
                        x 4                              = 0.08%

Combined decay + leakage + Vcc-pi error projection:
    Omega_c=200, Delta=250, K=0.97 at a=4 um, tau_r=83.9 us:
      decay floor          : 2*206/83900 = 0.49%
      blockade leakage     : ~1/M_2^2     ~ 0.16% (M_2=25)
      Vcc-pi-pulse error   : ~0.08%
      sub-total            : ~0.73%
      ->  F_bar ~ 0.9927    (small improvement over FINAL 0.99282)

The win is gate speed: 200 ns instead of 385 ns (3x faster cycle budget),
at comparable F_bar.

Lifetimes: ARC-verified (tau_r = 83.90 us, tau_R = 83.41 us).

Cases:
    REF   : FINAL config (a=5, Omega_c=50, Delta=500, K=0.95)  -- baseline
    Sweep : a=4, Delta=250, K=0.97, Omega_c in {50, 75, 100, 150, 200} MHz
    Final : best Omega_c with K refinement
"""
import math
import numpy as np

from triqg.atoms import CsAtom, RbAtom, composite_basis_state
from triqg.pulses import omega_gaussian
from triqg.hamiltonian import build_hamiltonian
from triqg.decoherence import build_collapse_operators
from triqg.solver import simulate
from triqg.analysis import state_fidelity, average_gate_fidelity
from level_parameters import C3_tilde, C6_RbRb, C6_CsCs

# ARC lifetimes
tau_r_us = 83.90
tau_R_us = 83.41
tau_P_us = 0.131
gamma_r = 1.0 / tau_r_us
gamma_R = 1.0 / tau_R_us
gamma_P = 1.0 / tau_P_us
c_ops = build_collapse_operators(gamma_r, gamma_R, gamma_P)

omega_p_MHz = 50.0
omega_R_MHz = 3.5 * omega_p_MHz
omega_p_amp = 2 * np.pi * omega_p_MHz
omega_R_amp = 2 * np.pi * omega_R_MHz
ALPHA = 4.0
I_inf = 2 * 0.92770 * 2**(-1.0 / 6.0)

cs = CsAtom(); rb = RbAtom()
ctrl_levels = [cs.level_index["0"], cs.level_index["1"]]
tgt_levels  = [rb.level_index["A"], rb.level_index["B"]]
basis_inputs = []
for c1 in range(2):
    for c2 in range(2):
        for t in range(2):
            psi_in = composite_basis_state(ctrl_levels[c1], ctrl_levels[c2], tgt_levels[t])
            psi_id = composite_basis_state(
                ctrl_levels[c1], ctrl_levels[c2],
                tgt_levels[1 - t if (c1 or c2) else t]
            )
            label = f"|{c1}{c2}{['A','B'][t]}>"
            basis_inputs.append((label, psi_in, psi_id))


def interactions_at(a_um):
    r_DA = a_um / math.sqrt(2)
    r_AA = a_um * math.sqrt(2)
    V_ct_MHz = 1000 * C3_tilde / r_DA**3
    V_cc_MHz = 1000 * C6_CsCs  / r_AA**6
    return V_ct_MHz, V_cc_MHz


def run_one(a_um, delta_MHz, K, omega_c_MHz, ntime=400):
    V_ct_MHz, V_cc_MHz = interactions_at(a_um)
    V_ct = 2 * np.pi * V_ct_MHz
    V_cc = 2 * np.pi * V_cc_MHz
    delta = 2 * np.pi * delta_MHz
    omega_c_amp = 2 * np.pi * omega_c_MHz
    T_c = np.pi / omega_c_amp
    sigma_pi4 = ((2 * np.pi * delta) / (omega_p_amp**2 * I_inf)) ** 3
    sigma = (K ** 3) * sigma_pi4
    T_f = (ALPHA * sigma) ** (1.0 / 3.0)
    t_total = 2 * T_c + 2 * T_f
    args = {"omega_c_amp": omega_c_amp, "omega_p_amp": omega_p_amp,
            "omega_R_amp": omega_R_amp, "T_c": T_c, "T_f": T_f, "sigma": sigma}
    H = build_hamiltonian(delta, V_ct, pulse_p=omega_gaussian, V_cc=V_cc)
    tlist = np.linspace(0, t_total, ntime)
    opts = {"store_final_state": True, "nsteps": 200000}
    pairs = []; per_F = []
    for lbl, psi_in, psi_id in basis_inputs:
        res = simulate(method="mesolve", H=H, psi0=psi_in, tlist=tlist,
                       c_ops=c_ops, e_ops=[], options=opts, args=args)
        pairs.append((res.final_state, psi_id))
        per_F.append(state_fidelity(res.final_state, psi_id))
    F_bar = average_gate_fidelity(pairs)
    M2 = V_ct_MHz / (omega_p_MHz * omega_R_MHz / (2 * delta_MHz))
    return F_bar, per_F, T_c * 1e3, T_f * 1e3, t_total * 1e3, M2


# =====================================================================
# REF: FINAL with ARC lifetimes  (a = 5, Omega_c = 50, Delta = 500, K = 0.95)
# =====================================================================
print("=" * 84)
print("REF:  FINAL config with ARC lifetimes")
print("      a=5 um, Omega_c=50 MHz, Delta=500 MHz, K=0.95")
print("=" * 84)
F_REF, perF_REF, Tc_REF, Tf_REF, gate_REF, M2_REF = run_one(5.00, 500.0, 0.95, 50.0)
print(f"T_c = {Tc_REF:.2f} ns,  T_f = {Tf_REF:.2f} ns,  gate = {gate_REF:.2f} ns,  "
      f"M_2 = {M2_REF:.2f}")
print(f"F_bar = {F_REF:.6f},  1-F = {1-F_REF:.3e}")
for (lbl, _, _), F in zip(basis_inputs, perF_REF):
    print(f"   {lbl:<8s}  F = {F:.6f}")

# =====================================================================
# Sweep Omega_c at (a=4, Delta=250, K=0.97)
# =====================================================================
print()
print("=" * 84)
print("Sweep Omega_c at  a = 4 um,  Delta = 250 MHz,  K = 0.97")
print("=" * 84)
omega_c_list = [50.0, 75.0, 100.0, 150.0, 200.0]
print(f"{'Oc':>5s} {'T_c':>7s} {'gate':>8s} {'M2':>7s} {'F_bar':>10s} {'1-F':>10s} "
      f"{'(Vcc/Oc)^2*4':>12s}")
print("-" * 78)
V_cc_a4 = 2.7771  # MHz
sweep = []
for oc in omega_c_list:
    F, perF, Tcn, Tfn, gn, M2 = run_one(4.00, 250.0, 0.97, oc)
    proj = (V_cc_a4 / oc)**2 * 4
    print(f"{oc:>5.0f} {Tcn:>7.2f} {gn:>8.2f} {M2:>7.2f} "
          f"{F:>10.6f} {1-F:>10.3e} {proj:>12.4e}")
    sweep.append((oc, F, perF, Tcn, Tfn, gn))

best = max(sweep, key=lambda x: x[1])
oc_best, F_best, perF_best, Tc_best, Tf_best, gate_best = best
print(f"\nBest Omega_c = {oc_best:.0f} MHz:  F_bar = {F_best:.6f}, gate = {gate_best:.2f} ns")
print("Per-input at best Omega_c:")
for (lbl, _, _), F in zip(basis_inputs, perF_best):
    print(f"   {lbl:<8s}  F = {F:.6f}")

# =====================================================================
# K refinement at the best Omega_c
# =====================================================================
print()
print("=" * 84)
print(f"K refinement at a = 4 um, Delta = 250 MHz, Omega_c = {oc_best:.0f} MHz")
print("=" * 84)
Ks = [0.94, 0.95, 0.96, 0.97, 0.98, 1.00]
print(f"{'K':>5s} {'gate':>8s} {'F_bar':>10s} {'1-F':>10s}")
print("-" * 45)
best_KF = F_best; best_K = 0.97; best_K_perF = perF_best; best_K_gate = gate_best
for K in Ks:
    F, perF, _, _, gn, _ = run_one(4.00, 250.0, K, oc_best)
    print(f"{K:>5.2f} {gn:>8.2f} {F:>10.6f} {1-F:>10.3e}")
    if F > best_KF:
        best_KF = F; best_K = K; best_K_perF = perF; best_K_gate = gn

print(f"\nBest K = {best_K:.2f}:  F_bar = {best_KF:.6f}, gate = {best_K_gate:.2f} ns")
for (lbl, _, _), F in zip(basis_inputs, best_K_perF):
    print(f"   {lbl:<8s}  F = {F:.6f}")

# =====================================================================
# Summary
# =====================================================================
print()
print("=" * 84)
print("SUMMARY  (ARC lifetimes throughout)")
print("=" * 84)
print(f"{'Config':<40s} {'gate':>8s} {'F_bar':>10s} {'1-F':>10s}")
print("-" * 75)
print(f"{'REF: a=5, Oc=50, D=500, K=0.95':<40s} {gate_REF:>8.2f} "
      f"{F_REF:>10.6f} {1-F_REF:>10.3e}")
print(f"{'Rev4 candidate: a=4, Oc='+str(int(oc_best))+', D=250, K='+f'{best_K:.2f}':<40s} "
      f"{best_K_gate:>8.2f} {best_KF:>10.6f} {1-best_KF:>10.3e}")

print()
print("Decay floor on |11,*>:")
print(f"   REF:           {2*gate_REF/1000/tau_r_us*100:.2f}%   "
      f"simulated 1-F_(11,*) = {(1-perF_REF[6])*100:.2f}%")
print(f"   Rev4 candidate: {2*best_K_gate/1000/tau_r_us*100:.2f}%   "
      f"simulated 1-F_(11,*) = {(1-best_K_perF[6])*100:.2f}%")

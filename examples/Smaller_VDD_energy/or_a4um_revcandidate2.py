"""
Second attempt at Rev4 candidate at a = 4 um.

Lesson from revcandidate.py:  smaller Delta degrades |00,*> EIT shielding
faster than it improves |11,*> decay.  Net F_bar drops.

Strategy now:  keep Delta = 500 MHz (preserving EIT) and sweep Omega_c to
target the residual |11,*> loss at a = 4 um, which is 0.6% above the FINAL
floor and presumably V_cc-driven.

Cases:
    a = 4 um, Delta = 500 MHz, K = 0.98 (Case B optimum from arclife scan)
    Omega_c in {50, 100, 150, 200} MHz

ARC lifetimes throughout.
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

tau_r_us, tau_R_us, tau_P_us = 83.90, 83.41, 0.131
gamma_r, gamma_R, gamma_P = 1/tau_r_us, 1/tau_R_us, 1/tau_P_us
c_ops = build_collapse_operators(gamma_r, gamma_R, gamma_P)

omega_p_amp = 2 * np.pi * 50.0
omega_R_amp = 2 * np.pi * 175.0
ALPHA = 4.0
I_inf = 2 * 0.92770 * 2**(-1/6)

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
            basis_inputs.append((f"|{c1}{c2}{['A','B'][t]}>", psi_in, psi_id))


def interactions_at(a_um):
    return (1000*C3_tilde/(a_um/math.sqrt(2))**3,
            1000*C6_CsCs/(a_um*math.sqrt(2))**6)


def run_one(a_um, delta_MHz, K, omega_c_MHz):
    V_ct_MHz, V_cc_MHz = interactions_at(a_um)
    V_ct = 2*np.pi*V_ct_MHz; V_cc = 2*np.pi*V_cc_MHz
    delta = 2*np.pi*delta_MHz
    omega_c_amp = 2*np.pi*omega_c_MHz
    T_c = np.pi/omega_c_amp
    sigma_pi4 = ((2*np.pi*delta)/(omega_p_amp**2 * I_inf))**3
    sigma = K**3 * sigma_pi4
    T_f = (ALPHA*sigma)**(1/3)
    t_total = 2*T_c + 2*T_f
    args = {"omega_c_amp": omega_c_amp, "omega_p_amp": omega_p_amp,
            "omega_R_amp": omega_R_amp, "T_c": T_c, "T_f": T_f, "sigma": sigma}
    H = build_hamiltonian(delta, V_ct, pulse_p=omega_gaussian, V_cc=V_cc)
    tlist = np.linspace(0, t_total, 400)
    opts = {"store_final_state": True, "nsteps": 200000}
    pairs = []; per_F = []
    for lbl, psi_in, psi_id in basis_inputs:
        res = simulate(method="mesolve", H=H, psi0=psi_in, tlist=tlist,
                       c_ops=c_ops, e_ops=[], options=opts, args=args)
        pairs.append((res.final_state, psi_id))
        per_F.append(state_fidelity(res.final_state, psi_id))
    F_bar = average_gate_fidelity(pairs)
    return F_bar, per_F, T_c*1e3, t_total*1e3


# REF
print("=" * 80)
print("REF: a=5, Omega_c=50, Delta=500, K=0.95 (FINAL with ARC lifetimes)")
F_REF, perF_REF, Tc_REF, gate_REF = run_one(5.00, 500.0, 0.95, 50.0)
print(f"gate = {gate_REF:.2f} ns,  F_bar = {F_REF:.6f}")
print("  per-input:", " ".join([f"{F:.4f}" for F in perF_REF]))

# Omega_c sweep at a=4, Delta=500, K=0.98
print()
print("=" * 80)
print("Sweep Omega_c at a = 4 um, Delta = 500 MHz, K = 0.98 (Case B optimum)")
print("=" * 80)
print(f"{'Oc':>5s} {'T_c':>7s} {'gate':>8s} {'F_bar':>10s} {'1-F':>10s} "
      f"{'F|00*>':>10s} {'F|01*>':>10s} {'F|11*>':>10s}")
print("-" * 90)
sweep = []
for oc in [50.0, 75.0, 100.0, 150.0, 200.0]:
    F, perF, Tcn, gn = run_one(4.00, 500.0, 0.98, oc)
    print(f"{oc:>5.0f} {Tcn:>7.2f} {gn:>8.2f} {F:>10.6f} {1-F:>10.3e} "
          f"{perF[0]:>10.6f} {perF[2]:>10.6f} {perF[6]:>10.6f}")
    sweep.append((oc, F, perF, gn))

# Also try (a=4, Delta=400) and (a=4, Delta=600) as middle-ground EIT/decay tradeoffs
print()
print("=" * 80)
print("Other (Delta, K) at a = 4 um, Omega_c = 50 MHz (cheap probe)")
print("=" * 80)
print(f"{'Delta':>7s} {'K':>6s} {'gate':>8s} {'F_bar':>10s} "
      f"{'F|00*>':>10s} {'F|01*>':>10s} {'F|11*>':>10s}")
print("-" * 80)
for d, k in [(400.0, 0.96), (500.0, 0.98), (600.0, 0.98), (700.0, 0.99), (800.0, 1.00)]:
    F, perF, _, gn = run_one(4.00, d, k, 50.0)
    print(f"{d:>7.0f} {k:>6.2f} {gn:>8.2f} {F:>10.6f} "
          f"{perF[0]:>10.6f} {perF[2]:>10.6f} {perF[6]:>10.6f}")

# Summary
print()
best = max(sweep, key=lambda x: x[1])
print("=" * 80)
print("SUMMARY")
print("=" * 80)
print(f"REF:                                              gate = {gate_REF:.1f} ns,  "
      f"F_bar = {F_REF:.6f}")
print(f"Best Omega_c (a=4, D=500, K=0.98):  Oc = {best[0]:.0f} MHz,  "
      f"gate = {best[3]:.1f} ns,  F_bar = {best[1]:.6f}")

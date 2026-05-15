"""
Re-run the a = 5 vs a = 4 comparison with two corrections:

  1.  ARC-verified lifetimes  (replacing the n^3-scaled estimates)
        tau_r (Cs 62 D_5/2) = 83.90 us   (was 77 us in level_parameters)
        tau_R (Rb 54 D_3/2) = 83.41 us   (was 74 us)
        tau_P (Rb 7 P_3/2)  =  0.131 us  (unchanged, paper value)

  2.  K refinement at every (a, Delta) point: K = 0.95 was tuned empirically
      at V_ct = 226 MHz (a = 5 um).  At V_ct = 442 MHz (a = 4 um) the
      protocol's optimal K may shift -- check by sweeping
      K in {0.92, 0.94, 0.95, 0.96, 0.97, 0.98, 1.00}.

Cases:
    A.  a = 5 um, Delta = 500 MHz, ARC lifetimes  -> updated FINAL
    B.  a = 4 um, Delta = 500 MHz, ARC lifetimes  -> "same Delta" test
    C.  a = 4 um, Delta = 250 MHz, ARC lifetimes  -> Delta retuned for M_2 ~ FINAL

For each case we report:
    - K refinement curve
    - Best F_bar and gate time
    - Per-input fidelities at the best K
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

# =====================================================================
# ARC-verified lifetimes (replace level_parameters' n^3 estimates)
# =====================================================================
tau_r_us = 83.90      # Cs 62 D_5/2  (ARC, 300K, BBR-incl, includeLevelsUpTo=92)
tau_R_us = 83.41      # Rb 54 D_3/2  (ARC, 300K, BBR-incl, includeLevelsUpTo=84)
tau_P_us = 0.131      # Rb 7 P_3/2   (paper, low-n; ARC agrees)
gamma_r  = 1.0 / tau_r_us
gamma_R  = 1.0 / tau_R_us
gamma_P  = 1.0 / tau_P_us

c_ops = build_collapse_operators(gamma_r, gamma_R, gamma_P)

# =====================================================================
# Pulse and geometry helpers
# =====================================================================
omega_p_MHz = 50.0
omega_R_MHz = 3.5 * omega_p_MHz
omega_c_MHz = 50.0
omega_p_amp = 2 * np.pi * omega_p_MHz
omega_R_amp = 2 * np.pi * omega_R_MHz
omega_c_amp = 2 * np.pi * omega_c_MHz
T_c = np.pi / omega_c_amp
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
    r_DD = a_um
    r_AA = a_um * math.sqrt(2)
    V_ct_MHz = 1000 * C3_tilde / r_DA**3
    V_DD_MHz = 1000 * C6_RbRb  / r_DD**6
    V_cc_MHz = 1000 * C6_CsCs  / r_AA**6
    return (V_ct_MHz, V_cc_MHz, V_DD_MHz, r_DA, r_AA)


def run_one(a_um, delta_MHz, K, ntime=350):
    V_ct_MHz, V_cc_MHz, V_DD_MHz, r_DA, r_AA = interactions_at(a_um)
    V_ct = 2 * np.pi * V_ct_MHz
    V_cc = 2 * np.pi * V_cc_MHz
    delta = 2 * np.pi * delta_MHz
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
    return F_bar, per_F, T_f * 1e3, t_total * 1e3, M2


def scan_K(a_um, delta_MHz, label):
    V_ct_MHz, V_cc_MHz, V_DD_MHz, r_DA, r_AA = interactions_at(a_um)
    print("=" * 80)
    print(f"{label}:  a = {a_um:.2f} um,  Delta = {delta_MHz:.0f} MHz")
    print(f"   V_ct = {V_ct_MHz:.2f} MHz,  V_cc = {V_cc_MHz:+.4f} MHz,  "
          f"R_AA = {V_ct_MHz/abs(V_cc_MHz):.1f}")
    print("=" * 80)
    Ks = [0.92, 0.94, 0.95, 0.96, 0.97, 0.98, 1.00]
    best_F = 0.0; best_K = None; best_perF = None; best_T = None; best_M2 = None
    print(f"{'K':>5s} {'gate(ns)':>10s} {'M2':>7s} {'F_bar':>10s} {'1-F':>10s}")
    print("-" * 50)
    for K in Ks:
        F, perF, Tfn, gn, M2 = run_one(a_um, delta_MHz, K)
        print(f"{K:>5.2f} {gn:>10.2f} {M2:>7.2f} {F:>10.6f} {1-F:>10.3e}")
        if F > best_F:
            best_F = F; best_K = K; best_perF = perF; best_T = gn; best_M2 = M2
    print(f"\nBest:  K = {best_K:.2f}, F_bar = {best_F:.6f}, "
          f"gate = {best_T:.2f} ns, M2 = {best_M2:.2f}")
    for (lbl, _, _), F in zip(basis_inputs, best_perF):
        print(f"   {lbl:<8s}  F = {F:.6f}")
    print()
    return {"a": a_um, "Delta": delta_MHz, "K": best_K, "F": best_F,
            "gate": best_T, "M2": best_M2, "perF": best_perF}


# =====================================================================
# Cases
# =====================================================================
print(f"Using ARC-verified lifetimes:  "
      f"tau_r = {tau_r_us:.2f} us, tau_R = {tau_R_us:.2f} us, "
      f"tau_P = {tau_P_us:.3f} us\n")

A = scan_K(5.00, 500.0, "Case A: FINAL transplanted with ARC lifetimes")
B = scan_K(4.00, 500.0, "Case B: a = 4 um, same Delta = 500 MHz")
C = scan_K(4.00, 250.0, "Case C: a = 4 um, Delta = 250 MHz (M_2 ~ FINAL)")

# =====================================================================
# Summary
# =====================================================================
print("=" * 80)
print("Summary  (all with ARC lifetimes tau_r = 83.90, tau_R = 83.41 us)")
print("=" * 80)
print(f"{'Case':<30s} {'a (um)':>8s} {'Delta':>8s} {'K':>6s} "
      f"{'gate(ns)':>10s} {'M2':>7s} {'F_bar':>10s} {'1-F':>10s}")
print("-" * 90)
for case, R in [("A: FINAL (a=5)", A),
                ("B: a=4, Delta=500", B),
                ("C: a=4, Delta=250", C)]:
    print(f"{case:<30s} {R['a']:>8.2f} {R['Delta']:>8.0f} {R['K']:>6.2f} "
          f"{R['gate']:>10.2f} {R['M2']:>7.2f} {R['F']:>10.6f} "
          f"{1-R['F']:>10.3e}")

# Decay-floor predictions
print()
print("Decay floor on |11,*>  (2 * T_gate / tau_r, with tau_r = 83.9 us):")
for case, R in [("A", A), ("B", B), ("C", C)]:
    floor = 2 * R["gate"] / 1000 / tau_r_us
    sim_inf_11 = 1 - R["perF"][6]
    print(f"   Case {case}:  T_gate = {R['gate']:6.2f} ns, "
          f"floor = {floor*100:5.2f}%,  simulated 1 - F_(11,*) = {sim_inf_11*100:5.2f}%, "
          f"residual = {(sim_inf_11 - floor)*100:5.2f}%")

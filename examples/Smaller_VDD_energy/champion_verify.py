"""
Detailed re-simulation of the brute-force champion configuration.

Most-robust cell found by brute_force_or.py:
    a            = 4.0 um
    Omega_p      = 2 pi * 60   MHz
    Omega_R      = 4.0 * Omega_p  =  2 pi * 240 MHz
    Delta        = 2 pi * 500  MHz
    Omega_c      = 2 pi * 70   MHz
    K            = 0.95
    alpha        = 4.0

Prints per-input fidelity and full population breakdown so the report can
quote the diagnostic table.
"""

import numpy as np

from triqg.atoms import CsAtom, RbAtom, composite_basis_state
from triqg.pulses import omega_gaussian, compute_pulse_area
from triqg.hamiltonian import build_hamiltonian
from triqg.decoherence import build_collapse_operators
from triqg.solver import simulate
from triqg.analysis import state_fidelity, average_gate_fidelity

from level_parameters import C3_tilde, C6_RbRb, C6_CsCs, gamma_r, gamma_R, gamma_P

# ---------------------------------------------------------------------
# Champion configuration (from brute_force_robust.csv, robust_score=1.0)
# ---------------------------------------------------------------------
a_um         = 4.00
omega_p_MHz  = 60.0
ratio        = 4.0
delta_MHz    = 500.0
omega_c_MHz  = 70.0

ALPHA = 4.0
K     = 0.95

# Derived geometry
r_DA = a_um / np.sqrt(2)
r_DD = a_um
r_AA = a_um * np.sqrt(2)
V_ct_MHz = 1000.0 * C3_tilde / r_DA ** 3
V_DD_MHz = 1000.0 * C6_RbRb  / r_DD ** 6
V_cc_MHz = 1000.0 * C6_CsCs  / r_AA ** 6
V_ct = 2 * np.pi * V_ct_MHz
V_cc = 2 * np.pi * V_cc_MHz

# Derived pulse parameters
omega_p_amp = 2 * np.pi * omega_p_MHz
omega_R_amp = 2 * np.pi * (ratio * omega_p_MHz)
omega_c_amp = 2 * np.pi * omega_c_MHz
delta       = 2 * np.pi * delta_MHz

T_c = np.pi / omega_c_amp
I_inf = 2 * 0.92770 * 2 ** (-1.0 / 6.0)
sigma_pi4 = ((2 * np.pi * delta) / (omega_p_amp ** 2 * I_inf)) ** 3
sigma = (K ** 3) * sigma_pi4
T_f   = (ALPHA * sigma) ** (1.0 / 3.0)
t_total = 2 * T_c + 2 * T_f

# Strong-blockade margins
one_photon_AC = omega_p_MHz ** 2 / (2 * delta_MHz)
two_photon_R  = omega_p_MHz * (ratio * omega_p_MHz) / (2 * delta_MHz)
M1 = V_ct_MHz / one_photon_AC
M2 = V_ct_MHz / two_photon_R

# ---------------------------------------------------------------------
# Verify pulse area
# ---------------------------------------------------------------------
args = {
    "omega_c_amp": omega_c_amp, "omega_p_amp": omega_p_amp,
    "omega_R_amp": omega_R_amp, "T_c": T_c, "T_f": T_f, "sigma": sigma,
}
area = compute_pulse_area(omega_gaussian, delta, T_c, T_c + 2 * T_f, args)
print(f"Effective two-photon pulse area: {area:.6f}  "
      f"(target: K*pi/4 = {K*np.pi/4:.6f}, K = {K})")

# ---------------------------------------------------------------------
# Header
# ---------------------------------------------------------------------
print("\nChampion configuration (brute_force_or.py rank #1)")
print("=" * 70)
print(f"  a        = {a_um:.2f} um")
print(f"  Omega_p  = 2 pi * {omega_p_MHz:.0f} MHz")
print(f"  Omega_R  = 2 pi * {ratio*omega_p_MHz:.0f} MHz  (ratio = {ratio})")
print(f"  Omega_c  = 2 pi * {omega_c_MHz:.0f} MHz  -> T_c = {T_c*1e3:.3f} ns")
print(f"  Delta    = 2 pi * {delta_MHz:.0f} MHz")
print(f"  K        = {K}  (area = {K}*pi/4)")
print(f"  alpha    = {ALPHA}  (smoothness)")
print(f"  sigma    = {sigma*1e3:.4f} ns")
print(f"  T_f      = {T_f*1e3:.3f} ns")
print(f"  Total    = {t_total*1e3:.2f} ns")
print()
print(f"  V_ct/(2pi) = {V_ct_MHz:>7.2f} MHz   (Rb-Cs Foerster)")
print(f"  V_cc/(2pi) = {V_cc_MHz:>7.4f} MHz   (Cs-Cs vdW at r_AA = {r_AA:.3f} um)")
print(f"  V_DD/(2pi) = {V_DD_MHz:>7.4f} MHz   (Rb-Rb vdW at r_DD = {a_um:.3f} um)")
print(f"  R_DD = V_ct / V_DD       = {V_ct_MHz/V_DD_MHz:>6.1f}")
print(f"  R_AA = V_ct / |V_cc|     = {V_ct_MHz/abs(V_cc_MHz):>6.1f}")
print()
print(f"  Strong-blockade margins (M >> 1 = safe):")
print(f"    M1 = V_ct / [Omega_p^2 / (2 Delta)]        = {M1:6.1f}")
print(f"    M2 = V_ct / [Omega_p Omega_R / (2 Delta)]  = {M2:6.1f}")

# ---------------------------------------------------------------------
# Hamiltonian + 8-input simulation
# ---------------------------------------------------------------------
H = build_hamiltonian(delta, V_ct, pulse_p=omega_gaussian, V_cc=V_cc)
c_ops = build_collapse_operators(gamma_r, gamma_R, gamma_P)

cs = CsAtom(); rb = RbAtom()
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
            label = f"|{c1},{c2},{['A','B'][t]}>"
            basis_inputs.append((label, c1, c2, t, psi_in, psi_id))

tlist = np.linspace(0, t_total, 500)
opts = {"store_final_state": True, "nsteps": 200000}

print()
print(f"{'Input':<14s}  {'F_k':>9}    final populations on target: A     B     P     R")
print("-" * 78)

fidelity_pairs = []
for label, c1, c2, t, psi_in, psi_id in basis_inputs:
    res = simulate(method="mesolve", H=H, psi0=psi_in, tlist=tlist,
                   c_ops=c_ops, e_ops=[], options=opts, args=args)
    F = state_fidelity(res.final_state, psi_id)
    fidelity_pairs.append((res.final_state, psi_id))

    pops = {}
    for lbl, idx in [("A", rb.level_index["A"]),
                     ("B", rb.level_index["B"]),
                     ("P", rb.level_index["P"]),
                     ("R", rb.level_index["R"])]:
        tgt = composite_basis_state(ctrl_levels[c1], ctrl_levels[c2], idx)
        pops[lbl] = state_fidelity(res.final_state, tgt)
    print(f"  {label:<11s} {F:>9.6f}    "
          f"{pops['A']:>6.4f} {pops['B']:>6.4f} "
          f"{pops['P']:>6.4f} {pops['R']:>6.4f}")

F_bar = average_gate_fidelity(fidelity_pairs)
print()
print(f"Average gate fidelity:  F_bar     = {F_bar:.6f}")
print(f"Gate infidelity:        1 - F_bar = {1 - F_bar:.3e}")

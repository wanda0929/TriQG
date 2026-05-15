"""
CCX (Toffoli) Average Gate Fidelity -- Smaller_V_DD redesign
========================================================================
Energy-level redesign:  Rb 54 D_{3/2} + Cs 62 D_{5/2}
    (replaces paper's Rb 66 D_{5/2} + Cs 76 D_{3/2})

Why this file exists
--------------------
The paper's Rb 66 D_5/2 + Cs 76 D_3/2 pair has a Rb-Rb same-species
blockade V_DD/(2 pi) ~ 91 MHz at r_DD = 5 um that was missed in
Appendix A.  Because the X-round Rb pulse is global, neighboring data
atoms blockade each other during the CCX correction, killing the
selectivity ratio R_DD = V_dd/V_DD down to ~5.7.  Moving Rb to
54 D_{3/2} parks the data atom on a same-species Foerster zero
crossing (|C_6^DD| -> 9.07 GHz um^6), lifting R_DD to ~389 at the
cost of a ~2.3x slower V_ct.

See:
    SelfCorrectingRydberg/plan/Rb_Cs_level_redesign.md (Sec. 4)
    examples/Smaller_VDD_energy/level_parameters.py

Pulse-sequence revision (this file)
-----------------------------------
Relative to examples/Average_fidelity/ccx_average_gate_fidelity.py
(Omega_cc = 2 pi * 100, Omega_t = 2 pi * 50):

  * Omega_cc_amp = 2 pi * 50  MHz   (was 2 pi * 100)
  * Omega_t_amp  = 2 pi * 20  MHz   (was 2 pi * 50)

The lower amplitudes lengthen each sub-pulse:
  T_cc = pi / Omega_cc =  10 ns  (was  5 ns)
  T_t  = pi / Omega_t  =  25 ns  (was 10 ns)

Total CCX time:  2*T_cc + 3*T_t  =  95 ns   (was 40 ns)

The motivation is the same as for the OR-gate revision:  smaller drive
amplitudes sit further inside the blockade window
(V_ct / Omega_t = 226/20 = 11.3 in this redesign vs. 517/50 = 10.3
in the paper), keeping leakage to the unwanted |R> branch suppressed
even though V_ct is now 2.3x smaller.

Truth table:
    |c1, c2, t>   ->   |c1, c2, t XOR (c1 AND c2)>
    target flips only when BOTH controls are |1>.
"""

import numpy as np

from triqg.atoms import CsAtom, RbAtom, composite_basis_state
from triqg.pulses import omega_cc, omega_t1, omega_t2
from triqg.hamiltonian import build_ccx_hamiltonian
from triqg.decoherence import build_collapse_operators
from triqg.solver import simulate
from triqg.analysis import state_fidelity, average_gate_fidelity

from level_parameters import (
    a_um, r_DA, r_AA,
    C3_tilde, C6_RbRb, C6_CsCs, Delta_F_MHz,
    V_ct, V_cc, V_ct_MHz, V_cc_MHz, V_DD_MHz,
    R_DD, R_AA,
    gamma_r, gamma_R, gamma_P,
    tau_r_us, tau_R_us, tau_P_us,
)

# =====================================================================
# Pulse-sequence parameters (REVISED for the smaller-V_DD redesign)
# =====================================================================
omega_cc_amp = 2 * np.pi * 100    # Cs control pi-pulse Rabi frequency [rad/us]
omega_t_amp  = 2 * np.pi * 40    # Rb target sub-pulse Rabi frequency [rad/us]

T_cc = np.pi / omega_cc_amp      # Control pi-pulse duration  = 10.0 ns
T_t  = np.pi / omega_t_amp       # Target  sub-pulse duration = 25.0 ns

args = {
    "omega_cc_amp": omega_cc_amp,
    "omega_t_amp": omega_t_amp,
    "T_cc": T_cc,
    "T_t": T_t,
}

# =====================================================================
# Build Hamiltonian and collapse operators
# =====================================================================
# V_ct, V_cc imported from level_parameters (angular MHz = rad/us).
# V_DD does not enter the 3-atom Hamiltonian -- it is a cross-gate
# global-drive effect.
H = build_ccx_hamiltonian(V_ct, V_cc=V_cc)
c_ops = build_collapse_operators(gamma_r, gamma_R, gamma_P)

# =====================================================================
# Time list and solver options
# =====================================================================
t_total = 2 * T_cc + 3 * T_t

# Step bound: keep 20 samples per drive period of the FASTEST tone.
max_freq = max(omega_cc_amp, omega_t_amp) / (2 * np.pi)
max_step = 1.0 / (20 * max_freq)

tlist = np.linspace(0, t_total, 500)
solver_opts = {
    "store_final_state": True,
    "nsteps": 100000,
    "max_step": max_step,
}

# =====================================================================
# Enumerate all 2^{n+1} = 8 computational basis states
# =====================================================================
cs = CsAtom()
rb = RbAtom()

n_controls = 2
ctrl_levels = [cs.level_index["0"], cs.level_index["1"]]
tgt_levels = [rb.level_index["A"], rb.level_index["B"]]


def ccx_ideal_output(c1_bit: int, c2_bit: int, t_bit: int):
    """Ideal CCX output: target flips iff c1 AND c2 are both 1."""
    if c1_bit == 1 and c2_bit == 1:
        t_out = 1 - t_bit
    else:
        t_out = t_bit
    return composite_basis_state(
        ctrl_levels[c1_bit], ctrl_levels[c2_bit], tgt_levels[t_out]
    )


basis_inputs = []
for c1_bit in range(2):
    for c2_bit in range(2):
        for t_bit in range(2):
            psi_in = composite_basis_state(
                ctrl_levels[c1_bit], ctrl_levels[c2_bit], tgt_levels[t_bit]
            )
            psi_ideal = ccx_ideal_output(c1_bit, c2_bit, t_bit)
            basis_inputs.append(
                (str(c1_bit), str(c2_bit), ["A", "B"][t_bit], psi_in, psi_ideal)
            )

# =====================================================================
# Header
# =====================================================================
print(f"CCX Average Gate Fidelity -- Smaller_V_DD redesign")
print(f"Energy levels: Rb 54 D_3/2 + Cs 62 D_5/2,  a = {a_um:.2f} um")
print(f"==================================================================")
print(f"n_controls = {n_controls},  basis states = {len(basis_inputs)}")
print(f"Pulse amplitudes:")
print(f"  Omega_cc / (2 pi) = {omega_cc_amp/(2*np.pi):.1f} MHz  -> T_cc = {T_cc*1e3:.4f} ns")
print(f"  Omega_t  / (2 pi) = {omega_t_amp /(2*np.pi):.1f} MHz  -> T_t  = {T_t *1e3:.4f} ns")
print(f"Interactions:")
print(f"  V_ct / (2 pi) = +{V_ct_MHz:7.3f} MHz   (Rb-Cs Foerster, |54D3/2;62D5/2>)")
print(f"  V_cc / (2 pi) = {V_cc_MHz:+8.4f} MHz   (Cs-Cs vdW at r_AA = {r_AA:.3f} um)")
print(f"  V_DD / (2 pi) = {V_DD_MHz:+8.4f} MHz   (Rb-Rb vdW at r_DD = {a_um:.3f} um) [info only]")
print(f"  Foerster defect Delta_F / (2 pi) = {Delta_F_MHz:+.2f} MHz")
print(f"  R_DD = {R_DD:.1f},  R_AA = {R_AA:.1f}   (both >= 100 OK)")
print(f"Ratios:")
print(f"  V_ct / Omega_t   = {V_ct_MHz/(omega_t_amp/(2*np.pi)):.3f}   "
      f"(target >> 1 for blockade)")
print(f"  |V_cc| / Omega_cc = {abs(V_cc_MHz)/(omega_cc_amp/(2*np.pi)):.4f}   "
      f"(target << 1 for clean controls)")
print(f"Decoherence (re-verify with ARC):")
print(f"  tau_r (Cs 62 D_5/2) = {tau_r_us:.2f} us")
print(f"  tau_R (Rb 54 D_3/2) = {tau_R_us:.2f} us")
print(f"  tau_P (Rb  7 P_3/2) = {tau_P_us:.3f} us")
print(f"Total CCX time = {t_total*1e3:.3f} ns  ({t_total:.4f} us)\n")

# =====================================================================
# Run mesolve for each input and collect fidelities
# =====================================================================
fidelity_pairs = []

for i, (c1_lbl, c2_lbl, t_lbl, psi_in, psi_ideal) in enumerate(basis_inputs):
    label = f"|{c1_lbl},{c2_lbl},{t_lbl}>"
    print(f"  [{i + 1}/8] Simulating input {label} ...", end=" ", flush=True)

    result = simulate(
        method="mesolve",
        H=H,
        psi0=psi_in,
        tlist=tlist,
        c_ops=c_ops,
        e_ops=[],
        options=solver_opts,
        args=args,
    )

    fid = state_fidelity(result.final_state, psi_ideal)
    fidelity_pairs.append((result.final_state, psi_ideal))
    print(f"F = {fid:.6f}")

# =====================================================================
# Compute average gate fidelity
# =====================================================================
F_bar = average_gate_fidelity(fidelity_pairs)

print(f"\n{'=' * 60}")
print(f"Average gate fidelity:  F_bar = {F_bar:.6f}")
print(f"Gate infidelity:        1 - F_bar = {1 - F_bar:.3e}")
print(f"{'=' * 60}")

# =====================================================================
# Diagnostic: final-state population breakdown
# =====================================================================
print("\nDiagnostic: final-state populations on each Rb target level")
print("-" * 70)
print(f"  {'Input':<12} {'P(A)':>9} {'P(B)':>9} {'P(P)':>9} {'P(R)':>9}")
print("-" * 70)

for i, (c1_lbl, c2_lbl, t_lbl, psi_in, psi_ideal) in enumerate(basis_inputs):
    rho_out = fidelity_pairs[i][0]
    pops = {}
    for lbl, idx in [
        ("A", rb.level_index["A"]),
        ("B", rb.level_index["B"]),
        ("P", rb.level_index["P"]),
        ("R", rb.level_index["R"]),
    ]:
        tgt = composite_basis_state(
            ctrl_levels[int(c1_lbl)], ctrl_levels[int(c2_lbl)], idx
        )
        pops[lbl] = state_fidelity(rho_out, tgt)
    label = f"|{c1_lbl},{c2_lbl},{t_lbl}>"
    print(
        f"  {label:<12} {pops['A']:>9.5f} {pops['B']:>9.5f} "
        f"{pops['P']:>9.5f} {pops['R']:>9.5f}"
    )

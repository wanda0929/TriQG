"""
OR Gate Average Gate Fidelity via Eq. (7) of Yu et al. (Gaussian Pulse)
------------------------------------------------------------------------
OPTION A: Rb 66 D_{5/2} + Cs 76 D_{3/2} d-state Foerster pair
          at a = 5.00 um,  Omega_R = 2.9 * Omega_p
------------------------------------------------------------------------

Computes the average gate fidelity of the three-qubit Rydberg OR gate
following the method in:

    D. Yu et al., "Multiqubit Toffoli gates and optimal geometry with
    Rydberg atoms", arXiv:2203.14302v2, Eq. (7).

    F_bar = (1 / 2^{n+1}) * sum_k  F(rho_out^(k), rho_et^(k))

The Rydberg pair is the 'Candidate #1' d-state Foerster pair from:

    B. J. Ireland, J. D. Pritchard, J. P. Shaffer,
    "Interspecies Foerster resonances of Rb-Cs Rydberg d-states
    for enhanced multi-qubit gate fidelities",
    Phys. Rev. Research 6, 013293 (2024), arXiv:2401.02308.

Configuration:
  * Rydberg levels:
        Rb |R> = |66 D_{5/2}>
        Cs |r> = |76 D_{3/2}>
        Rb intermediate |P> = |7 P_{3/2}>
  * Foerster channel (Ireland et al., Table I, row #1):
        |66 D_{5/2} ; 76 D_{3/2}>  <->  |67 P_{3/2} ; 74 F_{5/2}>
  * C_3, C_6 taken from Ireland et al. for this pair.
  * Lattice a = 5.00 um.
  * V_cc derived from C_6 and a_um (no override).
  * Omega_R / Omega_p = 2.9.
  * Rydberg lifetimes from ARC at T = 300 K (see below).

Interaction strengths at a = 5.00 um:
    V_ct / (2 pi) = +516.81 MHz
    V_cc / (2 pi) =   -5.54 MHz

The OR gate protocol (Farouk et al.):
    1. Excite controls in |1> to Rydberg |r> (positive pi-pulse)
    2. Two-photon Raman pulse on target via |P> (detuned by delta)
    3. De-excite controls (negative pi-pulse)

The target undergoes the two-photon process |A> -> |P> -> |R> -> |P> -> |A>
(a full Rabi cycle) and returns to its original computational state,
acquiring only a phase. Under blockade (any control in |r>), the |R>
level is shifted by V_ct, breaking the two-photon resonance -- the
target evolution is suppressed and it stays unchanged (no phase).

When both controls are |0>, the target completes the full cycle and
remains unchanged. When at least one control is |1> (OR condition),
the blockade prevents the return cycle and the target is flipped.

Truth table (computational populations):

    Input       Ideal output
    |0,0,A>  -> |0,0,A>   (both controls 0: target unchanged)
    |0,0,B>  -> |0,0,B>
    |0,1,A>  -> |0,1,B>   (one control in 1: target flips)
    |0,1,B>  -> |0,1,A>
    |1,0,A>  -> |1,0,B>   (one control in 1: target flips)
    |1,0,B>  -> |1,0,A>
    |1,1,A>  -> |1,1,B>   (both controls in 1: target flips)
    |1,1,B>  -> |1,1,A>

Outputs:
    Per-input state fidelities and the average gate fidelity F_bar.
"""

import numpy as np

from triqg.atoms import CsAtom, RbAtom, composite_basis_state
from triqg.pulses import omega_c, omega_gaussian, omega_R, compute_pulse_area
from triqg.hamiltonian import build_hamiltonian
from triqg.decoherence import build_collapse_operators
from triqg.solver import simulate
from triqg.analysis import state_fidelity, average_gate_fidelity

# =====================================================================
# Physical parameters (from or_gate_mesolve_gaussian.py)
# =====================================================================
# Pulse and interaction parameters from the SelfCorrectingRydberg paper
# (main.tex, Sec. III.A "Three-qubit OR gate (EIT + Rydberg blockade)").
omega_c_amp = 2 * np.pi * 50     # Cs control Rabi frequency Omega_c [MHz]
omega_p_amp = 2 * np.pi * 50.0   # Rb target two-photon probe amplitude
omega_R_amp = 2.9 * omega_p_amp  # Omega_R = 2.9 * Omega_p

delta = 2 * np.pi * 500  # Two-photon detuning Delta [MHz]

# ---------------------------------------------------------------------
# Lattice geometry
# ---------------------------------------------------------------------
a_um = 5.0                       # rotated-lattice spacing [um]
r_DA = a_um / np.sqrt(2)          # nearest data-ancilla distance ~ 3.5355 um
r_AA = a_um * np.sqrt(2)          # nearest same-type ancilla-ancilla ~ 7.0711 um

# ---------------------------------------------------------------------
# Rydberg interaction coefficients for the OPTION A Foerster pair
# Rb |66 D_{5/2}> + Cs |76 D_{3/2}>
# at theta = 90 deg (quantization axis perpendicular to the array).
#
# Source: Ireland, Pritchard & Shaffer, Phys. Rev. Research 6, 013293 (2024),
#         Table I row 1 (candidate #1), also tabulated locally in
#         reference/energy_level_inter.md.
#
# Foerster channel:
#     |Rb 66 D_{5/2}; Cs 76 D_{3/2}>  <->  |Rb 67 P_{3/2}; Cs 74 F_{5/2}>
# Blockade fidelity (single-control, Ireland et al.): P_1r = 0.9997
# ---------------------------------------------------------------------
C3_tilde = 22.84                  # Effective Forster C_3 [GHz * um^3]
C6_CsCs = -692.9                  # Cs-Cs vdW C_6 [GHz * um^6], SIGNED

# Derived blockade strengths (in angular MHz = rad / us since t is in us).
V_ct_MHz = 1000.0 * C3_tilde / r_DA**3   # ~ +516.92 MHz at a = 5.00 um
V_cc_MHz = 1000.0 * C6_CsCs  / r_AA**6   # ~   -5.543 MHz at a = 5.00 um
V_ct = 2 * np.pi * V_ct_MHz              # Rb-Cs dipole-dipole blockade
V_cc = 2 * np.pi * V_cc_MHz              # Cs-Cs van der Waals (signed; C_6 < 0)

T_c = np.pi / omega_c_amp  # Control pi-pulse duration = 10 ns (0.010 us)
T_f = 0.15                 # Target pulse half-window T_f = 150 ns
# Super-Gaussian width chosen so the effective two-photon pulse area
# integral Omega_p^2 / (2 Delta) dt equals pi/4.
sigma = 0.001771

# ---------------------------------------------------------------------
# Decoherence rates for the OPTION A level choice at T = 300 K,
# including blackbody radiation.
#
# Sources for the two Rydberg-state lifetimes:
#   [1] N. Sibalic, J. D. Pritchard, C. S. Adams, K. J. Weatherill,
#       "ARC: An open-source library for calculating properties of
#       alkali Rydberg atoms", Comp. Phys. Commun. 220, 319 (2017),
#       arXiv:1612.05529. ARC implements the
#   [2] I. I. Beterov, I. I. Ryabtsev, D. B. Tretyakov, V. M. Entin,
#       "Quasiclassical calculations of blackbody-radiation-induced
#       depopulation rates and effective lifetimes of Rydberg nS, nP,
#       and nD alkali-metal atoms with n <= 80",
#       Phys. Rev. A 79, 052504 (2009), arXiv:0902.4995
#   BBR formulas, combined with radiative rates from Einstein A
#   coefficients and standard quantum defects.
#
# Computed via:  atom.getStateLifetime(n, l, j, temperature=300,
#                                       includeLevelsUpTo=n+30)
#
# Values at T = 300 K from ARC 3.10.2:
#   Rb 66 D_{5/2}:  tau_rad = 294.31 us,  tau_BBR = 248.96 us,
#                   tau_total = 134.87 us      <-- used here
#   Cs 76 D_{3/2}:  tau_rad = 260.15 us,  tau_BBR = 316.22 us,
#                   tau_total = 142.73 us      <-- used here
#
# The Rb intermediate state |P> = |7 P_{3/2}> keeps the paper value
# tau_P = 0.131 us.
# Time unit throughout this script is microseconds.
gamma_r = 1.0 / 142.73   # Cs |r> = |76 D_{3/2}>, tau_r = 142.73 us (ARC, T=300 K)
gamma_R = 1.0 / 134.87   # Rb |R> = |66 D_{5/2}>, tau_R = 134.87 us (ARC, T=300 K)
gamma_P = 1.0 / 0.131    # Rb |P> = |7 P_{3/2}>,  tau_P = 0.131 us (paper)

args = {
    "omega_c_amp": omega_c_amp,
    "omega_p_amp": omega_p_amp,
    "omega_R_amp": omega_R_amp,
    "T_c": T_c,
    "T_f": T_f,
    "sigma": sigma,
}

# =====================================================================
# Verify pulse area
# =====================================================================
area = compute_pulse_area(omega_gaussian, delta, T_c, T_c + 2 * T_f, args)
print(f"Effective two-photon pulse area: {area:.4f} (target: pi = {np.pi:.4f})")

# =====================================================================
# Build Hamiltonian and collapse operators
# =====================================================================
H = build_hamiltonian(delta, V_ct, pulse_p=omega_gaussian, V_cc=V_cc)
c_ops = build_collapse_operators(gamma_r, gamma_R, gamma_P)

# =====================================================================
# Time list and solver options
# =====================================================================
t_total = 2 * T_c + 2 * T_f
tlist = np.linspace(0, t_total, 500)

solver_opts = {"store_final_state": True, "nsteps": 10000}

# =====================================================================
# Enumerate all 2^{n+1} = 8 computational basis states
# =====================================================================
cs = CsAtom()
rb = RbAtom()

n_controls = 2
ctrl_levels = [cs.level_index["0"], cs.level_index["1"]]
tgt_levels = [rb.level_index["A"], rb.level_index["B"]]


def or_gate_ideal_output(c1_bit: int, c2_bit: int, t_bit: int):
    """
    Return the ideal OR gate output state for given input bits.

    When both controls are in state |0>, the target is unchanged.
    When at least one control is in state |1> (OR condition), the
    target is flipped (A <-> B).
    """
    if c1_bit or c2_bit:
        t_out = 1 - t_bit  # target flips when OR condition is met
    else:
        t_out = t_bit  # target unchanged when both controls are 0
    return composite_basis_state(
        ctrl_levels[c1_bit], ctrl_levels[c2_bit], tgt_levels[t_out]
    )


# Build input/output table
basis_inputs = []
for c1_bit in range(2):
    for c2_bit in range(2):
        for t_bit in range(2):
            psi_in = composite_basis_state(
                ctrl_levels[c1_bit], ctrl_levels[c2_bit], tgt_levels[t_bit]
            )
            psi_ideal = or_gate_ideal_output(c1_bit, c2_bit, t_bit)

            c1_lbl = str(c1_bit)
            c2_lbl = str(c2_bit)
            t_lbl = ["A", "B"][t_bit]
            basis_inputs.append((c1_lbl, c2_lbl, t_lbl, psi_in, psi_ideal))

# =====================================================================
# Run mesolve for each input and collect fidelities
# =====================================================================
print(f"\nOR Gate Average Gate Fidelity [Yu et al., Eq. (7)] - Gaussian Pulse")
print(f"Option A: Rb 66 D_5/2 + Cs 76 D_3/2 at a = {a_um} um,  Omega_R/Omega_p = 2.9")
print(f"====================================================================")
print(f"n_controls = {n_controls},  basis states = {len(basis_inputs)}")
print(f"T_c = {T_c * 1e3:.3f} ns,  T_f = {T_f * 1e3:.3f} ns")
print(f"sigma = {sigma},  delta = {delta / (2 * np.pi):.1f} MHz")
print(f"Omega_R/Omega_p = {omega_R_amp / omega_p_amp:.2f}")
print(f"Lattice:  a = {a_um} um, r_DA = {r_DA:.4f} um, r_AA = {r_AA:.4f} um")
print(f"Foerster pair: Rb 66 D_5/2 + Cs 76 D_3/2 (Ireland et al. 2024)")
print(f"  C3_tilde = {C3_tilde} GHz*um^3,  C6_CsCs = {C6_CsCs} GHz*um^6")
print(f"V_ct/(2pi) = +{V_ct_MHz:7.3f} MHz  (Rb-Cs d-state Forster, derived)")
print(f"V_cc/(2pi) = {V_cc_MHz:+7.3f} MHz  (Cs-Cs vdW, derived)")
print(f"V_ct/Delta    = {V_ct_MHz/500.0:.3f}")
print(f"|V_cc|/Omega_c = {abs(V_cc_MHz)/50.0:.3f}")
print(f"Lifetimes (ARC, T=300 K):  tau_r = 142.73 us,  tau_R = 134.87 us,  tau_P = 0.131 us (paper)")
print(f"Total gate time = {t_total * 1e3:.3f} ns\n")

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

print(f"\n{'=' * 50}")
print(f"Average gate fidelity:  F_bar = {F_bar:.6f}")
print(f"Gate infidelity:        1 - F_bar = {1 - F_bar:.2e}")
print(f"{'=' * 50}")

# =====================================================================
# Diagnostic: population breakdown for each output
# =====================================================================
print("\nDiagnostic: output population breakdown")
print("-" * 70)
print(f"  {'Input':<12} {'P(A)':>8} {'P(B)':>8} {'P(P)':>8} {'P(R)':>8}")
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
        f"  {label:<12} {pops['A']:>8.4f} {pops['B']:>8.4f} "
        f"{pops['P']:>8.4f} {pops['R']:>8.4f}"
    )

"""
CCX (Toffoli) Average Gate Fidelity via Eq. (7) of Yu et al.
------------------------------------------------------------------------
OPTION A: Rb 66 D_{5/2} + Cs 76 D_{3/2} d-state Foerster pair
          at a = 5.00 um
------------------------------------------------------------------------

The average gate fidelity for an n-control Toffoli gate is:
F_bar_n = (1 / 2^{n+1}) * sum_k Tr{ [sqrt(rho_et) |Psi_out><Psi_out| sqrt(rho_et)]^{1/2} }
=============================================================

The Rydberg pair is the 'Candidate #1' d-state Foerster pair from:

    B. J. Ireland, J. D. Pritchard, J. P. Shaffer,
    "Interspecies Foerster resonances of Rb-Cs Rydberg d-states for
    enhanced multi-qubit gate fidelities",
    Phys. Rev. Research 6, 013293 (2024), arXiv:2401.02308,
    Table I row #1 (candidate #1).

Configuration:
  * Rb |R> = |66 D_{5/2}>
  * Cs |r> = |76 D_{3/2}>
  * Foerster channel: |66 D_{5/2}; 76 D_{3/2}> <-> |67 P_{3/2}; 74 F_{5/2}>
  * C_3, C_6 from Ireland et al. 2024 for this pair
  * Lattice a = 5.00 um
  * V_cc derived from C_6 and a_um (no override)
  * Rydberg lifetimes from ARC at T = 300 K

Interaction strengths at a = 5.00 um:
    V_ct / (2 pi) = +516.81 MHz
    V_cc / (2 pi) =   -5.54 MHz

Computes the average gate fidelity of the three-qubit Rydberg CCX gate
following the method in:

    D. Yu et al., "Multiqubit Toffoli gates and optimal geometry with
    Rydberg atoms", arXiv:2203.14302v2, Eq. (7).

The average gate fidelity is defined as:

    F_bar = (1 / 2^{n+1}) * sum_k  F(rho_out^(k), rho_et^(k))

where:
    - n = 2 is the number of control qubits
    - 2^{n+1} = 8 is the total number of computational basis inputs
    - rho_out^(k) is the simulated output for the k-th basis input
    - rho_et^(k)  is the ideal CCX output for that input
    - F is the quantum state fidelity (Uhlmann fidelity squared)

For the CCX gate the target qubit flips only when BOTH controls are |1>:

    |c1, c2, t>   ->   |c1, c2, t XOR (c1 AND c2)>

    Input       Ideal output
    |0,0,A>  -> |0,0,A>
    |0,0,B>  -> |0,0,B>
    |0,1,A>  -> |0,1,A>
    |0,1,B>  -> |0,1,B>
    |1,0,A>  -> |1,0,A>
    |1,0,B>  -> |1,0,B>
    |1,1,A>  -> |1,1,B>   (flipped)
    |1,1,B>  -> |1,1,A>   (flipped)

Outputs:
    Per-input state fidelities and the average gate fidelity F_bar.
"""

import numpy as np

from triqg.atoms import CsAtom, RbAtom, composite_basis_state
from triqg.pulses import omega_cc, omega_t1, omega_t2
from triqg.hamiltonian import build_ccx_hamiltonian
from triqg.decoherence import build_collapse_operators
from triqg.solver import simulate
from triqg.analysis import state_fidelity, average_gate_fidelity

# =====================================================================
# Physical parameters (same as ccx_gate_mesolve.py)
# =====================================================================
# Pulse and interaction parameters from the SelfCorrectingRydberg paper
# (main.tex, Sec. III.B "Three-qubit CCX (Toffoli) gate").
omega_cc_amp = 2 * np.pi * 100  # Cs control pi-pulse Rabi frequency [MHz]
omega_t_amp = 2 * np.pi * 50  # Rb target sub-pulse Rabi frequency [MHz]

T_cc = np.pi / omega_cc_amp  # Control pi-pulse duration  (= 5 ns)
T_t = np.pi / omega_t_amp  # Target sub-pulse duration (= 10 ns)

# ---------------------------------------------------------------------
# Lattice geometry
# ---------------------------------------------------------------------
a_um = 5.00                       # rotated-lattice spacing [um]
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
#   BBR formulas combined with radiative rates from Einstein A
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
# ---------------------------------------------------------------------
gamma_r = 1.0 / 142.73   # Cs |r> = |76 D_{3/2}>, tau_r = 142.73 us (ARC, T=300 K)
gamma_R = 1.0 / 134.87   # Rb |R> = |66 D_{5/2}>, tau_R = 134.87 us (ARC, T=300 K)
gamma_P = 1.0 / 0.131    # Rb |P> = |7 P_{3/2}>,  tau_P = 0.131 us (paper)

args = {
    "omega_cc_amp": omega_cc_amp,
    "omega_t_amp": omega_t_amp,
    "T_cc": T_cc,
    "T_t": T_t,
}

# =====================================================================
# Build Hamiltonian and collapse operators
# =====================================================================
H = build_ccx_hamiltonian(V_ct, V_cc=V_cc)
c_ops = build_collapse_operators(gamma_r, gamma_R, gamma_P)

# =====================================================================
# Time list and solver options
# =====================================================================
t_total = 2 * T_cc + 3 * T_t
tlist = np.linspace(0, t_total, 500)

max_freq = max(omega_cc_amp, omega_t_amp) / (2 * np.pi)
max_step = 1.0 / (20 * max_freq)

solver_opts = {"store_final_state": True, "nsteps": 100000, "max_step": max_step}

# =====================================================================
# Enumerate all 2^{n+1} = 8 computational basis states
# =====================================================================
cs = CsAtom()
rb = RbAtom()

n_controls = 2
# Computational levels for controls: |0> and |1>
ctrl_levels = [cs.level_index["0"], cs.level_index["1"]]
# Computational levels for target:   |A> (logical 0) and |B> (logical 1)
tgt_levels = [rb.level_index["A"], rb.level_index["B"]]


def ccx_ideal_output(c1_bit: int, c2_bit: int, t_bit: int):
    """Return the ideal CCX output state for given input bits."""
    if c1_bit == 1 and c2_bit == 1:
        t_out = 1 - t_bit  # flip target
    else:
        t_out = t_bit  # target unchanged
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
            psi_ideal = ccx_ideal_output(c1_bit, c2_bit, t_bit)

            c1_lbl = str(c1_bit)
            c2_lbl = str(c2_bit)
            t_lbl = ["A", "B"][t_bit]
            basis_inputs.append((c1_lbl, c2_lbl, t_lbl, psi_in, psi_ideal))

# =====================================================================
# Run mesolve for each input and collect fidelities
# =====================================================================
print(f"CCX Average Gate Fidelity [Yu et al., Eq. (7)]")
print(f"Option A: Rb 66 D_5/2 + Cs 76 D_3/2 at a = {a_um} um")
print(f"==============================================")
print(f"n_controls = {n_controls},  basis states = {len(basis_inputs)}")
print(f"T_cc = {T_cc * 1e3:.3f} ns,  T_t = {T_t * 1e3:.3f} ns")
print(f"Lattice:  a = {a_um} um, r_DA = {r_DA:.4f} um, r_AA = {r_AA:.4f} um")
print(f"Foerster pair: Rb 66 D_5/2 + Cs 76 D_3/2 (Ireland et al. 2024)")
print(f"  C3_tilde = {C3_tilde} GHz*um^3,  C6_CsCs = {C6_CsCs} GHz*um^6")
print(f"V_ct/(2pi) = +{V_ct_MHz:7.3f} MHz  (Rb-Cs d-state Forster, derived)")
print(f"V_cc/(2pi) = {V_cc_MHz:+7.3f} MHz  (Cs-Cs vdW, derived)")
print(f"V_ct/Omega_t   = {V_ct_MHz/50.0:.3f}")
print(f"|V_cc|/Omega_cc = {abs(V_cc_MHz)/100.0:.3f}")
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

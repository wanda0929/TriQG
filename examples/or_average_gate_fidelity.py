"""
OR Gate Average Gate Fidelity via Eq. (7) of Yu et al.
======================================================

Computes the average gate fidelity of the three-qubit Rydberg OR gate
following the method in:

    D. Yu et al., "Multiqubit Toffoli gates and optimal geometry with
    Rydberg atoms", arXiv:2203.14302v2, Eq. (7).

    F_bar = (1 / 2^{n+1}) * sum_k  F(rho_out^(k), rho_et^(k))

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
from triqg.pulses import omega_c, omega_p, omega_R, compute_pulse_area
from triqg.hamiltonian import build_hamiltonian
from triqg.decoherence import build_collapse_operators
from triqg.solver import simulate
from triqg.analysis import state_fidelity, average_gate_fidelity

# =====================================================================
# Physical parameters (same as or_gate_mesolve.py)
# =====================================================================
omega_c_amp = 2 * np.pi * 50  # Control pulse Rabi frequency [MHz]
omega_p_amp = 2 * np.pi * 70  # Target probe pulse amplitude [MHz]
omega_R_amp = 2.5 * omega_p_amp  # Target Rydberg coupling amplitude [MHz]

delta = 2 * np.pi * 1200  # Detuning [MHz]
V_ct = 2 * np.pi * 500  # Rydberg blockade strength [MHz]

T_c = np.pi / omega_c_amp  # Control pi-pulse duration [us]
T_f = 0.32653  # Target pulse window [us]

gamma_r = 1.0 / 548.0  # Cs |r> decay rate [MHz]
gamma_R = 1.0 / 505.0  # Rb |R> decay rate [MHz]
gamma_P = 1.0 / 0.131  # Rb |P> decay rate [MHz]

args = {
    "omega_c_amp": omega_c_amp,
    "omega_p_amp": omega_p_amp,
    "omega_R_amp": omega_R_amp,
    "T_c": T_c,
    "T_f": T_f,
}

# =====================================================================
# Verify pulse area
# =====================================================================
area = compute_pulse_area(omega_p, delta, T_c, T_c + 2 * T_f, args)
print(f"Effective two-photon pulse area: {area:.4f} (target: pi = {np.pi:.4f})")

# =====================================================================
# Build Hamiltonian and collapse operators
# =====================================================================
H = build_hamiltonian(delta, V_ct)
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
print(f"\nOR Gate Average Gate Fidelity [Yu et al., Eq. (7)]")
print(f"==================================================")
print(f"n_controls = {n_controls},  basis states = {len(basis_inputs)}")
print(f"T_c = {T_c * 1e3:.3f} ns,  T_f = {T_f * 1e3:.3f} ns")
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

"""
CCX (Toffoli) Average Gate Fidelity via Eq. (7) of Yu et al.
The average gate fidelity for an n-control Toffoli gate is:
F_bar_n = (1 / 2^{n+1}) * sum_k Tr{ [sqrt(rho_et) |Psi_out><Psi_out| sqrt(rho_et)]^{1/2} }
=============================================================

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

V_ct = 2 * np.pi * 593  # Rb-Cs Rydberg blockade strength [MHz] (Table I)

# Decoherence rates from main.tex, Sec. III.C "Decoherence channels and
# gate fidelity". Lifetimes are quoted at T = 300 K including blackbody
# radiation. Time throughout this script is in microseconds, so each
# rate is 1/tau with tau in us.
gamma_r = 1.0 / 340.0  # Cs |r> = |79 D_{5/2}>, tau_r ~ 340 us
gamma_R = 1.0 / 260.0  # Rb |R> = |69 D_{5/2}>, tau_R ~ 260 us
gamma_P = 1.0 / 0.131  # Rb |P> = |7 P_{3/2}>,  tau_P = 0.131 us

args = {
    "omega_cc_amp": omega_cc_amp,
    "omega_t_amp": omega_t_amp,
    "T_cc": T_cc,
    "T_t": T_t,
}

# =====================================================================
# Build Hamiltonian and collapse operators
# =====================================================================
H = build_ccx_hamiltonian(V_ct)
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
print(f"==============================================")
print(f"n_controls = {n_controls},  basis states = {len(basis_inputs)}")
print(f"T_cc = {T_cc * 1e3:.3f} ns,  T_t = {T_t * 1e3:.3f} ns")
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

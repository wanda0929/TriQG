"""
OR Gate Simulation using mcsolve (Monte Carlo Quantum Trajectories)
===================================================================

Simulates the three-qubit Rydberg OR gate with two 133-Cs control atoms
and one 87-Rb target atom. Uses QuTiP's mcsolve via the triqg library.

Reference: M. Farouk et al., "Parallel Implementation of CNOTN and C2NOT2
Gates via Homonuclear and Heteronuclear Förster Interactions of Rydberg Atoms."

Outputs:
  1. Pulse envelope plot
  2. Averaged population dynamics with ±1σ shaded bands
  3. Quantum jump statistics
  4. State fidelity from the averaged final state
"""

import numpy as np
import matplotlib

matplotlib.use("Agg")  # Use non-interactive backend; remove for interactive use
import matplotlib.pyplot as plt

from triqg.atoms import CsAtom, RbAtom, composite_basis_state, composite_projector
from triqg.pulses import omega_c, omega_p, omega_R, compute_pulse_area
from triqg.hamiltonian import build_hamiltonian
from triqg.decoherence import build_collapse_operators
from triqg.solver import simulate
from triqg.analysis import state_fidelity, extract_populations
from triqg.visualization import plot_pulses, plot_populations_mc

# =====================================================================
# Physical Parameters (Farouk et al.)
# =====================================================================

# Rabi frequencies [MHz, angular]
omega_c_amp = 2 * np.pi * 50  # Control pulse amplitude
omega_p_amp = 2 * np.pi * 50  # Target probe pulse amplitude
omega_R_amp = 2.5 * omega_p_amp  # Target Rydberg coupling amplitude

# Detuning [MHz, angular]
delta = 2 * np.pi * 1200

# Rydberg blockade interaction strength [MHz, angular]
# *** MODIFY THIS VALUE TO EXPLORE PARAMETER SPACE ***
V_ct = 2 * np.pi * 200

# Number of Monte Carlo trajectories
# *** Use ntraj >= 500 for publication; 50 for quick testing ***
NTRAJ = 500

# Timing
T_c = np.pi / omega_c_amp  # Duration of control pi-pulse
T_f = 0.5  # Duration of target pulse window [us]
sigma = 0.05  # Target pulse width parameter [us^3]

# Decay rates [MHz] (gamma = 1 / lifetime)
gamma_r = 1.0 / 548.0  # Cs Rydberg |r>, lifetime 548 us
gamma_R = 1.0 / 505.0  # Rb Rydberg |R>, lifetime 505 us
gamma_P = 1.0 / 0.131  # Rb intermediate |P>, lifetime 0.131 us

# =====================================================================
# Build args dict for pulse functions
# =====================================================================
args = {
    "omega_c_amp": omega_c_amp,
    "omega_p_amp": omega_p_amp,
    "omega_R_amp": omega_R_amp,
    "T_c": T_c,
    "T_f": T_f,
    "sigma": sigma,
}

# =====================================================================
# Verify pulse area condition
# =====================================================================
area = compute_pulse_area(omega_p, delta, T_c, T_c + T_f, args)
print(f"Effective two-photon pulse area: {area:.4f} (target: pi = {np.pi:.4f})")

# =====================================================================
# Build Hamiltonian and collapse operators
# =====================================================================
H = build_hamiltonian(delta, V_ct)
c_ops = build_collapse_operators(gamma_r, gamma_R, gamma_P)

# =====================================================================
# Initial state: |1, 1, A>
# =====================================================================
cs = CsAtom()
rb = RbAtom()
psi0 = composite_basis_state(
    cs.level_index["1"],
    cs.level_index["1"],
    rb.level_index["A"],
)

psi_target = psi0  # |1, 1, A> -> |1, 1, A> (blockade active)

# =====================================================================
# Time list and expectation operators
# =====================================================================
t_total = 2 * T_c + T_f
tlist = np.linspace(0, t_total, 500)

e_ops = [
    composite_projector(2, rb.level_index["A"]),
    composite_projector(2, rb.level_index["B"]),
    composite_projector(2, rb.level_index["P"]),
    composite_projector(2, rb.level_index["R"]),
]

# =====================================================================
# Run simulation
# =====================================================================
print(
    f"Running mcsolve (ntraj={NTRAJ}, t_total={t_total:.4f} us, {len(tlist)} points)..."
)
result = simulate(
    method="mcsolve",
    H=H,
    psi0=psi0,
    tlist=tlist,
    c_ops=c_ops,
    e_ops=e_ops,
    options={"store_final_state": True, "nsteps": 10000},
    args=args,
    ntraj=NTRAJ,
)
print(f"Done. {result.num_trajectories} trajectories completed.")

# =====================================================================
# Output 1: Pulse envelope plot
# =====================================================================
fig1, ax1 = plt.subplots(figsize=(10, 4))
plot_pulses(
    tlist,
    [omega_c, omega_p, omega_R],
    [r"$\Omega_c$", r"$\Omega_p$", r"$\Omega_R$"],
    args,
    ax=ax1,
)
ax1.set_title("OR Gate Pulse Sequence")
fig1.tight_layout()
fig1.savefig("or_gate_mcsolve_pulses.png", dpi=150)
print("Saved: or_gate_mcsolve_pulses.png")

# =====================================================================
# Output 2: Averaged population dynamics with ±1σ bands
# =====================================================================
fig2, ax2 = plt.subplots(figsize=(10, 5))
plot_populations_mc(result, ["|A>", "|B>", "|P>", "|R>"], ax=ax2)
ax2.set_title(f"Population Dynamics (mcsolve, ntraj={NTRAJ}) — Initial: |1, 1, A>")
fig2.tight_layout()
fig2.savefig("or_gate_mcsolve_populations.png", dpi=150)
print("Saved: or_gate_mcsolve_populations.png")

# =====================================================================
# Output 3: Quantum jump statistics
# =====================================================================
print("\n--- Quantum Jump Statistics ---")
col_times = result.col_times
col_which = result.col_which

jump_counts = [len(ct) for ct in col_times]
print(f"Total trajectories: {len(jump_counts)}")
print(
    f"Jumps per trajectory: min={min(jump_counts)}, max={max(jump_counts)}, "
    f"mean={np.mean(jump_counts):.1f}"
)

# Count which collapse operator fired most often
all_which = []
for cw in col_which:
    all_which.extend(cw)
if all_which:
    channel_labels = [
        "c1: |r>->|0>",
        "c1: |r>->|1>",
        "c2: |r>->|0>",
        "c2: |r>->|1>",
        "t: |R>->|A>",
        "t: |R>->|B>",
        "t: |P>->|A>",
        "t: |P>->|B>",
    ]
    unique, counts = np.unique(all_which, return_counts=True)
    print("\nCollapse channel breakdown:")
    for ch_idx, count in zip(unique, counts):
        label = (
            channel_labels[ch_idx] if ch_idx < len(channel_labels) else f"ch {ch_idx}"
        )
        print(f"  {label}: {count} jumps ({100 * count / sum(counts):.1f}%)")
else:
    print("No quantum jumps occurred.")

# =====================================================================
# Output 4: State fidelity
# =====================================================================
fid = state_fidelity(result.final_state, psi_target)
print(f"\nState fidelity F(averaged_final, target) = {fid:.6f}")
print(f"Target state: |1, 1, A> (blockade should prevent transition)")

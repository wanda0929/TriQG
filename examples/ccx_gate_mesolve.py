"""
CCX (Toffoli) Gate Simulation using mesolve (Lindblad Master Equation)
======================================================================

Simulates the three-qubit Rydberg CCX gate with two 133-Cs control atoms
and one 87-Rb target atom. Uses QuTiP's mesolve via the triqg library.

The CCX gate flips the target qubit only when both control qubits are
in |1>. The protocol:
  1. Excite controls in |0> to Rydberg |r> (positive pi-pulse)
  2. Three-sub-pulse NOT sequence on the target:
     - Sub-pulse 1: |B> <-> |R>
     - Sub-pulse 2: |A> <-> |R>
     - Sub-pulse 3: |B> <-> |R>
  3. De-excite controls (negative pi-pulse)

When both controls are in |1>, no Rydberg excitation occurs, so the target
undergoes the full NOT. When any control is in |0>, it is excited to |r>,
blocking the target via the V_ct interaction.

Outputs:
  1. Pulse envelope plot
  2. Population dynamics plot
  3. Per-state fidelity against the expected ideal output
"""

import numpy as np
import matplotlib

matplotlib.use("Agg")  # Use non-interactive backend; remove for interactive use
import matplotlib.pyplot as plt

from triqg.atoms import CsAtom, RbAtom, composite_basis_state, composite_projector
from triqg.pulses import omega_cc, omega_t1, omega_t2
from triqg.hamiltonian import build_ccx_hamiltonian
from triqg.decoherence import build_collapse_operators
from triqg.solver import simulate
from triqg.analysis import state_fidelity
from triqg.visualization import plot_pulses, plot_populations

# =====================================================================
# Physical Parameters
# =====================================================================

# Rabi frequencies [MHz, angular]
omega_cc_amp = 2 * np.pi * 100  # Control pulse amplitude
omega_t_amp = 2 * np.pi * 50  # Target pulse amplitude

# Timing (derived from amplitudes)
T_cc = np.pi / omega_cc_amp  # Duration of control pi-pulse
T_t = np.pi / omega_t_amp  # Duration of each target sub-pulse

# Rydberg blockade interaction strength [MHz, angular]
V_ct = 2 * np.pi * 200

# Decay rates [MHz] (gamma = 1 / lifetime)
gamma_r = 1.0 / 548.0  # Cs Rydberg |r>, lifetime 548 us
gamma_R = 1.0 / 505.0  # Rb Rydberg |R>, lifetime 505 us
gamma_P = 1.0 / 0.131  # Rb intermediate |P>, lifetime 0.131 us

# =====================================================================
# Build args dict for pulse functions
# =====================================================================
args = {
    "omega_cc_amp": omega_cc_amp,
    "omega_t_amp": omega_t_amp,
    "T_cc": T_cc,
    "T_t": T_t,
}

print(f"T_cc = {T_cc:.6f} us ({T_cc * 1e3:.3f} ns)")
print(f"T_t  = {T_t:.6f} us ({T_t * 1e3:.3f} ns)")
print(
    f"Total gate time = {2 * T_cc + 3 * T_t:.6f} us ({(2 * T_cc + 3 * T_t) * 1e3:.3f} ns)"
)

# =====================================================================
# Build Hamiltonian and collapse operators
# =====================================================================
H = build_ccx_hamiltonian(V_ct)
c_ops = build_collapse_operators(gamma_r, gamma_R, gamma_P)

# =====================================================================
# Initial state: |1, 1, A> (both controls in |1>, target in |A>)
# =====================================================================
cs = CsAtom()
rb = RbAtom()
psi0 = composite_basis_state(
    cs.level_index["1"],
    cs.level_index["1"],
    rb.level_index["A"],
)

# Expected ideal output for CCX gate with |1,1> controls:
# Both controls in |1> -> no Rydberg excitation -> no blockade
# -> target flips: |A> -> |B>
psi_target = composite_basis_state(
    cs.level_index["1"],
    cs.level_index["1"],
    rb.level_index["B"],
)

# =====================================================================
# Time list and expectation operators
# =====================================================================
t_total = 2 * T_cc + 3 * T_t
tlist = np.linspace(0, t_total, 500)

# Population projectors for the target atom levels
e_ops = [
    composite_projector(2, rb.level_index["A"]),  # P(|A>)
    composite_projector(2, rb.level_index["B"]),  # P(|B>)
    composite_projector(2, rb.level_index["P"]),  # P(|P>)
    composite_projector(2, rb.level_index["R"]),  # P(|R>)
]

# =====================================================================
# Run simulation
# =====================================================================
# max_step is required so the adaptive ODE solver resolves the piecewise
# pulse windows (QuTiP 5 can otherwise step over the short active regions).
max_freq = max(omega_cc_amp, omega_t_amp) / (2 * np.pi)
max_step = 1.0 / (20 * max_freq)

print(f"\nRunning mesolve (t_total = {t_total:.6f} us, {len(tlist)} points)...")
result = simulate(
    method="mesolve",
    H=H,
    psi0=psi0,
    tlist=tlist,
    c_ops=c_ops,
    e_ops=e_ops,
    options={"store_final_state": True, "nsteps": 100000, "max_step": max_step},
    args=args,
)
print("Done.")

# =====================================================================
# Output 1: Pulse envelope plot
# =====================================================================
fig1, ax1 = plt.subplots(figsize=(10, 4))
plot_pulses(
    tlist,
    [omega_cc, omega_t1, omega_t2],
    [r"$\Omega_{cc}$", r"$\Omega_{t1}$", r"$\Omega_{t2}$"],
    args,
    ax=ax1,
)
ax1.set_title("CCX Gate Pulse Sequence")
fig1.tight_layout()
fig1.savefig("ccx_gate_mesolve_pulses.png", dpi=150)
print("\nSaved: ccx_gate_mesolve_pulses.png")

# =====================================================================
# Output 2: Population dynamics plot
# =====================================================================
fig2, ax2 = plt.subplots(figsize=(10, 5))
plot_populations(result, ["|A>", "|B>", "|P>", "|R>"], ax=ax2)
ax2.set_title(r"Population Dynamics (mesolve) — Initial: |1, 1, A$\rangle$")
fig2.tight_layout()
fig2.savefig("ccx_gate_mesolve_populations.png", dpi=150)
print("Saved: ccx_gate_mesolve_populations.png")

# =====================================================================
# Output 3: State fidelity
# =====================================================================
fid = state_fidelity(result.final_state, psi_target)
print(f"\nState fidelity F(final, target) = {fid:.6f}")
print(f"Target state: |1, 1, B> (both controls |1> -> target should flip)")

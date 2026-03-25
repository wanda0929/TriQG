"""
OR Gate Simulation using mesolve with Gaussian Pulse
=====================================================

Simulates the three-qubit Rydberg OR gate with two 133-Cs control atoms
and one 87-Rb target atom. Uses QuTiP's mesolve via the triqg library.

This version uses a super-Gaussian pulse (order 6) instead of the Hanning pulse.

Reference: M. Farouk et al., "Parallel Implementation of CNOTN and C2NOT2
Gates via Homonuclear and Heteronuclear Förster Interactions of Rydberg Atoms."

Outputs:
  1. Pulse envelope plot
  2. Population dynamics plot
  3. State fidelity against the expected ideal output
"""

import numpy as np
import matplotlib

matplotlib.use("Agg")  # Use non-interactive backend; remove for interactive use
import matplotlib.pyplot as plt

from triqg.atoms import CsAtom, RbAtom, composite_basis_state, composite_projector
from triqg.pulses import omega_c, omega_gaussian, omega_R, compute_pulse_area
from triqg.hamiltonian import build_hamiltonian
from triqg.decoherence import build_collapse_operators
from triqg.solver import simulate
from triqg.analysis import state_fidelity
from triqg.visualization import plot_pulses, plot_populations

# =====================================================================
# Physical Parameters (Farouk et al.)
# =====================================================================

# Rabi frequencies [MHz, angular]
omega_c_amp = 2 * np.pi * 50  # Control pulse amplitude
omega_p_amp = 2 * np.pi * 50.0 * 1.039975   # Target probe pulse amplitude
omega_R_amp = 3.5 * omega_p_amp  # Target Rydberg coupling amplitude

# Detuning [MHz, angular]
delta = 2 * np.pi * 500

# Rydberg blockade interaction strength [MHz, angular]
V_ct = 2 * np.pi * 5000

# Timing
T_c = np.pi / omega_c_amp  # Duration of control pi-pulse
T_f = 0.15  # Duration of target pulse half-window [us]
# Super-Gaussian pulse: shape = exp(-((t - t_center)^3 / sigma)^2)
# The area is determined by numerical integration

# Gaussian pulse width parameter
sigma = 0.0014  # Width parameter for super-Gaussian pulse

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
area = compute_pulse_area(omega_gaussian, delta, T_c, T_c + 2 * T_f, args)
print(f"Effective two-photon pulse area: {area:.4f} (target: pi = {np.pi:.4f})")

# =====================================================================
# Build Hamiltonian and collapse operators
# =====================================================================
H = build_hamiltonian(delta, V_ct, pulse_p=omega_gaussian)
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


# Expected ideal output for OR gate with |1,1> controls:
# When both controls are in |1>, they get excited to |r>,
# the Rydberg blockade prevents the target transition,
# so the target should remain in |A>.
psi_target = composite_basis_state(
    cs.level_index["1"],
    cs.level_index["1"],
    rb.level_index["B"],
)  # |1, 1, A> -> |1, 1, B> (blockade active)

# =====================================================================
# Time list and expectation operators
# =====================================================================
t_total = 2 * T_c + 2 * T_f
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
print(f"Running mesolve (t_total = {t_total:.4f} us, {len(tlist)} points)...")
result = simulate(
    method="mesolve",
    H=H,
    psi0=psi0,
    tlist=tlist,
    c_ops=c_ops,
    e_ops=e_ops,
    options={"store_final_state": True, "nsteps": 10000},
    args=args,
)
print("Done.")


# =====================================================================
# Output 1: Pulse envelope plot
# =====================================================================
fig1, ax1 = plt.subplots(figsize=(10, 4))
plot_pulses(
    tlist,
    [omega_c, omega_gaussian, omega_R],
    [r"$\Omega_c$", r"$\Omega_p$ (Gaussian)", r"$\Omega_R$"],
    args,
    ax=ax1,
)
ax1.set_title("OR Gate Pulse Sequence (Gaussian Pulse)")
fig1.tight_layout()
fig1.savefig("or_gate_mesolve_gaussian_pulses.png", dpi=150)
print("Saved: or_gate_mesolve_gaussian_pulses.png")


# =====================================================================
# Output 2: Population dynamics plot
# =====================================================================
fig2, ax2 = plt.subplots(figsize=(10, 5))
plot_populations(result, ["|A>", "|B>", "|P>", "|R>"], ax=ax2)
ax2.set_title("Population Dynamics (mesolve, Gaussian) — Initial: |1, 1, A>")
fig2.tight_layout()
fig2.savefig("or_gate_mesolve_gaussian_populations.png", dpi=150)
print("Saved: or_gate_mesolve_gaussian_populations.png")

# =====================================================================
# Output 3: State fidelity
# =====================================================================
fid = state_fidelity(result.final_state, psi_target)
print(f"\nState fidelity F(final, target) = {fid:.6f}")
print(f"Target state: |1, 1, B> (blockade should prevent transition)")

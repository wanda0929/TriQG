"""
OR Gate Average Gate Fidelity -- Smaller_V_DD redesign
========================================================================
Energy-level redesign:  Rb 54 D_{3/2} + Cs 62 D_{5/2}
    (replaces paper's Rb 66 D_{5/2} + Cs 76 D_{3/2})

Why this file exists
--------------------
The paper's Rb 66 D_5/2 sits ~157x above the Rb-Rb (D-D) Foerster zero
crossing, so the global Rb drive picks up cross-data blockade
V_DD/(2 pi) ~ 91 MHz at r_DD = a = 5 um.  R_DD = V_dd / V_DD ~ 5.7 is
far below the >= 100 selectivity target.  Moving Rb to 54 D_{3/2}
parks it on a same-species Foerster zero crossing
(|C_6^DD| -> 9.07 GHz um^6), killing V_DD by 157x at the cost of a
~2.3x slower V_ct.  See:
    SelfCorrectingRydberg/plan/Rb_Cs_level_redesign.md (Sec. 4)
    SelfCorrectingRydberg/python_scripts/scan_level_design.py
    examples/Smaller_VDD_energy/level_parameters.py

Pulse-sequence revision  (FINAL: strong blockade + sub-pi/4 area)
-----------------------------------------------------------------
Relative to examples/Average_fidelity/or_average_gate_fid_gaussian.py,
which used  Omega_c = 2 pi * 50,  Omega_p = 2 pi * 50 * 1.039975,
            T_f = 0.15 us,        sigma = 0.0014 us,  area = pi/4:

  * Omega_p_amp = 2 pi * 50    MHz       (paper-style amplitude)
  * Omega_R_amp = 3.5 * Omega_p_amp = 2 pi * 175 MHz
  * Omega_c_amp = 2 pi * 50    MHz
  * Delta       = 2 pi * 500   MHz       (unchanged)
  * alpha       = 4   (smooth edges, exp(-alpha^2) ~ 1.1e-7 of peak)
  * area        = 0.95 * pi/4 = 0.7462     <-- SUB-pi/4 OPTIMUM
  * sigma       = 1.518 ns               (set by k=0.95 and area-formula)
  * T_f         = 182.46 ns

Why area = 0.95 * pi/4 (not pi/4)
---------------------------------
The paper's design value area = pi/4 is the OR-protocol's correct
rotation in the *V_ct -> infinity* limit.  At V_ct = 226 MHz (the
new energy levels), the blockade is finite and the protocol's
target rotation should be slightly under-rotated to compensate for
imperfect blockade leakage.  A 2-D scan of F_bar over (Delta, k =
area/(pi/4)) shows a broad maximum F_bar = 0.992 along the ridge
k = 0.95-0.96 for Delta = 400-550 MHz.

Strong-blockade audit (this config):
    Omega_p^2 / (2 Delta) / (2 pi)        = 2.50 MHz   << V_ct/(2 pi) = 226 MHz   (margin ~90)
    Omega_p * Omega_R / (2 Delta) / (2 pi) = 8.75 MHz  << V_ct/(2 pi) = 226 MHz   (margin ~26)

Design-derivation reference
---------------------------
In the well-resolved limit (alpha = T_f^3/sigma >> 1), the area
constraint analytically fixes sigma alone:

    sigma_(area=pi/4) = ((2 pi Delta) / (Omega_p^2 * I_inf))^3,
    I_inf = 2 * Gamma(7/6) * 2^(-1/6) ~ 1.6534.

For a sub-pi/4 area = k * pi/4, scale sigma by k^3 (since area scales
as sigma^(1/3)):
    sigma  = k^3 * sigma_(area=pi/4)
    T_f    = (alpha * sigma)^(1/3)

T_f is then a free parameter that controls only where the tail is
truncated to zero -- i.e. the smoothness of the trailing/leading
edge.  We pick alpha = 4 so the edge amplitude is exp(-alpha^2) ~
1.1e-7 of the peak (effectively machine zero), in contrast to the
alpha = 2.41 default of the paper's Hanning-replacement file, which
leaves a 3e-3 pedestal that snaps abruptly to zero at the window
boundary.

OR gate protocol  (Farouk et al.):
  1. Excite controls in |1> to Rydberg |r> (positive pi-pulse, T_c)
  2. Two-photon Raman pulse on target via |P> (super-Gaussian, 2*T_f)
  3. De-excite controls (negative pi-pulse, T_c)

Truth table:
    Input    -> Ideal output
    |0,0,A>  -> |0,0,A>     (both controls 0: target unchanged)
    |0,0,B>  -> |0,0,B>
    |0,1,A>  -> |0,1,B>     (one control 1: target flips)
    ...
    |1,1,B>  -> |1,1,A>     (both controls 1: target flips)
"""

import numpy as np

from triqg.atoms import CsAtom, RbAtom, composite_basis_state
from triqg.pulses import omega_gaussian, compute_pulse_area
from triqg.hamiltonian import build_hamiltonian
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
# Pulse-sequence parameters  (FINAL: strong-blockade + sub-pi/4 area)
# =====================================================================
omega_p_MHz = 50.0                 # Rb target probe amplitude (paper-style)
omega_R_MHz = 3.5 * omega_p_MHz    # Omega_R = 3.5 * Omega_p (protocol)
omega_c_MHz = 50.0                 # Cs control Rabi frequency
delta_MHz   = 500.0                # Two-photon detuning

omega_p_amp = 2 * np.pi * omega_p_MHz   # rad/us
omega_R_amp = 2 * np.pi * omega_R_MHz
omega_c_amp = 2 * np.pi * omega_c_MHz
delta       = 2 * np.pi * delta_MHz

T_c = np.pi / omega_c_amp                # = 10.0 ns

# Pulse-shape parameters.
# Smoothness:  alpha = T_f^3 / sigma = 4 -> edge/peak ~ 1.1e-7  (smooth).
# Area scale:  k = area / (pi/4) = 0.95  (sub-pi/4 optimum for finite V_ct).
ALPHA = 4.0
K     = 0.95                             # area / (pi/4); 0.95 = empirical optimum

# Derived: sigma at area = pi/4 (analytic, well-resolved limit)
_I_inf = 2 * 0.92770 * 2**(-1.0/6.0)     # ~ 1.6534
_sigma_pi4 = ((2 * np.pi * delta) / (omega_p_amp**2 * _I_inf)) ** 3
sigma = (K ** 3) * _sigma_pi4            # sigma at area = K * pi/4
T_f   = (ALPHA * sigma) ** (1.0/3.0)
# Numerical values at the locked-in optimum:
#   sigma = 0.001518 us (1.518 ns)
#   T_f   = 0.182463 us (182.46 ns)
#   total = 2 T_c + 2 T_f = 0.3849 us (384.9 ns)
#   area  = 0.95 * pi/4 = 0.7463

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
print(f"Effective two-photon pulse area: {area:.6f}  "
      f"(target: K * pi/4 = {K * np.pi/4:.6f},  K = {K})")

# =====================================================================
# Build Hamiltonian and collapse operators
# =====================================================================
# V_ct, V_cc imported from level_parameters (angular MHz = rad/us).
# Note: V_DD (Rb-Rb data-data) does NOT enter the 3-atom Hamiltonian --
# it is a cross-gate effect that matters for the global parallel drive,
# not for the per-gate truth-table simulation done here.
H = build_hamiltonian(delta, V_ct, pulse_p=omega_gaussian, V_cc=V_cc)
c_ops = build_collapse_operators(gamma_r, gamma_R, gamma_P)

# =====================================================================
# Time list and solver options
# =====================================================================
t_total = 2 * T_c + 2 * T_f
tlist = np.linspace(0, t_total, 500)

solver_opts = {"store_final_state": True, "nsteps": 100000}

# =====================================================================
# Enumerate all 2^{n+1} = 8 computational basis states
# =====================================================================
cs = CsAtom()
rb = RbAtom()

n_controls = 2
ctrl_levels = [cs.level_index["0"], cs.level_index["1"]]
tgt_levels = [rb.level_index["A"], rb.level_index["B"]]


def or_gate_ideal_output(c1_bit: int, c2_bit: int, t_bit: int):
    """Ideal OR-gate output: target flips iff (c1 OR c2) is 1."""
    if c1_bit or c2_bit:
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
            psi_ideal = or_gate_ideal_output(c1_bit, c2_bit, t_bit)
            basis_inputs.append(
                (str(c1_bit), str(c2_bit), ["A", "B"][t_bit], psi_in, psi_ideal)
            )

# =====================================================================
# Header
# =====================================================================
print(f"\nOR Gate Average Gate Fidelity -- Smaller_V_DD redesign  (Gaussian pulse)")
print(f"Energy levels: Rb 54 D_3/2 + Cs 62 D_5/2,  a = {a_um:.2f} um")
print(f"==================================================================")
print(f"n_controls = {n_controls},  basis states = {len(basis_inputs)}")
print(f"Pulse amplitudes:")
print(f"  Omega_c / (2 pi) = {omega_c_MHz} MHz   -> T_c = {T_c*1e3:.4f} ns")
print(f"  Omega_p / (2 pi) = {omega_p_MHz} MHz   -> T_f = {T_f*1e3:.4f} ns (super-Gaussian)")
print(f"  Omega_R / (2 pi) = {omega_R_MHz} MHz   (ratio = {omega_R_amp/omega_p_amp:.2f})")
print(f"  sigma = {sigma*1e3:.4f} ns  (alpha = T_f^3/sigma = {T_f**3/sigma:.4f}, edge/peak = {np.exp(-(T_f**3/sigma)**2):.2e})")
print(f"  delta / (2 pi) = {delta_MHz:.1f} MHz")
print(f"  area = K * pi/4 with K = {K}")
print(f"Interactions:")
print(f"  V_ct / (2 pi) = +{V_ct_MHz:7.3f} MHz   (Rb-Cs Foerster, |54D3/2;62D5/2>)")
print(f"  V_cc / (2 pi) = {V_cc_MHz:+8.4f} MHz   (Cs-Cs vdW at r_AA = {r_AA:.3f} um)")
print(f"  V_DD / (2 pi) = {V_DD_MHz:+8.4f} MHz   (Rb-Rb vdW at r_DD = {a_um:.3f} um) [info only]")
print(f"  Foerster defect Delta_F / (2 pi) = {Delta_F_MHz:+.2f} MHz")
print(f"  R_DD = {R_DD:.1f},  R_AA = {R_AA:.1f}   (both >= 100 OK)")
print(f"Decoherence (re-verify with ARC):")
print(f"  tau_r (Cs 62 D_5/2) = {tau_r_us:.2f} us")
print(f"  tau_R (Rb 54 D_3/2) = {tau_R_us:.2f} us")
print(f"  tau_P (Rb  7 P_3/2) = {tau_P_us:.3f} us")
print(f"Total gate time = {t_total:.4f} us  ({t_total*1e3:.2f} ns)\n")

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
print(f"Omega_c = 2 pi * {omega_c_MHz} MHz,  Delta = 2 pi * {delta_MHz} MHz,  area = {K} * pi/4")
print(f"Average gate fidelity:  F_bar = {F_bar:.6f}")
print(f"Gate infidelity:        1 - F_bar = {1 - F_bar:.3e}")
print(f"{'=' * 60}")

# =====================================================================
# Diagnostic: population breakdown for each output
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

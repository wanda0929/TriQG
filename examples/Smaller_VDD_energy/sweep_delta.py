"""
Strong-blockade sweep over Delta for the OR-gate Raman drive.

Holds fixed:
  Omega_p / (2 pi) = 50 MHz                  (paper-style amplitude)
  Omega_R / (2 pi) = 3.5 * Omega_p           (175 MHz)
  Omega_c / (2 pi) = 50 MHz
  alpha            = 4   (smooth edge, exp(-16) ~ 1e-7)
  Energy levels    = Rb 54 D_3/2 + Cs 62 D_5/2  -> V_ct/(2 pi) = 226.27 MHz

For each Delta, the area = pi/4 constraint sets sigma analytically,
and T_f follows from alpha:

    sigma^(1/3) = (2 pi Delta) / (Omega_p^2 * I_inf),  I_inf ~ 1.6534
    T_f         = (alpha * sigma)^(1/3)

Reports:
  - Two blockade margins:
      M1 = V_ct / [Omega_p^2 / (2 Delta)]          (AC Stark on |P>)
      M2 = V_ct / [Omega_p * Omega_R / (2 Delta)]  (two-photon Raman)
  - Average OR-gate fidelity F_bar from mesolve.
"""
import numpy as np
from triqg.atoms import CsAtom, RbAtom, composite_basis_state
from triqg.pulses import omega_gaussian, compute_pulse_area
from triqg.hamiltonian import build_hamiltonian
from triqg.decoherence import build_collapse_operators
from triqg.solver import simulate
from triqg.analysis import state_fidelity, average_gate_fidelity

from level_parameters import V_ct, V_cc, V_ct_MHz, gamma_r, gamma_R, gamma_P

# Fixed parameters
omega_p_MHz = 50.0
omega_R_MHz = 3.5 * omega_p_MHz
omega_c_MHz = 50.0
alpha       = 4.0

omega_p_amp = 2 * np.pi * omega_p_MHz
omega_R_amp = 2 * np.pi * omega_R_MHz
omega_c_amp = 2 * np.pi * omega_c_MHz
T_c = np.pi / omega_c_amp

I_inf = 2 * 0.92770 * 2**(-1/6)   # ~ 1.6534, integral of exp(-2 v^6) over R

def sigma_for_area(delta_MHz):
    delta = 2 * np.pi * delta_MHz
    s_third = (2 * np.pi * delta) / (omega_p_amp**2 * I_inf)
    return s_third**3   # in us

# Delta sweep (MHz)
deltas = [300.0, 400.0, 500.0, 600.0, 800.0, 1000.0, 1500.0]

# Solver setup
cs = CsAtom(); rb = RbAtom()
ctrl_levels = [cs.level_index["0"], cs.level_index["1"]]
tgt_levels  = [rb.level_index["A"], rb.level_index["B"]]

def or_ideal(c1, c2, t):
    return composite_basis_state(
        ctrl_levels[c1], ctrl_levels[c2],
        tgt_levels[1 - t if (c1 or c2) else t],
    )

basis_inputs = []
for c1 in range(2):
    for c2 in range(2):
        for t in range(2):
            psi_in = composite_basis_state(ctrl_levels[c1], ctrl_levels[c2], tgt_levels[t])
            basis_inputs.append((c1, c2, t, psi_in, or_ideal(c1, c2, t)))

c_ops = build_collapse_operators(gamma_r, gamma_R, gamma_P)

print(f"{'delta':>7s} {'M1':>6s} {'M2':>6s} {'sigma(ns)':>9s} {'T_f(ns)':>8s} "
      f"{'2Tf(ns)':>8s} {'area':>8s} {'F_bar':>9s} {'1-Fbar':>10s}")
print("-" * 88)

results = []
for delta_MHz in deltas:
    delta = 2 * np.pi * delta_MHz
    sigma = sigma_for_area(delta_MHz)
    T_f   = (alpha * sigma) ** (1.0/3.0)

    # Blockade margins (dimensionless)
    one_photon_AC  = omega_p_MHz**2 / (2 * delta_MHz)       # MHz
    two_photon_Rab = omega_p_MHz * omega_R_MHz / (2 * delta_MHz)  # MHz
    M1 = V_ct_MHz / one_photon_AC
    M2 = V_ct_MHz / two_photon_Rab

    args = {"omega_c_amp": omega_c_amp, "omega_p_amp": omega_p_amp,
            "omega_R_amp": omega_R_amp, "T_c": T_c, "T_f": T_f, "sigma": sigma}

    # Verify area numerically
    area = compute_pulse_area(omega_gaussian, delta, T_c, T_c + 2*T_f, args)

    # Build H, run all 8 inputs
    H = build_hamiltonian(delta, V_ct, pulse_p=omega_gaussian, V_cc=V_cc)
    t_total = 2 * T_c + 2 * T_f
    tlist = np.linspace(0, t_total, 400)
    solver_opts = {"store_final_state": True, "nsteps": 200000}

    pairs = []
    for c1, c2, t, psi_in, psi_ideal in basis_inputs:
        res = simulate(method="mesolve", H=H, psi0=psi_in, tlist=tlist,
                       c_ops=c_ops, e_ops=[], options=solver_opts, args=args)
        pairs.append((res.final_state, psi_ideal))

    F_bar = average_gate_fidelity(pairs)
    per_input = [state_fidelity(rho, psi) for (rho, psi) in pairs]

    print(f"{delta_MHz:7.1f} {M1:6.1f} {M2:6.1f} {sigma*1e3:9.3f} "
          f"{T_f*1e3:8.2f} {2*T_f*1e3:8.2f} {area:8.5f} {F_bar:9.6f} {1-F_bar:10.3e}")
    results.append((delta_MHz, M1, M2, sigma, T_f, area, F_bar, per_input))

print()
print("Per-input fidelities (rows: input; cols: delta values)")
print(f"{'input':<10s}" + "".join(f"{d:>10.0f} " for d in deltas))
labels = [f"|{c1}{c2}{['A','B'][t]}>" for c1,c2,t,_,_ in basis_inputs]
for i, lbl in enumerate(labels):
    row = "".join(f"{r[7][i]:>10.6f} " for r in results)
    print(f"{lbl:<10s}{row}")

best = max(results, key=lambda r: r[6])
print(f"\nBest F_bar = {best[6]:.6f} at Delta = {best[0]:.0f} MHz (M1={best[1]:.1f}, M2={best[2]:.1f})")

"""
Sweep Omega_p (paired with Delta at fixed blockade margin) for the
OR-gate Raman drive.

Constraint:  M2 = V_ct / (Omega_p Omega_R / (2 Delta)) = const = 30
             with Omega_R = 3.5 Omega_p  =>  Omega_p^2 / Delta = 2 V_ct / (3.5 * 30) = const

So Delta scales as Omega_p^2 at fixed margin.  T_f then depends only on V_ct
and M, not on Omega_p separately.  What DOES change with Omega_p:
    V_ct / Omega_R  (EIT cleanness; we want this >> 1)

Fixed:
    Omega_c / (2 pi) = 50 MHz
    Omega_R / Omega_p = 3.5
    alpha = 4 (smooth)
    Energy levels: Rb 54 D_3/2 + Cs 62 D_5/2  -> V_ct/(2 pi) = 226.27 MHz
    Blockade margin M2 = 30
"""
import numpy as np
from triqg.atoms import CsAtom, RbAtom, composite_basis_state
from triqg.pulses import omega_gaussian, compute_pulse_area
from triqg.hamiltonian import build_hamiltonian
from triqg.decoherence import build_collapse_operators
from triqg.solver import simulate
from triqg.analysis import state_fidelity, average_gate_fidelity

from level_parameters import V_ct, V_cc, V_ct_MHz, gamma_r, gamma_R, gamma_P

omega_c_MHz = 50.0
alpha       = 4.0
M_target    = 30.0    # blockade margin V_ct / (Omega_p Omega_R / (2 Delta))

omega_c_amp = 2 * np.pi * omega_c_MHz
T_c = np.pi / omega_c_amp

I_inf = 2 * 0.92770 * 2**(-1/6)

cs = CsAtom(); rb = RbAtom()
ctrl_levels = [cs.level_index["0"], cs.level_index["1"]]
tgt_levels  = [rb.level_index["A"], rb.level_index["B"]]

basis_inputs = []
for c1 in range(2):
    for c2 in range(2):
        for t in range(2):
            psi_in = composite_basis_state(ctrl_levels[c1], ctrl_levels[c2], tgt_levels[t])
            psi_id = composite_basis_state(
                ctrl_levels[c1], ctrl_levels[c2],
                tgt_levels[1 - t if (c1 or c2) else t]
            )
            basis_inputs.append((c1,c2,t, psi_in, psi_id))

c_ops = build_collapse_operators(gamma_r, gamma_R, gamma_P)

omega_p_values = [10.0, 15.0, 20.0, 25.0, 30.0, 40.0, 50.0]
print(f"V_ct/(2 pi) = {V_ct_MHz:.2f} MHz,  fixed M2 = {M_target}, alpha = {alpha}\n")
print(f"{'Op':>4s} {'OR':>5s} {'Delta':>6s} {'Vct/OR':>7s} {'sigma':>8s} "
      f"{'T_f':>7s} {'2Tf+Tc':>7s} {'M1':>5s} {'M2':>5s} {'F_bar':>9s} {'1-F':>9s}")
print("-" * 90)

results = []
for omega_p_MHz in omega_p_values:
    omega_R_MHz = 3.5 * omega_p_MHz
    omega_p_amp = 2 * np.pi * omega_p_MHz
    omega_R_amp = 2 * np.pi * omega_R_MHz

    # Delta from M_target = V_ct / (Omega_p Omega_R / (2 Delta))
    delta_MHz = (omega_p_MHz * omega_R_MHz) / (2 * V_ct_MHz / M_target)
    delta = 2 * np.pi * delta_MHz

    # Sigma from area = pi/4
    sigma_third = (2 * np.pi * delta) / (omega_p_amp**2 * I_inf)
    sigma = sigma_third**3
    T_f = (alpha * sigma) ** (1.0/3.0)

    one_photon_AC  = omega_p_MHz**2 / (2 * delta_MHz)
    two_photon_Rab = omega_p_MHz * omega_R_MHz / (2 * delta_MHz)
    M1 = V_ct_MHz / one_photon_AC
    M2 = V_ct_MHz / two_photon_Rab
    Vct_over_OR = V_ct_MHz / omega_R_MHz

    args = {"omega_c_amp": omega_c_amp, "omega_p_amp": omega_p_amp,
            "omega_R_amp": omega_R_amp, "T_c": T_c, "T_f": T_f, "sigma": sigma}
    area = compute_pulse_area(omega_gaussian, delta, T_c, T_c + 2*T_f, args)

    H = build_hamiltonian(delta, V_ct, pulse_p=omega_gaussian, V_cc=V_cc)
    t_total = 2*T_c + 2*T_f
    tlist = np.linspace(0, t_total, 400)
    opts = {"store_final_state": True, "nsteps": 200000}

    pairs = []
    for c1,c2,t,psi_in,psi_id in basis_inputs:
        res = simulate(method="mesolve", H=H, psi0=psi_in, tlist=tlist,
                       c_ops=c_ops, e_ops=[], options=opts, args=args)
        pairs.append((res.final_state, psi_id))
    F_bar = average_gate_fidelity(pairs)
    per_input = [state_fidelity(rho, psi) for (rho, psi) in pairs]

    print(f"{omega_p_MHz:4.0f} {omega_R_MHz:5.0f} {delta_MHz:6.0f} {Vct_over_OR:7.2f} "
          f"{sigma*1e3:8.3f} {T_f*1e3:7.1f} {(t_total)*1e3:7.1f} "
          f"{M1:5.1f} {M2:5.1f} {F_bar:9.6f} {1-F_bar:9.3e}")
    results.append((omega_p_MHz, omega_R_MHz, delta_MHz, sigma, T_f, F_bar, per_input))

print()
print("Per-input fidelities (rows: input; cols: Omega_p values)")
print(f"{'input':<8s}" + "".join(f"{op:>10.0f} " for op in omega_p_values))
labels = [f"|{c1}{c2}{['A','B'][t]}>" for c1,c2,t,_,_ in basis_inputs]
for i, lbl in enumerate(labels):
    print(f"{lbl:<8s}" + "".join(f"{r[6][i]:>10.6f} " for r in results))

best = max(results, key=lambda r: r[5])
print(f"\nBest F_bar = {best[5]:.6f} at Omega_p = {best[0]:.0f} MHz, "
      f"Delta = {best[2]:.0f} MHz, V_ct/Omega_R = {V_ct_MHz/best[1]:.2f}")

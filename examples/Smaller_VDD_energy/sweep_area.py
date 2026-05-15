"""
Sweep the two-photon pulse area to find the optimum at finite V_ct.

The paper's protocol assumes area = pi/4 in the V_ct -> infinity limit.
At V_ct = 226 MHz, the blockade is finite, and a different area may
give a cleaner target rotation in the blockaded branches.

Scan:  area_factor in [0.6, 1.6] x (pi/4),  by scaling sigma at fixed
       Omega_p, Delta, alpha.
       Then T_f = (alpha * sigma)^(1/3) follows.
"""
import numpy as np
from triqg.atoms import CsAtom, RbAtom, composite_basis_state
from triqg.pulses import omega_gaussian, compute_pulse_area
from triqg.hamiltonian import build_hamiltonian
from triqg.decoherence import build_collapse_operators
from triqg.solver import simulate
from triqg.analysis import state_fidelity, average_gate_fidelity

from level_parameters import V_ct, V_cc, V_ct_MHz, gamma_r, gamma_R, gamma_P

omega_p_MHz = 50.0
omega_R_MHz = 3.5 * omega_p_MHz
omega_c_MHz = 50.0
delta_MHz   = 700.0     # near the F_bar peak from sweep_delta.py
alpha       = 4.0

omega_p_amp = 2 * np.pi * omega_p_MHz
omega_R_amp = 2 * np.pi * omega_R_MHz
omega_c_amp = 2 * np.pi * omega_c_MHz
delta       = 2 * np.pi * delta_MHz
T_c         = np.pi / omega_c_amp

I_inf = 2 * 0.92770 * 2**(-1/6)
sigma_pi4_third = (2 * np.pi * delta) / (omega_p_amp**2 * I_inf)
sigma_pi4 = sigma_pi4_third**3   # value that gives area = pi/4 exactly

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
H = build_hamiltonian(delta, V_ct, pulse_p=omega_gaussian, V_cc=V_cc)

# Scan: sigma = k^3 * sigma_pi4  =>  area = k * pi/4 (in well-resolved limit)
ks = [0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15, 1.20, 1.30, 1.40]

print(f"Fixed: Op={omega_p_MHz}, OR={omega_R_MHz}, Oc={omega_c_MHz} MHz, "
      f"Delta={delta_MHz} MHz, alpha={alpha}")
print(f"V_ct/(2pi)={V_ct_MHz:.1f} MHz\n")
print(f"{'k':>6s} {'sigma':>8s} {'T_f':>7s} {'2Tf':>7s} {'area':>8s} {'area/pi4':>9s} "
      f"{'F_bar':>9s} {'F00':>8s} {'F01':>8s} {'F11':>8s}")
print("-" * 95)

results = []
for k in ks:
    sigma = k**3 * sigma_pi4
    T_f = (alpha * sigma) ** (1.0/3.0)

    args = {"omega_c_amp": omega_c_amp, "omega_p_amp": omega_p_amp,
            "omega_R_amp": omega_R_amp, "T_c": T_c, "T_f": T_f, "sigma": sigma}
    area = compute_pulse_area(omega_gaussian, delta, T_c, T_c + 2*T_f, args)

    tlist = np.linspace(0, 2*T_c + 2*T_f, 400)
    opts  = {"store_final_state": True, "nsteps": 200000}

    pairs = []
    fids = []
    for c1,c2,t,psi_in,psi_id in basis_inputs:
        res = simulate(method="mesolve", H=H, psi0=psi_in, tlist=tlist,
                       c_ops=c_ops, e_ops=[], options=opts, args=args)
        pairs.append((res.final_state, psi_id))
        fids.append(state_fidelity(res.final_state, psi_id))
    F_bar = average_gate_fidelity(pairs)

    print(f"{k:6.3f} {sigma*1e3:8.3f} {T_f*1e3:7.1f} {2*T_f*1e3:7.1f} "
          f"{area:8.5f} {area/(np.pi/4):9.5f} {F_bar:9.6f} "
          f"{fids[0]:8.5f} {fids[2]:8.5f} {fids[6]:8.5f}")
    results.append((k, sigma, T_f, area, F_bar, fids))

best = max(results, key=lambda r: r[4])
print(f"\nBest F_bar = {best[4]:.6f} at k = {best[0]:.3f}, "
      f"area = {best[3]:.5f} = {best[3]/(np.pi/4):.4f} pi/4")

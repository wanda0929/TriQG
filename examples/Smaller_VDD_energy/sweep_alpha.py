"""
Sweep alpha (= T_f^3 / sigma, the super-Gaussian shape factor) at
fixed Delta, to quantify the smoothness <-> decay tradeoff.
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
delta_MHz   = 500.0       # paper's original detuning

omega_p_amp = 2 * np.pi * omega_p_MHz
omega_R_amp = 2 * np.pi * omega_R_MHz
omega_c_amp = 2 * np.pi * omega_c_MHz
delta       = 2 * np.pi * delta_MHz
T_c         = np.pi / omega_c_amp

I_inf = 2 * 0.92770 * 2**(-1/6)
sigma_third = (2 * np.pi * delta) / (omega_p_amp**2 * I_inf)
sigma = sigma_third**3

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

alphas = [1.5, 2.0, 2.4107, 2.8, 3.0, 3.5, 4.0, 5.0, 6.0]
print(f"Fixed: Omega_p={omega_p_MHz} MHz, Omega_R={omega_R_MHz} MHz, "
      f"Omega_c={omega_c_MHz} MHz, Delta={delta_MHz} MHz")
print(f"       V_ct/(2pi)={V_ct_MHz:.1f} MHz, sigma={sigma*1e3:.3f} ns "
      f"(set by area=pi/4)\n")
print(f"{'alpha':>7s} {'edge/pk':>10s} {'T_f(ns)':>8s} {'2Tf(ns)':>8s} "
      f"{'area':>8s} {'F_bar':>9s} {'1-Fbar':>10s}")
print("-" * 72)

for alpha in alphas:
    T_f = (alpha * sigma) ** (1.0/3.0)
    args = {"omega_c_amp": omega_c_amp, "omega_p_amp": omega_p_amp,
            "omega_R_amp": omega_R_amp, "T_c": T_c, "T_f": T_f, "sigma": sigma}
    area = compute_pulse_area(omega_gaussian, delta, T_c, T_c + 2*T_f, args)
    edge = np.exp(-alpha**2)

    tlist = np.linspace(0, 2*T_c + 2*T_f, 400)
    opts  = {"store_final_state": True, "nsteps": 200000}

    pairs = []
    for c1,c2,t,psi_in,psi_id in basis_inputs:
        res = simulate(method="mesolve", H=H, psi0=psi_in, tlist=tlist,
                       c_ops=c_ops, e_ops=[], options=opts, args=args)
        pairs.append((res.final_state, psi_id))
    F_bar = average_gate_fidelity(pairs)

    print(f"{alpha:7.4f} {edge:10.2e} {T_f*1e3:8.2f} {2*T_f*1e3:8.2f} "
          f"{area:8.5f} {F_bar:9.6f} {1-F_bar:10.3e}")

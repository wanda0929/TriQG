"""
2-D scan over (Delta, area_factor k) to locate the joint optimum
for F_bar.  Same fixed parameters as sweep_omega_p.py + sweep_area.py:
  Omega_p = 2 pi * 50 MHz, Omega_R = 2 pi * 175 MHz, Omega_c = 2 pi * 50 MHz
  alpha = 4 (smooth edge)
  Energy levels: Rb 54 D_3/2 + Cs 62 D_5/2  -> V_ct = 226 MHz
"""
import numpy as np
from triqg.atoms import CsAtom, RbAtom, composite_basis_state
from triqg.pulses import omega_gaussian
from triqg.hamiltonian import build_hamiltonian
from triqg.decoherence import build_collapse_operators
from triqg.solver import simulate
from triqg.analysis import state_fidelity, average_gate_fidelity

from level_parameters import V_ct, V_cc, V_ct_MHz, gamma_r, gamma_R, gamma_P

omega_p_MHz = 50.0
omega_R_MHz = 175.0
omega_c_MHz = 50.0
alpha       = 4.0

omega_p_amp = 2 * np.pi * omega_p_MHz
omega_R_amp = 2 * np.pi * omega_R_MHz
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

deltas = [500.0, 600.0, 700.0, 800.0, 900.0]
ks     = [0.90, 0.93, 0.95, 0.96, 0.97, 0.99, 1.02]

# Header
print(f"V_ct/(2 pi) = {V_ct_MHz:.1f} MHz,  Omega_p=50, OR=175, Oc=50 MHz, alpha={alpha}\n")
print("F_bar matrix  (rows: Delta MHz, cols: k = area / (pi/4))")
print(f"{'Delta\\k':>8s}" + "".join(f"{k:>10.3f}" for k in ks))
print("-" * (8 + 10*len(ks)))

best_F = 0.0; best_cfg = None
for delta_MHz in deltas:
    delta = 2 * np.pi * delta_MHz
    sigma_pi4 = ((2*np.pi*delta)/(omega_p_amp**2 * I_inf))**3
    H = build_hamiltonian(delta, V_ct, pulse_p=omega_gaussian, V_cc=V_cc)
    row_results = []
    for k in ks:
        sigma = (k**3) * sigma_pi4
        T_f = (alpha * sigma) ** (1.0/3.0)
        args = {"omega_c_amp": omega_c_amp, "omega_p_amp": omega_p_amp,
                "omega_R_amp": omega_R_amp, "T_c": T_c, "T_f": T_f, "sigma": sigma}
        tlist = np.linspace(0, 2*T_c + 2*T_f, 350)
        opts  = {"store_final_state": True, "nsteps": 200000}
        pairs = []
        for c1,c2,t,psi_in,psi_id in basis_inputs:
            res = simulate(method="mesolve", H=H, psi0=psi_in, tlist=tlist,
                           c_ops=c_ops, e_ops=[], options=opts, args=args)
            pairs.append((res.final_state, psi_id))
        F = average_gate_fidelity(pairs)
        row_results.append(F)
        if F > best_F:
            best_F = F
            best_cfg = (delta_MHz, k, sigma, T_f)
    print(f"{delta_MHz:>8.0f}" + "".join(f"{F:>10.6f}" for F in row_results))

print()
print(f"Best: F_bar = {best_F:.6f} at Delta = {best_cfg[0]:.0f} MHz, k = {best_cfg[1]:.3f}")
print(f"      sigma = {best_cfg[2]*1e3:.4f} ns, T_f = {best_cfg[3]*1e3:.3f} ns, "
      f"2T_f = {2*best_cfg[3]*1e3:.3f} ns")
print(f"      total gate = {(2*np.pi/(2*np.pi*omega_c_MHz) + 2*best_cfg[3])*1e3:.3f} ns")

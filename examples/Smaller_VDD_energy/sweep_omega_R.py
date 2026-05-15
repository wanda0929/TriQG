"""
2-D sweep over (Delta, Omega_R) at fixed Omega_p, alpha, K.

Goal
----
Test whether the Farouk ratio Omega_R = 3.5 * Omega_p is actually optimal,
or whether relaxing it opens up a wider F_bar > 0.99 plateau (and in
particular unlocks the small-Delta regime, Delta < 400 MHz, that is
currently blockade-leakage-limited).

Physics summary
---------------
With Omega_p fixed, three things change as we vary (Delta, Omega_R):

  1. Two-photon blockade margin  M2 = V_ct / [Omega_p Omega_R / (2 Delta)]
     -> M2 ~ V_ct * 2 Delta / (Omega_p Omega_R).
     Smaller Omega_R restores the same M2 at smaller Delta.

  2. One-photon AC-Stark margin  M1 = V_ct / [Omega_p^2 / (2 Delta)]
     Independent of Omega_R.

  3. EIT dressing of |P>:  Omega_R splits |P> by ~Omega_R, providing
     the "dark state" that distinguishes blockaded vs unblockaded
     dynamics.  Too small -> P-state population builds up
     (P decays at tau_P = 131 ns, the dominant decoherence channel
     in this protocol).  Too large -> EIT shield is so strong that
     even the blockaded case suffers reduced Raman efficiency.

In the present codebase, ``compute_pulse_area`` integrates
Omega_p^2 / (2 Delta), so the sigma-formula
    sigma_pi4 = ((2 pi Delta) / (Omega_p^2 I_inf))^3
does NOT depend on Omega_R.  Pulse duration T_f and shape are
therefore held fixed across the Omega_R sweep at fixed (Delta, Omega_p).
Only Omega_R itself changes -> the EIT control field.

Fixed parameters
----------------
    Omega_p / (2 pi) = 50 MHz
    Omega_c / (2 pi) = 50 MHz
    alpha            = 4 (smooth edges)
    K                = 0.95  (sub-pi/4 optimum found earlier)
    Energy levels    : Rb 54 D_3/2 + Cs 62 D_5/2  ->  V_ct/(2 pi) = 226.27 MHz

Output
------
F_bar matrix over (Delta, ratio = Omega_R / Omega_p).  Per-input
breakdown printed for the best cell.  Best (Delta, ratio) recorded.
"""

import numpy as np

from triqg.atoms import CsAtom, RbAtom, composite_basis_state
from triqg.pulses import omega_gaussian
from triqg.hamiltonian import build_hamiltonian
from triqg.decoherence import build_collapse_operators
from triqg.solver import simulate
from triqg.analysis import state_fidelity, average_gate_fidelity

from level_parameters import V_ct, V_cc, V_ct_MHz, gamma_r, gamma_R, gamma_P

# -----------------------------------------------------------------
# Fixed parameters
# -----------------------------------------------------------------
omega_p_MHz = 50.0
omega_c_MHz = 50.0
alpha       = 4.0
K           = 0.95     # area / (pi/4)

omega_p_amp = 2 * np.pi * omega_p_MHz
omega_c_amp = 2 * np.pi * omega_c_MHz
T_c = np.pi / omega_c_amp

I_inf = 2 * 0.92770 * 2**(-1/6)   # ~ 1.6534

# -----------------------------------------------------------------
# Sweep grid
# -----------------------------------------------------------------
# Extend Delta DOWN to 200 MHz to test the small-Delta regime.
deltas_MHz = [200.0, 300.0, 400.0, 500.0, 700.0, 1000.0]
# Sweep Omega_R / Omega_p ratio from 1.0 to 5.0 (Farouk default = 3.5).
ratios     = [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 5.0]

# -----------------------------------------------------------------
# Decoherence + basis
# -----------------------------------------------------------------
c_ops = build_collapse_operators(gamma_r, gamma_R, gamma_P)

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
            label = f"|{c1}{c2}{['A','B'][t]}>"
            basis_inputs.append((label, psi_in, psi_id))

# -----------------------------------------------------------------
# Header
# -----------------------------------------------------------------
print(f"sweep_omega_R: V_ct/(2 pi) = {V_ct_MHz:.2f} MHz,  Omega_p = {omega_p_MHz} MHz,")
print(f"               Omega_c = {omega_c_MHz} MHz,  alpha = {alpha},  K = {K}")
print(f"               (Sigma fixed by area = K*pi/4 at each Delta; T_f from alpha.)\n")

# -----------------------------------------------------------------
# Sweep
# -----------------------------------------------------------------
F_grid     = np.zeros((len(deltas_MHz), len(ratios)))
inf_grid   = np.zeros_like(F_grid)
M2_grid    = np.zeros_like(F_grid)
Tf_grid    = np.zeros_like(F_grid)
perinput   = {}   # (Delta, ratio) -> [F per input]
best_F     = 0.0
best_cell  = None

print(f"{'Delta':>6s} {'OR/Op':>6s} {'OR':>5s} {'sigma':>7s} {'T_f':>7s} "
      f"{'2Tf+2Tc':>8s} {'M1':>6s} {'M2':>6s} {'F_bar':>9s} {'1-F':>9s}")
print("-" * 90)

for i, delta_MHz in enumerate(deltas_MHz):
    delta = 2 * np.pi * delta_MHz
    H = build_hamiltonian(delta, V_ct, pulse_p=omega_gaussian, V_cc=V_cc)

    # sigma fixed at K*pi/4 (depends only on Omega_p, Delta)
    sigma_pi4 = ((2 * np.pi * delta) / (omega_p_amp**2 * I_inf)) ** 3
    sigma     = (K ** 3) * sigma_pi4
    T_f       = (alpha * sigma) ** (1.0 / 3.0)
    t_total   = 2 * T_c + 2 * T_f
    tlist     = np.linspace(0, t_total, 350)
    opts      = {"store_final_state": True, "nsteps": 200000}

    one_photon_AC = omega_p_MHz**2 / (2 * delta_MHz)
    M1 = V_ct_MHz / one_photon_AC

    for j, ratio in enumerate(ratios):
        omega_R_MHz = ratio * omega_p_MHz
        omega_R_amp = 2 * np.pi * omega_R_MHz
        two_photon_R = omega_p_MHz * omega_R_MHz / (2 * delta_MHz)
        M2 = V_ct_MHz / two_photon_R

        args = {"omega_c_amp": omega_c_amp, "omega_p_amp": omega_p_amp,
                "omega_R_amp": omega_R_amp, "T_c": T_c, "T_f": T_f, "sigma": sigma}

        pairs = []
        per_F = []
        for lbl, psi_in, psi_id in basis_inputs:
            res = simulate(method="mesolve", H=H, psi0=psi_in, tlist=tlist,
                           c_ops=c_ops, e_ops=[], options=opts, args=args)
            pairs.append((res.final_state, psi_id))
            per_F.append(state_fidelity(res.final_state, psi_id))
        F_bar = average_gate_fidelity(pairs)

        F_grid[i, j]   = F_bar
        inf_grid[i, j] = 1.0 - F_bar
        M2_grid[i, j]  = M2
        Tf_grid[i, j]  = T_f * 1e3
        perinput[(delta_MHz, ratio)] = per_F

        print(f"{delta_MHz:>6.0f} {ratio:>6.2f} {omega_R_MHz:>5.0f} "
              f"{sigma*1e3:>7.3f} {T_f*1e3:>7.2f} {t_total*1e3:>8.2f} "
              f"{M1:>6.1f} {M2:>6.2f} {F_bar:>9.6f} {1-F_bar:>9.3e}")

        if F_bar > best_F:
            best_F = F_bar
            best_cell = (delta_MHz, ratio, sigma, T_f, M1, M2, omega_R_MHz)
    print()

# -----------------------------------------------------------------
# F_bar matrix view
# -----------------------------------------------------------------
print("=" * 90)
print("F_bar matrix  (rows: Delta MHz; cols: ratio = Omega_R / Omega_p)")
print(f"{'Delta\\r':>8s}" + "".join(f"{r:>10.2f}" for r in ratios))
print("-" * (8 + 10 * len(ratios)))
for i, d in enumerate(deltas_MHz):
    print(f"{d:>8.0f}" + "".join(f"{F_grid[i,j]:>10.6f}" for j in range(len(ratios))))

print()
print("1 - F_bar matrix  (infidelity)")
print(f"{'Delta\\r':>8s}" + "".join(f"{r:>10.2f}" for r in ratios))
print("-" * (8 + 10 * len(ratios)))
for i, d in enumerate(deltas_MHz):
    print(f"{d:>8.0f}" + "".join(f"{inf_grid[i,j]:>10.3e}" for j in range(len(ratios))))

print()
print("M2 (two-photon blockade margin)  V_ct / [Omega_p Omega_R / (2 Delta)]")
print(f"{'Delta\\r':>8s}" + "".join(f"{r:>10.2f}" for r in ratios))
print("-" * (8 + 10 * len(ratios)))
for i, d in enumerate(deltas_MHz):
    print(f"{d:>8.0f}" + "".join(f"{M2_grid[i,j]:>10.2f}" for j in range(len(ratios))))

# -----------------------------------------------------------------
# Best cell + per-input breakdown
# -----------------------------------------------------------------
print()
print("=" * 90)
delta_b, r_b, sigma_b, Tf_b, M1_b, M2_b, OR_b = best_cell
print(f"Best F_bar = {best_F:.6f}  at  Delta = {delta_b:.0f} MHz,  "
      f"Omega_R/Omega_p = {r_b:.2f}  (Omega_R = {OR_b:.0f} MHz)")
print(f"   sigma = {sigma_b*1e3:.4f} ns,  T_f = {Tf_b*1e3:.3f} ns,  "
      f"total gate = {(2*T_c + 2*Tf_b)*1e3:.2f} ns")
print(f"   M1 = {M1_b:.1f},  M2 = {M2_b:.2f}")
print()
print("Per-input fidelities at the best cell:")
for (lbl, _, _), F in zip(basis_inputs, perinput[(delta_b, r_b)]):
    print(f"   {lbl:<8s}  F = {F:.6f}")

# -----------------------------------------------------------------
# Region with F_bar > 0.99 (robustness check)
# -----------------------------------------------------------------
print()
print("F_bar > 0.99 cells:")
n_robust = 0
for i, d in enumerate(deltas_MHz):
    for j, r in enumerate(ratios):
        if F_grid[i, j] > 0.99:
            print(f"   Delta = {d:>4.0f} MHz,  Omega_R/Omega_p = {r:.2f}   F_bar = {F_grid[i,j]:.6f}")
            n_robust += 1
print(f"Total cells with F_bar > 0.99: {n_robust} / {F_grid.size}")

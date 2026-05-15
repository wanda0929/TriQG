"""
Refinement scan around the surprising small-Delta region found by
sweep_omega_R.py:

    (Delta = 200 MHz, Omega_R/Omega_p = 2.5)  ->  F_bar = 0.9915
    gate time = 166 ns  (vs 385 ns at the Farouk-region optimum)

Three things to check:

  1. Is this a knife-edge?  Map a finer (Delta, ratio) grid around it.
  2. Does the K = 0.95 optimum from the large-Delta regime transfer?
     Re-optimize K at the best cell.
  3. Is the ratio=2.0 dip real (worth understanding) or a numeric artifact?

Pass A:  fine (Delta, ratio) grid at K = 0.95
         Delta  in {150, 175, 200, 225, 250, 275, 300}
         ratio  in {2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5}

Pass B:  K scan at the best (Delta, ratio) of Pass A
         K in {0.80, 0.85, 0.90, 0.93, 0.95, 0.97, 1.00, 1.03}

Pass C:  (ratio, K) refinement at the best Delta of Pass A,
         to check robustness of the optimum.
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
# Common setup
# -----------------------------------------------------------------
omega_p_MHz = 50.0
omega_c_MHz = 50.0
alpha       = 4.0

omega_p_amp = 2 * np.pi * omega_p_MHz
omega_c_amp = 2 * np.pi * omega_c_MHz
T_c = np.pi / omega_c_amp
I_inf = 2 * 0.92770 * 2**(-1/6)

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


def run_one(delta_MHz, omega_R_MHz, K, ntime=300):
    """Run F_bar for one cell, return (F_bar, per_input, T_f_ns, M2)."""
    delta = 2 * np.pi * delta_MHz
    omega_R_amp = 2 * np.pi * omega_R_MHz
    sigma_pi4 = ((2 * np.pi * delta) / (omega_p_amp**2 * I_inf)) ** 3
    sigma = (K ** 3) * sigma_pi4
    T_f = (alpha * sigma) ** (1.0 / 3.0)
    args = {"omega_c_amp": omega_c_amp, "omega_p_amp": omega_p_amp,
            "omega_R_amp": omega_R_amp, "T_c": T_c, "T_f": T_f, "sigma": sigma}
    H = build_hamiltonian(delta, V_ct, pulse_p=omega_gaussian, V_cc=V_cc)
    t_total = 2 * T_c + 2 * T_f
    tlist = np.linspace(0, t_total, ntime)
    opts = {"store_final_state": True, "nsteps": 200000}
    pairs = []
    per_F = []
    for lbl, psi_in, psi_id in basis_inputs:
        res = simulate(method="mesolve", H=H, psi0=psi_in, tlist=tlist,
                       c_ops=c_ops, e_ops=[], options=opts, args=args)
        pairs.append((res.final_state, psi_id))
        per_F.append(state_fidelity(res.final_state, psi_id))
    F_bar = average_gate_fidelity(pairs)
    M2 = V_ct_MHz / (omega_p_MHz * omega_R_MHz / (2 * delta_MHz))
    return F_bar, per_F, T_f * 1e3, M2, t_total * 1e3


# =================================================================
# Pass A: fine (Delta, ratio) grid at K = 0.95
# =================================================================
print("=" * 90)
print("Pass A: fine (Delta, Omega_R/Omega_p) grid at K = 0.95")
print("=" * 90)
deltas_A = [150.0, 175.0, 200.0, 225.0, 250.0, 275.0, 300.0]
ratios_A = [2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5]
K_A = 0.95

F_A = np.zeros((len(deltas_A), len(ratios_A)))
print(f"{'Delta':>6s} {'ratio':>6s} {'OR':>5s} {'T_f':>7s} {'gate':>7s} "
      f"{'M2':>6s} {'F_bar':>9s} {'1-F':>9s}")
print("-" * 78)
best_A_F = 0.0; best_A_cell = None
for i, dm in enumerate(deltas_A):
    for j, r in enumerate(ratios_A):
        ORm = r * omega_p_MHz
        F, perF, Tfn, M2, gn = run_one(dm, ORm, K_A)
        F_A[i, j] = F
        print(f"{dm:>6.0f} {r:>6.2f} {ORm:>5.0f} {Tfn:>7.2f} {gn:>7.2f} "
              f"{M2:>6.2f} {F:>9.6f} {1-F:>9.3e}")
        if F > best_A_F:
            best_A_F = F
            best_A_cell = (dm, r, ORm, perF, gn)
    print()

print("F_bar matrix (rows: Delta, cols: ratio)")
print(f"{'D\\r':>6s}" + "".join(f"{r:>10.2f}" for r in ratios_A))
print("-" * (6 + 10 * len(ratios_A)))
for i, d in enumerate(deltas_A):
    print(f"{d:>6.0f}" + "".join(f"{F_A[i,j]:>10.6f}" for j in range(len(ratios_A))))

dm_b, r_b, OR_b, perF_b, gate_b = best_A_cell
print()
print(f"Pass A best: F_bar = {best_A_F:.6f} at Delta = {dm_b:.0f} MHz, "
      f"ratio = {r_b:.2f}, gate = {gate_b:.2f} ns")
for (lbl, _, _), F in zip(basis_inputs, perF_b):
    print(f"   {lbl:<8s}  F = {F:.6f}")

# =================================================================
# Pass B: K refinement at the best (Delta, ratio) of Pass A
# =================================================================
print()
print("=" * 90)
print(f"Pass B: K refinement at Delta = {dm_b:.0f} MHz, ratio = {r_b:.2f}")
print("=" * 90)
Ks = [0.80, 0.85, 0.90, 0.93, 0.95, 0.97, 1.00, 1.03]
print(f"{'K':>5s} {'T_f':>7s} {'gate':>7s} {'F_bar':>9s} {'1-F':>9s}")
print("-" * 50)
best_K = K_A; best_KF = best_A_F; best_K_perF = perF_b
for K in Ks:
    F, perF, Tfn, M2, gn = run_one(dm_b, OR_b, K)
    print(f"{K:>5.2f} {Tfn:>7.2f} {gn:>7.2f} {F:>9.6f} {1-F:>9.3e}")
    if F > best_KF:
        best_KF = F; best_K = K; best_K_perF = perF

print(f"\nPass B best: F_bar = {best_KF:.6f} at K = {best_K:.2f}")
for (lbl, _, _), F in zip(basis_inputs, best_K_perF):
    print(f"   {lbl:<8s}  F = {F:.6f}")

# =================================================================
# Pass C: (ratio, K) refinement at Delta = dm_b to test robustness
# =================================================================
print()
print("=" * 90)
print(f"Pass C: (ratio, K) refinement at Delta = {dm_b:.0f} MHz")
print("=" * 90)
ratios_C = [2.25, 2.4, 2.5, 2.6, 2.75, 3.0]
Ks_C     = [0.85, 0.90, 0.93, 0.95, 0.97, 1.00]
F_C = np.zeros((len(ratios_C), len(Ks_C)))
print(f"{'ratio\\K':>8s}" + "".join(f"{K:>10.2f}" for K in Ks_C))
print("-" * (8 + 10 * len(Ks_C)))
for i, r in enumerate(ratios_C):
    ORm = r * omega_p_MHz
    row = []
    for K in Ks_C:
        F, _, _, _, _ = run_one(dm_b, ORm, K, ntime=250)
        row.append(F)
    F_C[i, :] = row
    print(f"{r:>8.2f}" + "".join(f"{F:>10.6f}" for F in row))

# Count F > 0.99 plateau
n99 = int(np.sum(F_C > 0.99))
print(f"\nPass C: cells with F_bar > 0.99: {n99} / {F_C.size}")
ii, jj = np.unravel_index(np.argmax(F_C), F_C.shape)
print(f"Pass C best: F_bar = {F_C[ii, jj]:.6f} at ratio = {ratios_C[ii]:.2f}, K = {Ks_C[jj]:.2f}")

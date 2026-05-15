"""
Sweep at Omega_p = 20 MHz (vs the FINAL Omega_p = 50 MHz).

Motivation
----------
The user observed that the adiabatic-elimination small parameter
Omega_p / (2 Delta) hits 12.5% at the small-Delta peak
(Delta = 200 MHz, Omega_p = 50 MHz), which is borderline.
Dropping Omega_p to 20 MHz cuts the small parameter to
20 / (2*200) = 5% -- the same value the FINAL config achieves at
Delta = 500.

Trade-off
---------
At fixed area K*pi/4:
    sigma  prop Delta^3 / Omega_p^6
    T_f    prop Delta   / Omega_p^2
So dropping Omega_p from 50 -> 20 at fixed Delta stretches T_f by
(50/20)^2 = 6.25x.  At Delta = 200:
    Omega_p = 50:  T_f = 73 ns,    gate = 166 ns
    Omega_p = 20:  T_f = 456 ns,   gate = 932 ns
So the win must come from blockade-margin gain at the cost of
Rydberg-decay loss over the longer gate.

Blockade margins at Omega_p = 20, Delta = 200, ratio = 3.5:
    M1 = V_ct / [Omega_p^2/(2 Delta)]       = 226.27 / 1.00   = 226
    M2 = V_ct / [Omega_p Omega_R / (2 Delta)] = 226.27 / 3.50   = 64.6
Both *much* deeper in strong-blockade than the FINAL
(M1=90, M2=26). So small-Omega_p + small-Delta should be
blockade-clean if it works at all.

Decay-error estimate (rough):
    Cs |r> dwell time ~ T_gate ~ 900 ns,  tau_r = 77 us
    -> per-Cs decay = T_gate/tau_r ~ 1.2%
    -> |11*> branch loses 2 * 1.2% ~ 2.4%, set the bar low.

Scan plan
---------
Pass A:  fix Omega_p = 20, K = 0.95, sweep (Delta, ratio)
Pass B:  K refinement at the best (Delta, ratio) of Pass A
Pass C:  side-by-side comparison of Omega_p in {15, 20, 25, 30, 50}
         at the best (Delta, ratio) of Pass A
"""
import numpy as np

from triqg.atoms import CsAtom, RbAtom, composite_basis_state
from triqg.pulses import omega_gaussian
from triqg.hamiltonian import build_hamiltonian
from triqg.decoherence import build_collapse_operators
from triqg.solver import simulate
from triqg.analysis import state_fidelity, average_gate_fidelity

from level_parameters import V_ct, V_cc, V_ct_MHz, gamma_r, gamma_R, gamma_P

omega_c_MHz = 50.0
omega_c_amp = 2 * np.pi * omega_c_MHz
T_c = np.pi / omega_c_amp
alpha = 4.0
I_inf = 2 * 0.92770 * 2**(-1 / 6)

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


def run_one(omega_p_MHz, delta_MHz, ratio, K, ntime=350):
    omega_p_amp = 2 * np.pi * omega_p_MHz
    omega_R_MHz = ratio * omega_p_MHz
    omega_R_amp = 2 * np.pi * omega_R_MHz
    delta = 2 * np.pi * delta_MHz

    sigma_pi4 = ((2 * np.pi * delta) / (omega_p_amp**2 * I_inf)) ** 3
    sigma = (K ** 3) * sigma_pi4
    T_f = (alpha * sigma) ** (1 / 3)
    t_total = 2 * T_c + 2 * T_f

    args = {"omega_c_amp": omega_c_amp, "omega_p_amp": omega_p_amp,
            "omega_R_amp": omega_R_amp, "T_c": T_c, "T_f": T_f, "sigma": sigma}
    H = build_hamiltonian(delta, V_ct, pulse_p=omega_gaussian, V_cc=V_cc)
    tlist = np.linspace(0, t_total, ntime)
    opts = {"store_final_state": True, "nsteps": 400000}

    pairs = []; per_F = []
    for lbl, psi_in, psi_id in basis_inputs:
        res = simulate(method="mesolve", H=H, psi0=psi_in, tlist=tlist,
                       c_ops=c_ops, e_ops=[], options=opts, args=args)
        pairs.append((res.final_state, psi_id))
        per_F.append(state_fidelity(res.final_state, psi_id))
    F_bar = average_gate_fidelity(pairs)
    M1 = V_ct_MHz / (omega_p_MHz ** 2 / (2 * delta_MHz))
    M2 = V_ct_MHz / (omega_p_MHz * omega_R_MHz / (2 * delta_MHz))
    return F_bar, per_F, T_f * 1e3, t_total * 1e3, M1, M2


# =====================================================================
# Pass A:  Omega_p = 20 MHz fixed, scan (Delta, ratio)
# =====================================================================
print("=" * 96)
print("Pass A:  Omega_p = 20 MHz,  K = 0.95,  alpha = 4  (sweep Delta, ratio)")
print("=" * 96)
omega_p_A = 20.0
deltas_A  = [100.0, 150.0, 200.0, 300.0, 500.0, 700.0, 1000.0]
ratios_A  = [2.0, 2.5, 3.0, 3.5, 4.5, 6.0]
K_A = 0.95

F_A = np.zeros((len(deltas_A), len(ratios_A)))
print(f"{'Delta':>6s} {'ratio':>6s} {'OR':>5s} {'T_f':>8s} {'gate':>8s} "
      f"{'M1':>7s} {'M2':>7s} {'Op/2D':>7s} {'F_bar':>9s} {'1-F':>9s}")
print("-" * 96)
best_A_F = 0.0; best_A_cell = None
for i, dm in enumerate(deltas_A):
    for j, r in enumerate(ratios_A):
        F, perF, Tfn, gn, M1, M2 = run_one(omega_p_A, dm, r, K_A)
        F_A[i, j] = F
        Op_2D = omega_p_A / (2 * dm)
        print(f"{dm:>6.0f} {r:>6.2f} {r*omega_p_A:>5.0f} {Tfn:>8.2f} {gn:>8.2f} "
              f"{M1:>7.1f} {M2:>7.2f} {Op_2D:>7.4f} {F:>9.6f} {1-F:>9.3e}")
        if F > best_A_F:
            best_A_F = F; best_A_cell = (dm, r, perF, gn, M1, M2)
    print()

print("F_bar matrix  (rows: Delta MHz, cols: ratio)")
print(f"{'D\\r':>6s}" + "".join(f"{r:>10.2f}" for r in ratios_A))
print("-" * (6 + 10 * len(ratios_A)))
for i, d in enumerate(deltas_A):
    print(f"{d:>6.0f}" + "".join(f"{F_A[i,j]:>10.6f}" for j in range(len(ratios_A))))

print()
print(f"Cells with F_bar > 0.99: {int(np.sum(F_A > 0.99))} / {F_A.size}")
print(f"Cells with F_bar > 0.992: {int(np.sum(F_A > 0.992))} / {F_A.size}")
dm_b, r_b, perF_b, gate_b, M1_b, M2_b = best_A_cell
print(f"Pass A best:  F_bar = {best_A_F:.6f}  at Delta = {dm_b:.0f} MHz, ratio = {r_b:.2f}")
print(f"               gate = {gate_b:.2f} ns,  M1 = {M1_b:.1f},  M2 = {M2_b:.2f}")
for (lbl, _, _), F in zip(basis_inputs, perF_b):
    print(f"   {lbl:<8s}  F = {F:.6f}")

# =====================================================================
# Pass B:  K refinement at best (Delta, ratio) of Pass A
# =====================================================================
print()
print("=" * 96)
print(f"Pass B:  K refinement at Omega_p = 20 MHz, Delta = {dm_b:.0f} MHz, ratio = {r_b:.2f}")
print("=" * 96)
Ks = [0.85, 0.90, 0.93, 0.95, 0.97, 1.00, 1.03, 1.05]
print(f"{'K':>5s} {'T_f':>8s} {'gate':>8s} {'F_bar':>9s} {'1-F':>9s}")
print("-" * 50)
best_KF = best_A_F; best_K = K_A; best_K_perF = perF_b
for K in Ks:
    F, perF, Tfn, gn, M1, M2 = run_one(omega_p_A, dm_b, r_b, K)
    print(f"{K:>5.2f} {Tfn:>8.2f} {gn:>8.2f} {F:>9.6f} {1-F:>9.3e}")
    if F > best_KF:
        best_KF = F; best_K = K; best_K_perF = perF

print(f"\nPass B best:  F_bar = {best_KF:.6f}  at K = {best_K:.2f}")
for (lbl, _, _), F in zip(basis_inputs, best_K_perF):
    print(f"   {lbl:<8s}  F = {F:.6f}")

# =====================================================================
# Pass C:  Omega_p comparison at best (Delta, ratio) of Pass A, K = best_K
# =====================================================================
print()
print("=" * 96)
print(f"Pass C:  Omega_p comparison at Delta = {dm_b:.0f} MHz, ratio = {r_b:.2f}, K = {best_K:.2f}")
print("=" * 96)
omega_p_C = [15.0, 20.0, 25.0, 30.0, 40.0, 50.0]
print(f"{'Op':>4s} {'OR':>5s} {'Op/2D':>7s} {'T_f':>8s} {'gate':>8s} "
      f"{'M1':>7s} {'M2':>7s} {'F_bar':>9s} {'1-F':>9s}")
print("-" * 80)
results_C = []
for op_MHz in omega_p_C:
    F, perF, Tfn, gn, M1, M2 = run_one(op_MHz, dm_b, r_b, best_K)
    Op_2D = op_MHz / (2 * dm_b)
    print(f"{op_MHz:>4.0f} {r_b*op_MHz:>5.0f} {Op_2D:>7.4f} {Tfn:>8.2f} {gn:>8.2f} "
          f"{M1:>7.1f} {M2:>7.2f} {F:>9.6f} {1-F:>9.3e}")
    results_C.append((op_MHz, F, perF, gn))

best_C = max(results_C, key=lambda x: x[1])
print(f"\nPass C best: F_bar = {best_C[1]:.6f}  at Omega_p = {best_C[0]:.0f} MHz "
      f"(gate = {best_C[3]:.2f} ns)")

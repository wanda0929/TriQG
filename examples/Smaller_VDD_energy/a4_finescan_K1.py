"""
Focused fine-grid sweep at a = 4.0 um, K = 1, with ARC-verified lifetimes.

Goal
----
Find the truly robust K = 1 champion at the a = 4 um plateau identified
by `brute_force_or_K1_smallA.py`: that earlier coarse run showed the
broadest pass-rate plateau (56 %) lives at a = 4.0 um, but the
"highest-F_bar" cell at a = 3.5 had robust_score = 0.57 (knife-edge).
Here we zoom into a = 4.0 um at finer resolution on the remaining
four axes to find the cell whose entire 4-D Manhattan-1 neighbourhood
also exceeds 0.99 (robust_score = 1.0).

Lifetime upgrade  (revised: T = 0 K, BBR removed)
-------------------------------------------------
The brute-force runs used n^3-scaled estimates:
    tau_R (Rb 54 D_3/2) = 74 us   (n^3 from paper's n=66, 134.87 us)
    tau_r (Cs 62 D_5/2) = 77 us   (n^3 from paper's n=76, 142.73 us)

ARC 3.10.2 (Sibalic et al. 2017) with includeLevelsUpTo = n+30 at
T = 0 K (spontaneous emission only, no blackbody contribution):
    tau_R (Rb 54 D_3/2) = 164.55 us   (2.22 x longer than estimate)
    tau_r (Cs 62 D_5/2) = 138.87 us   (1.80 x longer)
    tau_P (Rb  7 P_3/2) =   0.270 us  (ARC; 0 K and 300 K differ < 0.5 % at n=7)
This is the intrinsic spontaneous-emission-limited regime; the 300 K
BBR contribution would shorten tau_R, tau_r by roughly 2 x.  At n = 7
BBR is negligible so tau_P is essentially the same as at 300 K.

Grid
----
  Omega_p_MHz  : [50, 55, 60, 65, 70]              (5)
  ratio        : [2.8, 3.0, 3.2]                    (3)
  Delta_MHz    : [400, 450, 500, 550, 600]          (5)
  Omega_c_MHz  : [50, 60, 70]                       (3)

Total 5 * 3 * 5 * 3 = 225 cells, ~ 4 min on Apple Silicon.

Outputs
-------
  a4_finescan_K1.csv               -- all 225 cells
  a4_finescan_K1_robust.csv        -- cells with F_bar > 0.99
  a4_finescan_K1.log               -- this log
"""

import csv
import time
from itertools import product

import numpy as np

from triqg.atoms import CsAtom, RbAtom, composite_basis_state
from triqg.pulses import omega_gaussian
from triqg.hamiltonian import build_hamiltonian
from triqg.decoherence import build_collapse_operators
from triqg.solver import simulate
from triqg.analysis import state_fidelity, average_gate_fidelity

from level_parameters import C3_tilde, C6_RbRb, C6_CsCs

# =====================================================================
# ARC-verified lifetimes (T = 0 K, spontaneous emission only,
#                         includeLevelsUpTo = n + 30)
# =====================================================================
tau_R_us = 164.55   # Rb 54 D_3/2  (ARC, 0 K; was 83.41 us at 300 K)
tau_r_us = 138.87   # Cs 62 D_5/2  (ARC, 0 K; was 83.90 us at 300 K)
tau_P_us = 0.270    # Rb 7 P_3/2   (ARC, 0 K; BBR negligible at n = 7)

gamma_r = 1.0 / tau_r_us
gamma_R = 1.0 / tau_R_us
gamma_P = 1.0 / tau_P_us

# =====================================================================
# Fixed pulse-shape parameters
# =====================================================================
ALPHA = 4.0
K     = 1.00
I_inf = 2 * 0.92770 * 2 ** (-1.0 / 6.0)

# =====================================================================
# Fixed lattice: a = 4.0 um exactly
# =====================================================================
a_um = 4.0
r_DA = a_um / np.sqrt(2)
r_DD = a_um
r_AA = a_um * np.sqrt(2)
V_ct_MHz = 1000.0 * C3_tilde / r_DA ** 3
V_DD_MHz = 1000.0 * C6_RbRb  / r_DD ** 6
V_cc_MHz = 1000.0 * C6_CsCs  / r_AA ** 6
V_ct = 2 * np.pi * V_ct_MHz
V_cc = 2 * np.pi * V_cc_MHz
R_DD = V_ct_MHz / V_DD_MHz
R_AA = V_ct_MHz / abs(V_cc_MHz)

# =====================================================================
# Grids
# =====================================================================
omega_p_grid = [50.0, 55.0, 60.0, 65.0, 70.0]
ratio_grid   = [2.8, 3.0, 3.2]
delta_grid   = [400.0, 450.0, 500.0, 550.0, 600.0]
omega_c_grid = [50.0, 60.0, 70.0]

grids = {
    "omega_p_MHz":  omega_p_grid,
    "ratio":        ratio_grid,
    "delta_MHz":    delta_grid,
    "omega_c_MHz":  omega_c_grid,
}
AXES = list(grids.keys())
SHAPE = tuple(len(grids[a]) for a in AXES)
TOTAL_CELLS = int(np.prod(SHAPE))

# =====================================================================
# Basis states
# =====================================================================
cs = CsAtom()
rb = RbAtom()
ctrl_levels = [cs.level_index["0"], cs.level_index["1"]]
tgt_levels  = [rb.level_index["A"], rb.level_index["B"]]

basis_inputs = []
for c1 in range(2):
    for c2 in range(2):
        for t in range(2):
            psi_in = composite_basis_state(
                ctrl_levels[c1], ctrl_levels[c2], tgt_levels[t]
            )
            psi_id = composite_basis_state(
                ctrl_levels[c1], ctrl_levels[c2],
                tgt_levels[1 - t if (c1 or c2) else t],
            )
            basis_inputs.append((c1, c2, t, psi_in, psi_id))

c_ops = build_collapse_operators(gamma_r, gamma_R, gamma_P)
solver_opts = {"store_final_state": True, "nsteps": 200000}

# =====================================================================
# Header
# =====================================================================
print("=" * 90)
print("a = 4.0 um fine-scan, K = 1, ARC-verified lifetimes")
print("=" * 90)
print(f"Fixed:  a = {a_um:.2f} um,  K = {K},  alpha = {ALPHA}")
print(f"V_ct / (2 pi) = {V_ct_MHz:7.2f} MHz")
print(f"R_DD = {R_DD:5.1f},  R_AA = {R_AA:5.1f}   (both >> 100, selectivity OK)")
print(f"Lifetimes (ARC, 0 K):    tau_R = {tau_R_us:6.2f} us,  "
      f"tau_r = {tau_r_us:6.2f} us,  tau_P = {tau_P_us:6.3f} us")
print(f"Grids:")
for ax in AXES:
    print(f"    {ax:<14s} {grids[ax]}")
print(f"Total cells: {TOTAL_CELLS}")
print()
print(f"{'#':>4} {'Op':>5} {'r':>5} {'D':>5} {'Oc':>5} "
      f"{'Tc':>5} {'Tf':>6} {'Ttot':>6} {'M2':>5} {'F_bar':>9} {'eta':>4}")
print("-" * 80)

# =====================================================================
# Sweep
# =====================================================================
H = build_hamiltonian(2 * np.pi * 500.0, V_ct, pulse_p=omega_gaussian, V_cc=V_cc)
# H is rebuilt per cell because delta changes; placeholder above for clarity

results = []
t0 = time.time()
n = 0

for (op_MHz, ratio, d_MHz, oc_MHz) in product(
    omega_p_grid, ratio_grid, delta_grid, omega_c_grid
):
    n += 1

    omega_p_amp = 2 * np.pi * op_MHz
    omega_R_MHz = ratio * op_MHz
    omega_R_amp = 2 * np.pi * omega_R_MHz
    omega_c_amp = 2 * np.pi * oc_MHz
    delta       = 2 * np.pi * d_MHz

    T_c = np.pi / omega_c_amp
    sigma_pi4 = ((2 * np.pi * delta) / (omega_p_amp ** 2 * I_inf)) ** 3
    sigma = (K ** 3) * sigma_pi4
    T_f   = (ALPHA * sigma) ** (1.0 / 3.0)
    t_total = 2 * T_c + 2 * T_f

    one_photon_AC = op_MHz ** 2 / (2 * d_MHz)
    two_photon_R  = op_MHz * omega_R_MHz / (2 * d_MHz)
    M1 = V_ct_MHz / one_photon_AC
    M2 = V_ct_MHz / two_photon_R

    H = build_hamiltonian(delta, V_ct, pulse_p=omega_gaussian, V_cc=V_cc)
    args = {
        "omega_c_amp": omega_c_amp,
        "omega_p_amp": omega_p_amp,
        "omega_R_amp": omega_R_amp,
        "T_c": T_c, "T_f": T_f, "sigma": sigma,
    }
    tlist = np.linspace(0, t_total, 350)

    pairs = []
    per_F = []
    for c1, c2, t, psi_in, psi_id in basis_inputs:
        res = simulate(
            method="mesolve", H=H, psi0=psi_in, tlist=tlist,
            c_ops=c_ops, e_ops=[], options=solver_opts, args=args,
        )
        pairs.append((res.final_state, psi_id))
        per_F.append(state_fidelity(res.final_state, psi_id))
    F_bar = average_gate_fidelity(pairs)

    F_00 = 0.5 * (per_F[0] + per_F[1])
    F_01 = 0.5 * (per_F[2] + per_F[3])
    F_10 = 0.5 * (per_F[4] + per_F[5])
    F_11 = 0.5 * (per_F[6] + per_F[7])

    row = {
        "omega_p_MHz": op_MHz, "ratio": ratio,
        "delta_MHz": d_MHz, "omega_c_MHz": oc_MHz,
        "omega_R_MHz": omega_R_MHz,
        "T_c_ns":   T_c   * 1e3,
        "T_f_ns":   T_f   * 1e3,
        "sigma_ns": sigma * 1e3,
        "t_total_ns": t_total * 1e3,
        "M1": M1, "M2": M2,
        "F_bar": F_bar,
        "F_00": F_00, "F_01": F_01, "F_10": F_10, "F_11": F_11,
    }
    results.append(row)

    elapsed = time.time() - t0
    eta = elapsed / n * (TOTAL_CELLS - n)
    print(
        f"{n:>4d} {op_MHz:>5.0f} {ratio:>5.2f} {d_MHz:>5.0f} {oc_MHz:>5.0f} "
        f"{T_c*1e3:>5.1f} {T_f*1e3:>6.1f} {t_total*1e3:>6.1f} "
        f"{M2:>5.1f} {F_bar:>9.6f} {eta:>4.0f}s",
        flush=True,
    )

# =====================================================================
# Persist
# =====================================================================
with open("a4_finescan_K1.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(results[0].keys()))
    w.writeheader()
    w.writerows(results)
print(f"\nWrote a4_finescan_K1.csv ({len(results)} rows).")

# =====================================================================
# Robustness scoring
# =====================================================================
def to_idx(row):
    return tuple(grids[a].index(row[a]) for a in AXES)

idx_to_F = {to_idx(r): r["F_bar"] for r in results}

def neighbors(idx):
    out = []
    for i, _ in enumerate(AXES):
        for d in (-1, +1):
            j = idx[i] + d
            if 0 <= j < SHAPE[i]:
                ni = list(idx); ni[i] = j
                out.append(tuple(ni))
    return out

THRESH = 0.99
robust_cells = []
for r in results:
    idx = to_idx(r)
    nbrs = neighbors(idx)
    n_pass_nbrs = sum(1 for ni in nbrs if idx_to_F[ni] > THRESH)
    r["n_neighbors_total"] = len(nbrs)
    r["n_neighbors_pass"]  = n_pass_nbrs
    r["robust_score"]      = (n_pass_nbrs / len(nbrs)) if nbrs else 0.0
    r["pass"]              = int(r["F_bar"] > THRESH)
    if r["pass"]:
        robust_cells.append(r)

robust_cells.sort(
    key=lambda r: (r["robust_score"], r["F_bar"]),
    reverse=True,
)

# =====================================================================
# Summary
# =====================================================================
print("\n" + "=" * 90)
print(f"K = 1, a = 4.0 um:  cells with F_bar > {THRESH}:  "
      f"{len(robust_cells)} / {TOTAL_CELLS}  ({100*len(robust_cells)/TOTAL_CELLS:.1f} %)")
print("=" * 90)

print("\nTop 30 robust cells (ranked by robust_score then F_bar):")
print(f"{'rank':>4} {'Op':>5} {'r':>5} {'D':>5} {'Oc':>5} "
      f"{'F_bar':>9} {'F_11':>9} {'robust':>8} {'Ttot':>6} {'M2':>5}")
print("-" * 90)
for k, r in enumerate(robust_cells[:30], start=1):
    print(f"{k:>4} {r['omega_p_MHz']:>5.0f} {r['ratio']:>5.2f} "
          f"{r['delta_MHz']:>5.0f} {r['omega_c_MHz']:>5.0f} "
          f"{r['F_bar']:>9.6f} {r['F_11']:>9.6f} "
          f"{r['n_neighbors_pass']:>3} / {r['n_neighbors_total']:>3}  "
          f"{r['t_total_ns']:>6.1f} {r['M2']:>5.1f}")

with open("a4_finescan_K1_robust.csv", "w", newline="") as f:
    if robust_cells:
        w = csv.DictWriter(f, fieldnames=list(robust_cells[0].keys()))
        w.writeheader()
        w.writerows(robust_cells)
print(f"\nWrote a4_finescan_K1_robust.csv ({len(robust_cells)} rows).")

# =====================================================================
# Champions
# =====================================================================
if robust_cells:
    print("\n" + "=" * 90)
    print("CHAMPION (highest robust_score, then highest F_bar):")
    print("=" * 90)
    champ = robust_cells[0]
    for k, v in champ.items():
        if isinstance(v, float):
            print(f"   {k:<22s} = {v:.6f}")
        else:
            print(f"   {k:<22s} = {v}")

    print("\nHIGHEST F_bar:")
    bestF = max(robust_cells, key=lambda r: r["F_bar"])
    for k, v in bestF.items():
        if isinstance(v, float):
            print(f"   {k:<22s} = {v:.6f}")
        else:
            print(f"   {k:<22s} = {v}")

    print("\nFASTEST (lowest T_total) with F_bar > 0.99:")
    fastest = min(robust_cells, key=lambda r: r["t_total_ns"])
    for k, v in fastest.items():
        if isinstance(v, float):
            print(f"   {k:<22s} = {v:.6f}")
        else:
            print(f"   {k:<22s} = {v}")

elapsed_total = time.time() - t0
print(f"\nTotal wall time: {elapsed_total:.1f} s ({elapsed_total/60:.1f} min).")

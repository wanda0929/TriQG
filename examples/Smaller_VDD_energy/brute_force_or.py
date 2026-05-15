"""
Brute-force 5-D parameter search for OR-gate average fidelity > 0.99
=====================================================================
Energy levels:  Rb 54 D_{3/2} + Cs 62 D_{5/2}  (level_parameters.py)

Axes searched
-------------
  1.  a_um            -- lattice spacing (mu m)
                          -> V_ct  ~ 1/a^3,  V_DD ~ 1/a^6,  V_cc ~ 1/a^6
  2.  Omega_p_MHz     -- Rb target probe amplitude / (2 pi)
  3.  Omega_R / Omega_p_amp    -- EIT shielding ratio
  4.  Delta_MHz       -- two-photon detuning / (2 pi)
  5.  Omega_c_MHz     -- Cs control Rabi / (2 pi)  -> T_c = pi / Omega_c_amp

Held FIXED at the Smaller_VDD_energy report values (rev. 4 FINAL):
  alpha = T_f^3 / sigma = 4   (smooth super-Gaussian edges, edge/peak ~ 1.1e-7)
  K     = area / (pi/4) = 0.95   (sub-pi/4 area is empirically optimal
                                  at finite V_ct, see rev. 3 report Sec. 4)

What we actually compute per cell
---------------------------------
  Pulse area is set by Omega_p^2 (via I_inf integral, alpha=4 well-resolved
  limit), so sigma is closed-form from (Omega_p, Delta, K) and T_f from
  (sigma, alpha).  T_c is the Cs pi-pulse at omega_c_amp.

  The 8-state computational-basis F_bar is then computed by mesolve.

Robustness scoring
------------------
For each cell C with F_bar(C) > 0.99 we define
    robust_score(C) = (# 1-step nearest neighbors in the 5-D grid
                      that also exceed 0.99) / (# neighbors).
A neighbor is +/- one grid step on exactly ONE axis (Manhattan-1).
The top-N cells are reported, sorted first by pass (F > 0.99),
then by robust_score, then by F_bar.

Cost
----
With the grids below (4 * 3 * 3 * 4 * 3 = 432 cells, 8 mesolve runs
per cell, ~0.15 s each) the total is roughly 432 * 1.2 s ~ 9 minutes
on an Apple Silicon laptop.

Outputs
-------
  brute_force_results.csv     -- every cell, all columns
  brute_force_robust.csv      -- cells with F_bar > 0.99, ranked by robustness
  brute_force_or.log          -- printed progress + summary
"""

import csv
import time
import sys
from itertools import product

import numpy as np

from triqg.atoms import CsAtom, RbAtom, composite_basis_state
from triqg.pulses import omega_gaussian
from triqg.hamiltonian import build_hamiltonian
from triqg.decoherence import build_collapse_operators
from triqg.solver import simulate
from triqg.analysis import state_fidelity, average_gate_fidelity

from level_parameters import C3_tilde, C6_RbRb, C6_CsCs, gamma_r, gamma_R, gamma_P

# =====================================================================
# Fixed pulse-shape parameters (chosen from rev. 3 / rev. 4 analysis)
# =====================================================================
ALPHA = 4.0          # super-Gaussian smoothness  (T_f^3 / sigma)
K     = 0.95         # area / (pi/4)              (sub-pi/4 optimum)
I_inf = 2 * 0.92770 * 2 ** (-1.0 / 6.0)     # ~ 1.6534

# =====================================================================
# Brute-force grids
# =====================================================================
a_um_grid       = [4.0, 4.5, 5.0, 5.5]                    # mu m
omega_p_grid    = [40.0, 50.0, 60.0]                       # MHz / (2 pi)
ratio_grid      = [3.0, 3.5, 4.0]                          # Omega_R / Omega_p
delta_grid      = [400.0, 500.0, 600.0, 800.0]             # MHz / (2 pi)
omega_c_grid    = [30.0, 50.0, 70.0]                       # MHz / (2 pi)

grids = {
    "a_um":         a_um_grid,
    "omega_p_MHz":  omega_p_grid,
    "ratio":        ratio_grid,
    "delta_MHz":    delta_grid,
    "omega_c_MHz":  omega_c_grid,
}
AXES = list(grids.keys())
SHAPE = tuple(len(grids[a]) for a in AXES)
TOTAL_CELLS = int(np.prod(SHAPE))

# =====================================================================
# Basis states (built once)
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
            label = f"|{c1}{c2}{['A','B'][t]}>"
            basis_inputs.append((label, c1, c2, t, psi_in, psi_id))

c_ops = build_collapse_operators(gamma_r, gamma_R, gamma_P)
solver_opts = {"store_final_state": True, "nsteps": 200000}

# =====================================================================
# Header
# =====================================================================
print("=" * 80)
print("Brute-force OR-gate sweep over (a, Omega_p, Omega_R/Omega_p, Delta, Omega_c)")
print("=" * 80)
print(f"Fixed:   K = {K},  alpha = {ALPHA}  (super-Gaussian, edge/peak = "
      f"{np.exp(-ALPHA**2):.2e})")
print(f"Grids:")
for a in AXES:
    print(f"    {a:<14s} {grids[a]}")
print(f"Total cells: {TOTAL_CELLS}   (8 mesolve runs each)")
print(f"Pulse area:  area = K * pi/4 = {K * np.pi / 4:.6f}\n")
print(f"{'#':>5} {'a':>5} {'Op':>5} {'r':>5} {'D':>5} {'Oc':>5} "
      f"{'Vct':>6} {'Vcc':>7} {'Tc':>5} {'Tf':>6} {'Ttot':>6} "
      f"{'M1':>5} {'M2':>5} {'F_bar':>9} {'eta':>5}")
print("-" * 110)

# =====================================================================
# Sweep
# =====================================================================
results = []
t0 = time.time()
n = 0

for (a_um, op_MHz, ratio, d_MHz, oc_MHz) in product(
    a_um_grid, omega_p_grid, ratio_grid, delta_grid, omega_c_grid
):
    n += 1

    # ----- lattice-dependent interactions -----------------------------
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

    # ----- pulse parameters -------------------------------------------
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

    # ----- strong-blockade margins ------------------------------------
    one_photon_AC = op_MHz ** 2 / (2 * d_MHz)
    two_photon_R  = op_MHz * omega_R_MHz / (2 * d_MHz)
    M1 = V_ct_MHz / one_photon_AC
    M2 = V_ct_MHz / two_photon_R

    # ----- Hamiltonian + simulate -------------------------------------
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
    for lbl, c1, c2, t, psi_in, psi_id in basis_inputs:
        res = simulate(
            method="mesolve", H=H, psi0=psi_in, tlist=tlist,
            c_ops=c_ops, e_ops=[], options=solver_opts, args=args,
        )
        pairs.append((res.final_state, psi_id))
        per_F.append(state_fidelity(res.final_state, psi_id))
    F_bar = average_gate_fidelity(pairs)

    # Per-branch (degenerate over target A/B):
    F_00 = 0.5 * (per_F[0] + per_F[1])          # |0,0,*>
    F_01 = 0.5 * (per_F[2] + per_F[3])          # |0,1,*>
    F_10 = 0.5 * (per_F[4] + per_F[5])          # |1,0,*>
    F_11 = 0.5 * (per_F[6] + per_F[7])          # |1,1,*>

    row = {
        "a_um": a_um, "omega_p_MHz": op_MHz, "ratio": ratio,
        "delta_MHz": d_MHz, "omega_c_MHz": oc_MHz,
        "omega_R_MHz": omega_R_MHz,
        "V_ct_MHz": V_ct_MHz, "V_DD_MHz": V_DD_MHz, "V_cc_MHz": V_cc_MHz,
        "R_DD": R_DD, "R_AA": R_AA,
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
        f"{n:>5d} {a_um:>5.2f} {op_MHz:>5.0f} {ratio:>5.2f} {d_MHz:>5.0f} "
        f"{oc_MHz:>5.0f} {V_ct_MHz:>6.1f} {V_cc_MHz:>7.4f} "
        f"{T_c*1e3:>5.1f} {T_f*1e3:>6.1f} {t_total*1e3:>6.1f} "
        f"{M1:>5.1f} {M2:>5.1f} {F_bar:>9.6f} {eta:>5.0f}s",
        flush=True,
    )

# =====================================================================
# Persist raw results
# =====================================================================
csv_cols = list(results[0].keys())
with open("brute_force_results.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=csv_cols)
    w.writeheader()
    w.writerows(results)
print(f"\nWrote brute_force_results.csv ({len(results)} rows).")

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
# Print summary
# =====================================================================
print("\n" + "=" * 90)
print(f"Cells with F_bar > {THRESH}:  {len(robust_cells)} / {TOTAL_CELLS}")
print("=" * 90)
print(f"{'rank':>4} {'a':>5} {'Op':>5} {'r':>5} {'D':>5} {'Oc':>5} "
      f"{'F_bar':>9} {'F_11':>9} {'robust':>7} {'Ttot':>6} {'M2':>5}")
print("-" * 90)
for k, r in enumerate(robust_cells[:30], start=1):
    print(f"{k:>4} {r['a_um']:>5.2f} {r['omega_p_MHz']:>5.0f} "
          f"{r['ratio']:>5.2f} {r['delta_MHz']:>5.0f} {r['omega_c_MHz']:>5.0f} "
          f"{r['F_bar']:>9.6f} {r['F_11']:>9.6f} "
          f"{r['robust_score']:>7.2f} {r['t_total_ns']:>6.1f} {r['M2']:>5.1f}")

with open("brute_force_robust.csv", "w", newline="") as f:
    if robust_cells:
        w = csv.DictWriter(f, fieldnames=list(robust_cells[0].keys()))
        w.writeheader()
        w.writerows(robust_cells)
print(f"\nWrote brute_force_robust.csv ({len(robust_cells)} rows).")

# =====================================================================
# Champions
# =====================================================================
if robust_cells:
    most_robust = robust_cells[0]
    print("\n" + "=" * 90)
    print("MOST ROBUST cell (highest fraction of >0.99 neighbors):")
    for k, v in most_robust.items():
        if isinstance(v, float):
            print(f"   {k:<22s} = {v:.6f}")
        else:
            print(f"   {k:<22s} = {v}")

    fastest = min(
        [r for r in robust_cells],
        key=lambda r: r["t_total_ns"],
    )
    print("\nFASTEST cell with F_bar > 0.99:")
    for k, v in fastest.items():
        if isinstance(v, float):
            print(f"   {k:<22s} = {v:.6f}")
        else:
            print(f"   {k:<22s} = {v}")

    highest_F = max(robust_cells, key=lambda r: r["F_bar"])
    print("\nHIGHEST F_bar cell:")
    for k, v in highest_F.items():
        if isinstance(v, float):
            print(f"   {k:<22s} = {v:.6f}")
        else:
            print(f"   {k:<22s} = {v}")

elapsed_total = time.time() - t0
print(f"\nTotal wall time: {elapsed_total:.1f} s "
      f"({elapsed_total/60:.1f} min)  over {TOTAL_CELLS} cells.")

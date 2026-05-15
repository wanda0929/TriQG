"""
Brute-force 5-D parameter search for OR-gate F_bar > 0.99  -- K = 1, small-a variant
=====================================================================================
Same script structure as brute_force_or_K1.py, but with the lattice
spacing pushed DOWN into the tighter range a in [3.5, 4.5] um instead
of [4.0, 5.5] um.

Motivation
----------
The K = 1 addendum (`Brute_force_OR_K1_addendum.typ`) showed:
  - K = 1 already works at a = 4.0 um (34 / 108 pass) and a = 4.5 um (27 / 108 pass).
  - K = 1 fails completely at a >= 5.0 um (0 / 108 pass).
  - The dominant K = 1 error is OR-protocol leakage ~ 1 / M_2^2 ~ a^6 / V_ct^2.

The physics prediction: tightening a from 4 -> 3.5 um cuts leakage by
(3.5/4)^6 ~ 0.59x while keeping the same-species selectivity budget
R_DD = V_ct/V_DD ~ 3.12 * a^3 safely above 100 (R_DD = 134 at a=3.5).
We expect K = 1 to recover or exceed the K = 0.95 champion fidelity
within this window without paying the selectivity tax that would kick
in at a < 3.5 um.

Grids
-----
  a_um         : [3.5, 3.75, 4.0, 4.25, 4.5]    (5 values)
  Omega_p_MHz  : [50, 60, 70]                    (3 values, biased to high end)
  ratio        : [2.5, 3.0, 3.5]                 (3 values, biased to low end
                                                  per K = 1 addendum finding)
  Delta_MHz    : [400, 500, 600, 800]            (4 values, unchanged)
  Omega_c_MHz  : [50, 70]                        (2 values, drop the bad 30 MHz)

Total 5 * 3 * 3 * 4 * 2 = 360 cells, ~ 7.5 min on Apple Silicon.

Selectivity flag
----------------
Each row records R_DD and R_AA, and the "selectivity_ok" boolean is
True iff *both* R_DD >= 100 and R_AA >= 100.  All five a-values in the
grid above are safely OK, but the flag is logged anyway so any future
extension to a < 3.5 um surfaces the budget violation explicitly.

Outputs:  brute_force_results_K1_smallA.csv,
          brute_force_robust_K1_smallA.csv,
          brute_force_or_K1_smallA.log
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

from level_parameters import C3_tilde, C6_RbRb, C6_CsCs, gamma_r, gamma_R, gamma_P

# =====================================================================
# Fixed pulse-shape parameters
# =====================================================================
ALPHA = 4.0          # super-Gaussian smoothness  (T_f^3 / sigma)
K     = 1.00         # area / (pi/4)  -- canonical paper value
I_inf = 2 * 0.92770 * 2 ** (-1.0 / 6.0)     # ~ 1.6534

# =====================================================================
# Grids:  small-a window, K = 1 biased preferences
# =====================================================================
a_um_grid    = [3.5, 3.75, 4.0, 4.25, 4.5]
omega_p_grid = [50.0, 60.0, 70.0]
ratio_grid   = [2.5, 3.0, 3.5]
delta_grid   = [400.0, 500.0, 600.0, 800.0]
omega_c_grid = [50.0, 70.0]

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
print("Brute-force OR sweep, K = 1, SMALL-a variant   (a in [3.5, 4.5] um)")
print("=" * 90)
print(f"Fixed:   K = {K},  alpha = {ALPHA}")
print(f"Pulse area:  area = K * pi/4 = {K * np.pi/4:.6f}  (canonical paper value)")
print("Grids:")
for ax in AXES:
    print(f"    {ax:<14s} {grids[ax]}")
print(f"Total cells: {TOTAL_CELLS}")
print()
print("Selectivity preview (R_DD = 3.12 a^3, R_AA = 2.49 a^3; budget >= 100):")
for au in a_um_grid:
    rdd = 3.12 * au**3
    raa = 2.49 * au**3
    ok  = "OK" if (rdd >= 100 and raa >= 100) else "FAIL"
    vct = 28284.27 / au**3
    print(f"    a = {au:.2f} um:  V_ct = {vct:7.1f} MHz,  "
          f"R_DD = {rdd:6.1f},  R_AA = {raa:6.1f}   {ok}")
print()
print(f"{'#':>5} {'a':>5} {'Op':>5} {'r':>5} {'D':>5} {'Oc':>5} "
      f"{'Vct':>6} {'RDD':>6} {'RAA':>6} {'Tc':>5} {'Tf':>6} {'Ttot':>6} "
      f"{'M2':>5} {'F_bar':>9} {'eta':>5}")
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
    selectivity_ok = int((R_DD >= 100) and (R_AA >= 100))

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

    # ----- blockade margins -------------------------------------------
    one_photon_AC = op_MHz ** 2 / (2 * d_MHz)
    two_photon_R  = op_MHz * omega_R_MHz / (2 * d_MHz)
    M1 = V_ct_MHz / one_photon_AC
    M2 = V_ct_MHz / two_photon_R

    # ----- simulate ---------------------------------------------------
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
        "a_um": a_um, "omega_p_MHz": op_MHz, "ratio": ratio,
        "delta_MHz": d_MHz, "omega_c_MHz": oc_MHz,
        "omega_R_MHz": omega_R_MHz,
        "V_ct_MHz": V_ct_MHz, "V_DD_MHz": V_DD_MHz, "V_cc_MHz": V_cc_MHz,
        "R_DD": R_DD, "R_AA": R_AA,
        "selectivity_ok": selectivity_ok,
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
        f"{oc_MHz:>5.0f} {V_ct_MHz:>6.1f} {R_DD:>6.1f} {R_AA:>6.1f} "
        f"{T_c*1e3:>5.1f} {T_f*1e3:>6.1f} {t_total*1e3:>6.1f} "
        f"{M2:>5.1f} {F_bar:>9.6f} {eta:>5.0f}s",
        flush=True,
    )

# =====================================================================
# Persist
# =====================================================================
with open("brute_force_results_K1_smallA.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(results[0].keys()))
    w.writeheader()
    w.writerows(results)
print(f"\nWrote brute_force_results_K1_smallA.csv ({len(results)} rows).")

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
# Per-axis pass-rate marginals
# =====================================================================
print("\n" + "=" * 90)
print(f"K = 1, a in [3.5, 4.5]:  cells with F_bar > {THRESH}:  "
      f"{len(robust_cells)} / {TOTAL_CELLS}  ({100*len(robust_cells)/TOTAL_CELLS:.1f} %)")
print("=" * 90)

print("\nPer-axis pass-rate marginal (count / total) and max F_bar in marginal:")
for ax in AXES:
    print(f"\n  Axis: {ax}")
    for val in grids[ax]:
        cells = [r for r in results if r[ax] == val]
        passing = [r for r in cells if r["F_bar"] > THRESH]
        max_F = max(r["F_bar"] for r in cells)
        print(f"    {ax} = {val:<6}  pass: {len(passing):>3} / {len(cells):>3}  "
              f"({100*len(passing)/len(cells):>5.1f} %)  max F_bar = {max_F:.4f}")

# =====================================================================
# Top cells
# =====================================================================
print("\n" + "=" * 100)
print("Top 30 robust cells (ranked by robust_score then F_bar):")
print("=" * 100)
print(f"{'rank':>4} {'a':>5} {'Op':>5} {'r':>5} {'D':>5} {'Oc':>5} "
      f"{'F_bar':>9} {'F_11':>9} {'robust':>7} {'Ttot':>6} {'M2':>5} {'RDD':>5}")
print("-" * 100)
for k, r in enumerate(robust_cells[:30], start=1):
    print(f"{k:>4} {r['a_um']:>5.2f} {r['omega_p_MHz']:>5.0f} "
          f"{r['ratio']:>5.2f} {r['delta_MHz']:>5.0f} {r['omega_c_MHz']:>5.0f} "
          f"{r['F_bar']:>9.6f} {r['F_11']:>9.6f} "
          f"{r['robust_score']:>7.2f} {r['t_total_ns']:>6.1f} {r['M2']:>5.1f} "
          f"{r['R_DD']:>5.0f}")

with open("brute_force_robust_K1_smallA.csv", "w", newline="") as f:
    if robust_cells:
        w = csv.DictWriter(f, fieldnames=list(robust_cells[0].keys()))
        w.writeheader()
        w.writerows(robust_cells)
print(f"\nWrote brute_force_robust_K1_smallA.csv ({len(robust_cells)} rows).")

# =====================================================================
# Champions
# =====================================================================
if robust_cells:
    print("\n" + "=" * 90)
    print("HIGHEST F_bar cell at K = 1 (a in [3.5, 4.5]):")
    print("=" * 90)
    bestF_cell = max(robust_cells, key=lambda r: r["F_bar"])
    for k, v in bestF_cell.items():
        if isinstance(v, float):
            print(f"   {k:<22s} = {v:.6f}")
        else:
            print(f"   {k:<22s} = {v}")

    print("\nFASTEST cell with F_bar > 0.99:")
    fastest = min(robust_cells, key=lambda r: r["t_total_ns"])
    for k, v in fastest.items():
        if isinstance(v, float):
            print(f"   {k:<22s} = {v:.6f}")
        else:
            print(f"   {k:<22s} = {v}")

    print("\nMOST ROBUST cell:")
    most_robust = robust_cells[0]
    for k, v in most_robust.items():
        if isinstance(v, float):
            print(f"   {k:<22s} = {v:.6f}")
        else:
            print(f"   {k:<22s} = {v}")

    # Pareto front
    print("\nPareto front (low T_total, high F_bar) for F_bar > 0.99:")
    sorted_by_T = sorted(robust_cells, key=lambda r: r["t_total_ns"])
    pareto = []
    bestF = 0.0
    for r in sorted_by_T:
        if r["F_bar"] > bestF:
            pareto.append(r); bestF = r["F_bar"]
    for r in pareto:
        print(f"  Tg={r['t_total_ns']:>6.1f}ns  F={r['F_bar']:.4f}  "
              f"a={r['a_um']:.2f} Op={r['omega_p_MHz']:.0f} "
              f"r={r['ratio']:.1f} D={r['delta_MHz']:.0f} "
              f"Oc={r['omega_c_MHz']:.0f}  robust={r['robust_score']:.2f}")
else:
    print("\nNO cell at K = 1 in a in [3.5, 4.5] exceeded F_bar > 0.99!")
    bestF_cell = max(results, key=lambda r: r["F_bar"])
    print(f"\nHighest F_bar achieved: {bestF_cell['F_bar']:.6f}")
    for k, v in bestF_cell.items():
        if isinstance(v, float):
            print(f"   {k:<22s} = {v:.6f}")
        else:
            print(f"   {k:<22s} = {v}")

elapsed_total = time.time() - t0
print(f"\nTotal wall time: {elapsed_total:.1f} s ({elapsed_total/60:.1f} min).")

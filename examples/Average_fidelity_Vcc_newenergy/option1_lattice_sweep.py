"""
Option 1 (magic-spacing) lattice sweep
======================================

Sweep the inter-atomic lattice spacing ``a_um`` from 2.5 to 6.0 um in
steps of 0.1 um and, at each value, recompute the Pedersen Haar-averaged
gate fidelity (raw, *containing the V_cc-induced phase error*) for both
the CCX and OR gates. We also tabulate the phase-corrected and
basis-averaged fidelities and the computational-subspace survival as a
sanity check.

Why this scan:

    The branch-phase analysis predicts that ``F_raw`` should *oscillate*
    with the spacing ``a``, with peaks at the magic condition
    ``|V_cc(a)| * T_wait = 2 pi N``. F_PC (which absorbs all coherent
    diagonal phases) should remain flat at ~0.997, dipping only where
    population leakage from V_cc-induced detuning during the control
    pulses kicks in (the new-finding regime, V_cc ~ Omega_cc).

Output
------
    option1_lattice_sweep.json  -- per-point numerics, both gates
    option1_lattice_sweep.png   -- F_raw and F_PC vs. a, both gates,
                                   with predicted magic-spacing peaks

The script writes the JSON after each completed point so a Ctrl-C
mid-run still leaves usable partial data.
"""

import json
import os
import time
from typing import Callable, Dict, List, Tuple

import numpy as np
from scipy.optimize import minimize

import qutip
from triqg.atoms import CsAtom, RbAtom, DIMS, composite_basis_state
from triqg.pulses import omega_gaussian
from triqg.hamiltonian import build_ccx_hamiltonian, build_hamiltonian
from triqg.decoherence import build_collapse_operators
from triqg.pedersen import (
    build_permutation_unitary,
    choi_on_subspace_via_mesolve,
    pedersen_average_fidelity_from_choi,
)

HERE = os.path.dirname(os.path.abspath(__file__))
JSON_PATH = os.path.join(HERE, "option1_lattice_sweep.json")
PNG_PATH = os.path.join(HERE, "option1_lattice_sweep.png")


# -------------------------------------------------------------------------
# Sweep grid
# -------------------------------------------------------------------------
A_MIN_UM = 2.5
A_MAX_UM = 6.0
A_STEP_UM = 0.1
A_GRID = np.round(np.arange(A_MIN_UM, A_MAX_UM + 0.5 * A_STEP_UM, A_STEP_UM), 3)


# -------------------------------------------------------------------------
# Foerster coefficients (Ireland et al. 2024, Rb 66 D_5/2 + Cs 76 D_3/2)
# -------------------------------------------------------------------------
C3_TILDE = 22.84       # GHz * um^3       (Cs-Rb dipole-dipole, |c1=1>,|t=A> <-> |R>,|r>)
C6_CSCS = -692.9       # GHz * um^6       (Cs-Cs van der Waals, signed)


# -------------------------------------------------------------------------
# Decoherence (ARC, T = 300 K)
# -------------------------------------------------------------------------
GAMMA_R = 1.0 / 142.73    # Cs |r> = |76 D_3/2>
GAMMA_R_RB = 1.0 / 134.87 # Rb |R> = |66 D_5/2>
GAMMA_P = 1.0 / 0.131     # Rb |P> = |7 P_3/2>


# -------------------------------------------------------------------------
# Computational subspace and target unitaries (8-dim, three-qubit)
# -------------------------------------------------------------------------
def _build_subspace() -> Tuple[List[qutip.Qobj], List[str]]:
    cs = CsAtom()
    rb = RbAtom()
    ctrl_levels = [cs.level_index["0"], cs.level_index["1"]]
    tgt_levels = [rb.level_index["A"], rb.level_index["B"]]
    comp_kets, labels = [], []
    for c1 in range(2):
        for c2 in range(2):
            for t in range(2):
                comp_kets.append(
                    composite_basis_state(
                        ctrl_levels[c1], ctrl_levels[c2], tgt_levels[t]
                    )
                )
                labels.append(f"|{c1},{c2},{['A','B'][t]}>")
    return comp_kets, labels


COMP_KETS, LABELS = _build_subspace()


def ccx_perm(idx: int) -> int:
    """CCX (Toffoli): target flips iff both controls are |1>."""
    c1, c2, t = (idx >> 2) & 1, (idx >> 1) & 1, idx & 1
    return (c1 << 2) | (c2 << 1) | (t ^ (c1 & c2))


def or_perm(idx: int) -> int:
    """OR: target flips iff at least one control is |1>."""
    c1, c2, t = (idx >> 2) & 1, (idx >> 1) & 1, idx & 1
    return (c1 << 2) | (c2 << 1) | (t ^ (c1 | c2))


U0_CCX = build_permutation_unitary(COMP_KETS, ccx_perm)
U0_OR = build_permutation_unitary(COMP_KETS, or_perm)


# -------------------------------------------------------------------------
# Phase-corrected Pedersen: max over diagonal D
# -------------------------------------------------------------------------
def _Fpro_with_phases(choi: np.ndarray, U0: np.ndarray, phis: np.ndarray) -> float:
    d = U0.shape[0]
    D = np.diag(np.exp(1j * phis))
    U_target = D @ U0
    F = np.einsum("ki,ijkl,lj->", U_target.conj(), choi, U_target) / (d * d)
    return float(np.real(F))


def phase_corrected_pedersen(
    choi: np.ndarray, U0: np.ndarray,
    n_restarts: int = 8, seed: int = 0,
) -> Tuple[float, float, np.ndarray]:
    d = U0.shape[0]
    rng = np.random.default_rng(seed)
    best_phis = np.zeros(d)
    best_F_pro = _Fpro_with_phases(choi, U0, best_phis)
    for _ in range(n_restarts):
        x0 = rng.uniform(-np.pi, np.pi, size=d)
        res = minimize(
            lambda phis: -_Fpro_with_phases(choi, U0, phis),
            x0, method="L-BFGS-B",
        )
        if -res.fun > best_F_pro:
            best_F_pro = float(-res.fun)
            best_phis = res.x
    F_bar = (d * best_F_pro + 1) / (d + 1)
    return best_F_pro, F_bar, best_phis


# -------------------------------------------------------------------------
# Per-gate Hamiltonian builders. Each returns (H, args, t_total, U0,
# T_wait) so the magic-spacing condition |V_cc| T_wait = 2 pi N can be
# evaluated downstream.
# -------------------------------------------------------------------------
def ccx_problem(V_ct: float, V_cc: float):
    omega_cc_amp = 2 * np.pi * 100
    omega_t_amp = 2 * np.pi * 50
    T_cc = np.pi / omega_cc_amp        # 5 ns
    T_t = np.pi / omega_t_amp          # 10 ns
    H = build_ccx_hamiltonian(V_ct, V_cc=V_cc)
    args = {
        "omega_cc_amp": omega_cc_amp,
        "omega_t_amp": omega_t_amp,
        "T_cc": T_cc,
        "T_t": T_t,
    }
    t_total = 2 * T_cc + 3 * T_t       # 40 ns
    T_wait = 3 * T_t                   # window where both controls sit in |r>
    max_freq = max(omega_cc_amp, omega_t_amp) / (2 * np.pi)
    options = {
        "nsteps": 100000,
        "max_step": 1.0 / (20 * max_freq),
        "atol": 1e-10,
        "rtol": 1e-8,
    }
    return H, args, t_total, U0_CCX, T_wait, options


def or_problem(V_ct: float, V_cc: float):
    omega_c_amp = 2 * np.pi * 50
    omega_p_amp = 2 * np.pi * 50
    omega_R_amp = 3.5 * omega_p_amp
    delta = 2 * np.pi * 500
    T_c = np.pi / omega_c_amp          # 10 ns
    T_f = 0.15                         # 150 ns super-Gaussian half-window
    sigma = 0.001771
    H = build_hamiltonian(delta, V_ct, pulse_p=omega_gaussian, V_cc=V_cc)
    args = {
        "omega_c_amp": omega_c_amp,
        "omega_p_amp": omega_p_amp,
        "omega_R_amp": omega_R_amp,
        "T_c": T_c,
        "T_f": T_f,
        "sigma": sigma,
    }
    t_total = 2 * T_c + 2 * T_f        # 320 ns
    T_wait = 2 * T_f                   # window where both controls sit in |r>
    options = {"nsteps": 100000, "atol": 1e-10, "rtol": 1e-8}
    return H, args, t_total, U0_OR, T_wait, options


GATES: Dict[str, Callable[[float, float], tuple]] = {
    "ccx": ccx_problem,
    "or": or_problem,
}


# -------------------------------------------------------------------------
# Single-point evaluation
# -------------------------------------------------------------------------
def evaluate_point(gate_name: str, a_um: float) -> Dict:
    r_DA = a_um / np.sqrt(2)
    r_AA = a_um * np.sqrt(2)
    V_ct_MHz = 1000.0 * C3_TILDE / r_DA ** 3
    V_cc_MHz = 1000.0 * C6_CSCS / r_AA ** 6
    V_ct = 2 * np.pi * V_ct_MHz
    V_cc = 2 * np.pi * V_cc_MHz

    H, args, t_total, U0, T_wait, options = GATES[gate_name](V_ct, V_cc)
    c_ops = build_collapse_operators(GAMMA_R, GAMMA_R_RB, GAMMA_P)

    t0 = time.time()
    choi = choi_on_subspace_via_mesolve(
        H, c_ops, t_total, COMP_KETS, args=args, options=options, verbose=False,
    )
    dt = time.time() - t0

    F_raw, F_pro_raw, T_P = pedersen_average_fidelity_from_choi(choi, U0)

    # Per-input basis fidelities (Yu eq. 7) read off the Choi diagonal.
    if gate_name == "ccx":
        gate_perm = ccx_perm
    else:
        gate_perm = or_perm
    F_k = [float(np.real(choi[i, i, gate_perm(i), gate_perm(i)])) for i in range(8)]
    F_basis = float(np.mean(F_k))

    F_pro_pc, F_PC, best_phis = phase_corrected_pedersen(choi, U0)

    # Magic-spacing diagnostic: phase accumulated by V_cc on |r,r>
    phi_vcc_pi = (V_cc * T_wait) / np.pi   # in units of pi

    return {
        "gate": gate_name,
        "a_um": float(a_um),
        "r_AA_um": float(r_AA),
        "r_DA_um": float(r_DA),
        "V_ct_MHz": float(V_ct_MHz),
        "V_cc_MHz": float(V_cc_MHz),
        "T_wait_ns": float(T_wait * 1e3),
        "phi_Vcc_over_pi": float(phi_vcc_pi),
        "T_P": float(T_P),
        "F_basis": F_basis,
        "F_k": F_k,
        "F_raw": float(F_raw),
        "F_pro_raw": float(F_pro_raw),
        "F_PC": float(F_PC),
        "F_pro_PC": float(F_pro_pc),
        "phase_correction_over_pi": (best_phis / np.pi).tolist(),
        "wall_time_s": float(dt),
    }


# -------------------------------------------------------------------------
# Magic-spacing predictions (independent, for plot annotation)
# -------------------------------------------------------------------------
def predicted_magic_spacings(T_wait_us: float, n_max: int = 12) -> List[Tuple[int, float]]:
    """Return (N, a_um) pairs satisfying |V_cc(a)| * T_wait = 2 pi N
    that fall inside the swept range [A_MIN_UM, A_MAX_UM]."""
    out = []
    for N in range(1, n_max + 1):
        # |V_cc_MHz| = N / T_wait_us, V_cc_MHz = |C6_MHz_um6| / r_AA^6
        V_cc_MHz_target = N / T_wait_us
        r_AA = (abs(C6_CSCS) * 1000.0 / V_cc_MHz_target) ** (1.0 / 6.0)
        a_um = r_AA / np.sqrt(2)
        if A_MIN_UM <= a_um <= A_MAX_UM:
            out.append((N, float(a_um)))
    return out


# -------------------------------------------------------------------------
# Main sweep, with checkpointing
# -------------------------------------------------------------------------
def load_checkpoint() -> Dict:
    if os.path.exists(JSON_PATH):
        with open(JSON_PATH, "r") as f:
            return json.load(f)
    return {"a_grid": [], "ccx": [], "or": []}


def save_checkpoint(data: Dict) -> None:
    tmp = JSON_PATH + ".tmp"
    with open(tmp, "w") as f:
        json.dump(data, f, indent=2)
    os.replace(tmp, JSON_PATH)


def main():
    print("=" * 76)
    print("Option 1 magic-spacing lattice sweep")
    print(f"  a in [{A_MIN_UM}, {A_MAX_UM}] step {A_STEP_UM} um  "
          f"({len(A_GRID)} points x 2 gates = {2*len(A_GRID)} runs)")
    print(f"  output: {JSON_PATH}")
    print(f"          {PNG_PATH}")
    print("=" * 76)

    data = load_checkpoint()
    if data["a_grid"] and data["a_grid"] != A_GRID.tolist():
        print("[checkpoint] grid changed; starting fresh.")
        data = {"a_grid": [], "ccx": [], "or": []}
    data["a_grid"] = A_GRID.tolist()
    data.setdefault("ccx", [])
    data.setdefault("or", [])

    done_ccx = {round(p["a_um"], 3) for p in data["ccx"]}
    done_or = {round(p["a_um"], 3) for p in data["or"]}

    overall_t0 = time.time()
    for idx, a in enumerate(A_GRID):
        a_key = round(float(a), 3)
        for gate in ("ccx", "or"):
            done_set = done_ccx if gate == "ccx" else done_or
            if a_key in done_set:
                continue
            t0 = time.time()
            res = evaluate_point(gate, a)
            data[gate].append(res)
            data[gate].sort(key=lambda p: p["a_um"])
            save_checkpoint(data)
            print(
                f"[{idx+1:>2}/{len(A_GRID)}] a={a:.2f} um  {gate.upper():<3}  "
                f"V_cc/(2pi)={res['V_cc_MHz']:+8.3f} MHz  "
                f"phi_Vcc={res['phi_Vcc_over_pi']:+6.3f} pi  "
                f"T_P={res['T_P']:.5f}  "
                f"F_raw={res['F_raw']:.5f}  F_PC={res['F_PC']:.5f}  "
                f"({time.time()-t0:.1f} s)"
            )

    total = time.time() - overall_t0
    print(f"\nSweep complete in {total:.1f} s = {total/60:.1f} min")

    # ---------------------------------------------------------------------
    # Plot
    # ---------------------------------------------------------------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, (ax_ccx, ax_or) = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
    for ax, gate, color in [(ax_ccx, "ccx", "C0"), (ax_or, "or", "C1")]:
        pts = sorted(data[gate], key=lambda p: p["a_um"])
        a = np.array([p["a_um"] for p in pts])
        F_raw = np.array([p["F_raw"] for p in pts])
        F_PC = np.array([p["F_PC"] for p in pts])
        F_basis = np.array([p["F_basis"] for p in pts])
        T_P = np.array([p["T_P"] for p in pts])
        ax.plot(a, F_raw, "-o", color=color, lw=1.5, ms=4,
                label=r"$\overline{F}_{\rm raw}$ (with phase error)")
        ax.plot(a, F_PC, "--s", color="black", lw=1, ms=3, alpha=0.8,
                label=r"$\overline{F}_{\rm PC}$ (phase-corrected)")
        ax.plot(a, F_basis, ":", color="grey", lw=1, alpha=0.7,
                label=r"$\overline{F}_{\rm basis}$")
        ax.plot(a, T_P, ":", color="purple", lw=0.7, alpha=0.5,
                label=r"$T_P$")

        # Magic-spacing predictions
        if pts:
            T_wait_us = pts[0]["T_wait_ns"] * 1e-3
            for N, a_pred in predicted_magic_spacings(T_wait_us):
                ax.axvline(a_pred, color="red", lw=0.6, alpha=0.4, ls="-")
                ax.text(a_pred, ax.get_ylim()[1] if False else 0.05,
                        f"N={N}", fontsize=7, ha="center", va="bottom",
                        color="red", alpha=0.7,
                        transform=ax.get_xaxis_transform())

        ax.axhline(1.0, color="grey", lw=0.4, ls="-", alpha=0.3)
        gate_label = "CCX (Toffoli, 40 ns)" if gate == "ccx" else "OR (320 ns)"
        ax.set_title(gate_label)
        ax.set_ylabel("Fidelity")
        ax.grid(True, alpha=0.3)
        ax.set_ylim(-0.05, 1.05)
        ax.legend(loc="lower right", fontsize=8, ncol=2)

    ax_or.set_xlabel(r"Lattice spacing $a$ (μm)")
    fig.suptitle(
        "Option 1 magic-spacing scan: "
        r"raw Pedersen fidelity oscillates with $|V_{cc}|\,T_{\rm wait}=2\pi N$",
        fontsize=11,
    )
    fig.tight_layout()
    fig.savefig(PNG_PATH, dpi=150)
    print(f"plot saved to {PNG_PATH}")


if __name__ == "__main__":
    main()

"""
Hamiltonian construction for Rydberg three-qubit gates.

Builds the full 36-dim time-dependent Hamiltonian in the rotating frame
(RWA) as a QuTiP-compatible list ``[H_static, [H1, c1], ...]``.

Provides:
    build_hamiltonian     -- OR gate Hamiltonian
    build_ccx_hamiltonian -- CCX (Toffoli) gate Hamiltonian
"""

from __future__ import annotations

from typing import List

import qutip

from .atoms import CsAtom, RbAtom, DIMS, composite_projector
from .pulses import omega_c, omega_p, omega_R, omega_cc, omega_t1, omega_t2


def _transition_op(subsystem: int, i: int, j: int) -> qutip.Qobj:
    """Build |i><j| on one subsystem, tensored with identity on the others."""
    ops = []
    for s, dim in enumerate(DIMS):
        if s == subsystem:
            ops.append(qutip.basis(dim, i) * qutip.basis(dim, j).dag())
        else:
            ops.append(qutip.qeye(dim))
    return qutip.tensor(ops)


def build_hamiltonian(delta: float, V_ct: float) -> list:
    """
    Build the time-dependent Hamiltonian for the OR gate.

    Parameters
    ----------
    delta : float
        Detuning of the |A>,|B> -> |P> transition on the target.
    V_ct : float
        Rydberg blockade interaction strength between each control
        and the target.

    Returns
    -------
    list
        QuTiP-compatible Hamiltonian:
        ``[H_static, [H_c1, omega_c], [H_c2, omega_c],
           [H_p, omega_p], [H_R, omega_R]]``
    """
    cs = CsAtom()
    rb = RbAtom()

    # Level indices
    idx_1 = cs.level_index["1"]  # 1
    idx_r = cs.level_index["r"]  # 2
    idx_A = rb.level_index["A"]  # 0
    idx_B = rb.level_index["B"]  # 1
    idx_P = rb.level_index["P"]  # 2
    idx_R = rb.level_index["R"]  # 3

    # --- Static Hamiltonian ---
    # Detuning: delta * I_c1 x I_c2 x |P><P|
    H_detuning = delta * composite_projector(2, idx_P)

    # Rydberg interaction: V_ct * |r><r|_ci x |R><R|_t  for each control
    # composite_projector gives |level><level| x I_others on one subsystem,
    # but here we need the product on two subsystems simultaneously.
    proj_r_c1 = composite_projector(0, idx_r)  # |r><r|_c1 x I_c2 x I_t
    proj_r_c2 = composite_projector(1, idx_r)  # I_c1 x |r><r|_c2 x I_t
    proj_R_t = composite_projector(2, idx_R)  # I_c1 x I_c2 x |R><R|_t

    H_rydberg_c1 = V_ct * proj_r_c1 * proj_R_t
    H_rydberg_c2 = V_ct * proj_r_c2 * proj_R_t

    H_static = H_detuning + H_rydberg_c1 + H_rydberg_c2

    # --- Drive operators (Hermitian coupling terms) ---
    # Control 1: |1><r| + |r><1| on subsystem 0
    H_drive_c1 = _transition_op(0, idx_1, idx_r) + _transition_op(0, idx_r, idx_1)

    # Control 2: |1><r| + |r><1| on subsystem 1
    H_drive_c2 = _transition_op(1, idx_1, idx_r) + _transition_op(1, idx_r, idx_1)

    # Target probe: (|A><P| + |P><A|) + (|B><P| + |P><B|) on subsystem 2
    H_drive_p = (
        _transition_op(2, idx_A, idx_P)
        + _transition_op(2, idx_P, idx_A)
        + _transition_op(2, idx_B, idx_P)
        + _transition_op(2, idx_P, idx_B)
    )

    # Target Rydberg coupling: |P><R| + |R><P| on subsystem 2
    H_drive_R = _transition_op(2, idx_P, idx_R) + _transition_op(2, idx_R, idx_P)

    return [
        H_static,
        [H_drive_c1, omega_c],
        [H_drive_c2, omega_c],
        [H_drive_p, omega_p],
        [H_drive_R, omega_R],
    ]


def build_ccx_hamiltonian(V_ct: float) -> list:
    """
    Build the time-dependent Hamiltonian for the CCX (Toffoli) gate.

    Parameters
    ----------
    V_ct : float
        Rydberg blockade interaction strength between each control
        and the target.

    Returns
    -------
    list
        QuTiP-compatible Hamiltonian:
        ``[H_static, [H_cc, omega_cc], [H_t1, omega_t1], [H_t2, omega_t2]]``
    """
    cs = CsAtom()
    rb = RbAtom()

    # Level indices
    idx_0 = cs.level_index["0"]  # 0
    idx_r = cs.level_index["r"]  # 2
    idx_A = rb.level_index["A"]  # 0
    idx_B = rb.level_index["B"]  # 1
    idx_R = rb.level_index["R"]  # 3

    # --- Static Hamiltonian: Rydberg blockade only (no detuning) ---
    proj_r_c1 = composite_projector(0, idx_r)
    proj_r_c2 = composite_projector(1, idx_r)
    proj_R_t = composite_projector(2, idx_R)

    H_static = V_ct * proj_r_c1 * proj_R_t + V_ct * proj_r_c2 * proj_R_t

    # --- Drive operators ---
    # Control: |0><r| + |r><0| on each control subsystem
    H_drive_cc = (
        _transition_op(0, idx_0, idx_r)
        + _transition_op(0, idx_r, idx_0)
        + _transition_op(1, idx_0, idx_r)
        + _transition_op(1, idx_r, idx_0)
    )

    # Target drive 1: |B><R| + |R><B| on subsystem 2
    H_drive_t1 = _transition_op(2, idx_B, idx_R) + _transition_op(2, idx_R, idx_B)

    # Target drive 2: |A><R| + |R><A| on subsystem 2
    H_drive_t2 = _transition_op(2, idx_A, idx_R) + _transition_op(2, idx_R, idx_A)

    return [
        H_static,
        [H_drive_cc, omega_cc],
        [H_drive_t1, omega_t1],
        [H_drive_t2, omega_t2],
    ]

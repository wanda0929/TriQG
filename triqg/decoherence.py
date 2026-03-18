"""
Collapse operators for decay channels in the Rydberg gate simulation.

Builds a flat list of collapse operators embedded in the 36-dim
composite Hilbert space, suitable for passing to ``mesolve`` / ``mcsolve``.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Tuple

import numpy as np
import qutip

from .atoms import CsAtom, RbAtom, DIMS


def _collapse_op(
    subsystem: int, final: int, initial: int, rate: float, branching: float
) -> qutip.Qobj:
    """
    Build sqrt(rate * branching) * |final><initial| on one subsystem,
    tensored with identity on the others.
    """
    ops = []
    for s, dim in enumerate(DIMS):
        if s == subsystem:
            ops.append(qutip.basis(dim, final) * qutip.basis(dim, initial).dag())
        else:
            ops.append(qutip.qeye(dim))
    return np.sqrt(rate * branching) * qutip.tensor(ops)


def build_collapse_operators(
    gamma_r: float,
    gamma_R: float,
    gamma_P: float,
    branching: Optional[Dict[str, Tuple[float, float]]] = None,
) -> List[qutip.Qobj]:
    """
    Build all collapse operators for the three-atom system.

    Parameters
    ----------
    gamma_r : float
        Decay rate of control Rydberg state |r> (1 / lifetime).
    gamma_R : float
        Decay rate of target Rydberg state |R>.
    gamma_P : float
        Decay rate of target intermediate state |P>.
    branching : dict, optional
        Override branching ratios. Keys: "r", "R", "P".
        Values: (ratio_to_first_ground, ratio_to_second_ground).
        Default: equal branching (0.5, 0.5) for all channels.

    Returns
    -------
    list of qutip.Qobj
        8 collapse operators (36x36 each).
    """
    cs = CsAtom()
    rb = RbAtom()

    idx_0 = cs.level_index["0"]
    idx_1 = cs.level_index["1"]
    idx_r = cs.level_index["r"]
    idx_A = rb.level_index["A"]
    idx_B = rb.level_index["B"]
    idx_P = rb.level_index["P"]
    idx_R = rb.level_index["R"]

    br = branching or {}
    br_r = br.get("r", (0.5, 0.5))
    br_R = br.get("R", (0.5, 0.5))
    br_P = br.get("P", (0.5, 0.5))

    c_ops = []

    # Control 1: |r> -> |0>, |r> -> |1>
    c_ops.append(_collapse_op(0, idx_0, idx_r, gamma_r, br_r[0]))
    c_ops.append(_collapse_op(0, idx_1, idx_r, gamma_r, br_r[1]))

    # Control 2: |r> -> |0>, |r> -> |1>
    c_ops.append(_collapse_op(1, idx_0, idx_r, gamma_r, br_r[0]))
    c_ops.append(_collapse_op(1, idx_1, idx_r, gamma_r, br_r[1]))

    # Target: |R> -> |A>, |R> -> |B>
    c_ops.append(_collapse_op(2, idx_A, idx_R, gamma_R, br_R[0]))
    c_ops.append(_collapse_op(2, idx_B, idx_R, gamma_R, br_R[1]))

    # Target: |P> -> |A>, |P> -> |B>
    c_ops.append(_collapse_op(2, idx_A, idx_P, gamma_P, br_P[0]))
    c_ops.append(_collapse_op(2, idx_B, idx_P, gamma_P, br_P[1]))

    return c_ops

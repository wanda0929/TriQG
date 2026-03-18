"""
Atomic level structures and composite Hilbert space for Rydberg gate simulations.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List

import qutip


@dataclass(frozen=True)
class Atom:
    """An atom with a fixed number of energy levels."""

    num_levels: int
    labels: List[str]
    level_index: Dict[str, int] = field(init=False, repr=False)

    def __post_init__(self):
        # frozen=True requires object.__setattr__ for init
        object.__setattr__(
            self,
            "level_index",
            {label: i for i, label in enumerate(self.labels)},
        )


def CsAtom() -> Atom:
    """133-Cs control atom: |0>, |1>, |r>."""
    return Atom(num_levels=3, labels=["0", "1", "r"])


def RbAtom() -> Atom:
    """87-Rb target atom: |A>, |B>, |P>, |R>."""
    return Atom(num_levels=4, labels=["A", "B", "P", "R"])


# Composite Hilbert space dimensions: control1 x control2 x target
DIMS = [3, 3, 4]


def composite_basis_state(c1: int, c2: int, target: int) -> qutip.Qobj:
    """
    Construct a basis state in the composite 3x3x4 Hilbert space.

    Parameters
    ----------
    c1 : int
        Level index for control atom 1 (0, 1, or 2).
    c2 : int
        Level index for control atom 2 (0, 1, or 2).
    target : int
        Level index for the target atom (0, 1, 2, or 3).

    Returns
    -------
    qutip.Qobj
        A ket in the 36-dimensional composite space.
    """
    return qutip.tensor(
        qutip.basis(DIMS[0], c1),
        qutip.basis(DIMS[1], c2),
        qutip.basis(DIMS[2], target),
    )


def composite_projector(subsystem: int, level: int) -> qutip.Qobj:
    """
    Construct |level><level| for one subsystem, tensored with identity on others.

    Parameters
    ----------
    subsystem : int
        0 = control1, 1 = control2, 2 = target.
    level : int
        Level index within the subsystem.

    Returns
    -------
    qutip.Qobj
        A 36x36 projector operator in the composite space.
    """
    ops = []
    for i, dim in enumerate(DIMS):
        if i == subsystem:
            ket = qutip.basis(dim, level)
            ops.append(ket * ket.dag())
        else:
            ops.append(qutip.qeye(dim))
    return qutip.tensor(ops)

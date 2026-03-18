"""
Analysis utilities for Rydberg gate simulations.

Provides state fidelity computation and population extraction
from QuTiP solver results.
"""

from __future__ import annotations

from typing import List

import numpy as np
import qutip


def state_fidelity(final: qutip.Qobj, target: qutip.Qobj) -> float:
    """
    Compute the fidelity between two quantum states.

    For kets: F = |<target|final>|^2.
    For density matrices: F = (Tr sqrt(sqrt(rho_t) rho_f sqrt(rho_t)))^2,
    computed via ``qutip.fidelity``.

    Parameters
    ----------
    final : qutip.Qobj
        The state to evaluate (ket or density matrix).
    target : qutip.Qobj
        The reference state (ket or density matrix).

    Returns
    -------
    float
        Fidelity in [0, 1].
    """
    # Convert kets to density matrices for uniform handling
    rho_f = qutip.ket2dm(final) if final.isket else final
    rho_t = qutip.ket2dm(target) if target.isket else target

    return float(qutip.fidelity(rho_f, rho_t) ** 2)


def extract_populations(result, num_levels: int) -> np.ndarray:
    """
    Extract population data from a QuTiP solver result.

    Assumes the first ``num_levels`` entries in ``result.expect`` are
    expectation values of population projectors (diagonal elements).

    Parameters
    ----------
    result : QuTiP Result or any object with ``.expect`` attribute
        Solver result. ``result.expect`` should be a list of arrays,
        where each array has length ``len(tlist)``.
    num_levels : int
        Number of population levels to extract.

    Returns
    -------
    np.ndarray
        Shape ``(num_levels, len(tlist))`` array of real populations.
    """
    return np.array([np.real(result.expect[i]) for i in range(num_levels)])

"""
Unified solver interface wrapping QuTiP's mesolve and mcsolve.

Provides a single ``simulate()`` entry point that dispatches to the
appropriate solver and returns a ``SimulationResult`` with a common API.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Union

import numpy as np
import qutip


class SimulationResult:
    """
    Thin wrapper around QuTiP solver results with a uniform interface.

    Attributes
    ----------
    times : np.ndarray
    expect : list of np.ndarray
        For mesolve: direct expect. For mcsolve: average_expect.
    final_state : qutip.Qobj or None
    method : str
        "mesolve" or "mcsolve".

    MC-specific (None for mesolve):
    std_expect, num_trajectories, col_times, col_which
    """

    def __init__(self, raw_result, method: str):
        self._raw = raw_result
        self.method = method

    @property
    def times(self) -> np.ndarray:
        return np.asarray(self._raw.times)

    @property
    def expect(self) -> list:
        if self.method == "mcsolve":
            return self._raw.average_expect
        return self._raw.expect

    @property
    def final_state(self) -> Optional[qutip.Qobj]:
        if self.method == "mcsolve":
            return self._raw.average_final_state
        return self._raw.final_state

    # -- MC-specific attributes --

    @property
    def std_expect(self) -> Optional[list]:
        if self.method == "mcsolve":
            return self._raw.std_expect
        return None

    @property
    def num_trajectories(self) -> Optional[int]:
        if self.method == "mcsolve":
            return self._raw.num_trajectories
        return None

    @property
    def col_times(self) -> Optional[list]:
        if self.method == "mcsolve":
            return self._raw.col_times
        return None

    @property
    def col_which(self) -> Optional[list]:
        if self.method == "mcsolve":
            return self._raw.col_which
        return None

    @property
    def raw(self):
        """Access the underlying QuTiP result object."""
        return self._raw


def simulate(
    method: str,
    H,
    psi0: qutip.Qobj,
    tlist: np.ndarray,
    c_ops: Optional[list] = None,
    e_ops: Optional[list] = None,
    options: Optional[dict] = None,
    args: Optional[dict] = None,
    **kwargs,
) -> SimulationResult:
    """
    Run a simulation using mesolve or mcsolve.

    Parameters
    ----------
    method : str
        "mesolve" or "mcsolve".
    H : Qobj or list
        Hamiltonian (static or time-dependent list).
    psi0 : Qobj
        Initial state (ket or density matrix).
    tlist : array-like
        Time points.
    c_ops : list, optional
        Collapse operators.
    e_ops : list, optional
        Expectation value operators.
    options : dict, optional
        Solver options forwarded to QuTiP.
    args : dict, optional
        Arguments forwarded to time-dependent coefficient functions.
    **kwargs
        Additional arguments. For mcsolve: ``ntraj``, ``seeds``,
        ``target_tol``, ``timeout``.

    Returns
    -------
    SimulationResult
    """
    c_ops = c_ops or []
    e_ops = e_ops or []

    if method == "mesolve":
        raw = qutip.mesolve(
            H, psi0, tlist, c_ops=c_ops, e_ops=e_ops, options=options, args=args
        )
    elif method == "mcsolve":
        ntraj = kwargs.pop("ntraj", 100)
        raw = qutip.mcsolve(
            H,
            psi0,
            tlist,
            c_ops=c_ops,
            e_ops=e_ops,
            ntraj=ntraj,
            options=options,
            args=args,
            **kwargs,
        )
    else:
        raise ValueError(f"Unknown method '{method}'. Use 'mesolve' or 'mcsolve'.")

    return SimulationResult(raw, method=method)

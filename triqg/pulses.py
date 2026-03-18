"""
Pulse shape functions for Rydberg gate simulations.

All pulse functions have signature ``f(t, args) -> float``, compatible
with QuTiP's ``[Qobj, func]`` time-dependent coefficient format.

Expected keys in ``args``:
    omega_c_amp, omega_p_amp, omega_R_amp, T_c, T_f, sigma
"""

from __future__ import annotations

from typing import Callable

import numpy as np
from scipy import integrate


def omega_c(t: float, args: dict) -> float:
    """
    Piecewise square control pulse on |1> <-> |r>.

    Segments:
        [0, T_c)            : +omega_c_amp / 2
        [T_c, T_c + T_f)    : 0
        [T_c+T_f, 2*T_c+T_f): -omega_c_amp / 2
        otherwise            : 0
    """
    amp = args["omega_c_amp"]
    T_c = args["T_c"]
    T_f = args["T_f"]

    if 0 <= t < T_c:
        return amp / 2
    elif T_c + T_f <= t < 2 * T_c + T_f:
        return -amp / 2
    else:
        return 0.0


def omega_p(t: float, args: dict) -> float:
    """
    Cubic super-Gaussian target pulse on |A>,|B> <-> |P>.

    Active during [T_c, T_c + T_f).
    Shape: (omega_p_amp / 2) * exp(-((t - center)^3 / sigma)^2)
    where center = T_c + T_f / 2.
    """
    amp = args["omega_p_amp"]
    T_c = args["T_c"]
    T_f = args["T_f"]
    sigma = args["sigma"]

    if T_c <= t < T_c + T_f:
        center = T_c + T_f / 2
        return (amp / 2) * np.exp(-((((t - center) ** 3) / sigma) ** 2))
    else:
        return 0.0


def omega_R(t: float, args: dict) -> float:
    """
    Constant resonant coupling on |P> <-> |R> during the target window.

    Active during [T_c, T_c + T_f), returns omega_R_amp / 2.
    """
    T_c = args["T_c"]
    T_f = args["T_f"]
    amp = args["omega_R_amp"]

    if T_c <= t < T_c + T_f:
        return amp / 2
    else:
        return 0.0


def compute_pulse_area(
    pulse_func: Callable,
    delta: float,
    t_start: float,
    t_end: float,
    args: dict,
) -> float:
    """
    Numerically integrate pulse_func(t, args)^2 / (2 * delta) over [t_start, t_end].

    This computes the effective two-photon Rabi angle for an off-resonant
    transition with detuning ``delta``.

    Parameters
    ----------
    pulse_func : callable
        Pulse function with signature ``f(t, args) -> float``.
    delta : float
        Detuning (same units as pulse amplitude).
    t_start, t_end : float
        Integration bounds.
    args : dict
        Arguments forwarded to ``pulse_func``.

    Returns
    -------
    float
        The integrated area (should equal pi for a pi-pulse).
    """

    def integrand(t):
        val = pulse_func(t, args)
        return val**2 / (2 * delta)

    result, _ = integrate.quad(integrand, t_start, t_end)
    return result

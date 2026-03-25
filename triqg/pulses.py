"""
Pulse shape functions for Rydberg gate simulations.

All pulse functions have signature ``f(t, args) -> float``, compatible
with QuTiP's ``[Qobj, func]`` time-dependent coefficient format.

Expected keys in ``args``:
    OR gate:  omega_c_amp, omega_p_amp, omega_R_amp, T_c, T_f
    CCX gate: omega_cc_amp, omega_t_amp, T_cc, T_t
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
    elif T_c + 2 * T_f <= t < 2 * T_c + 2 * T_f:
        return -amp / 2
    else:
        return 0.0


def omega_p(t: float, args: dict) -> float:
    """
    Hanning (sin²) target pulse on |A>,|B> <-> |P>.

    Active during [T_c, T_c + 2*T_f].
    Shape: omega_p_amp * sin²(pi * (t - T_c) / (2 * T_f))

    Starts and ends exactly at zero, symmetric, smooth.
    Area efficiency: 3/8 of the rectangular-pulse area.
    """
    amp = args["omega_p_amp"]
    T_c = args["T_c"]
    T_f = args["T_f"]

    if T_c <= t <= T_c + 2 * T_f and T_f > 0:
        return (amp / 2) * np.sin(np.pi * (t - T_c) / (2 * T_f)) ** 2
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

    if 0 <= t < (T_c + T_f) * 2:
        return amp / 2
    else:
        return 0.0


def omega_gaussian(t: float, args: dict) -> float:
    """
    Gaussian target pulse on |A>,|B> <-> |P> with super-Gaussian profile.

    Active during [T_c, T_c + 2*T_f].
    Shape: omega_p_amp * exp(-((t - T_c - T_f)^3 / sigma)^2) / 2

    This is a super-Gaussian of order 6, centered at T_c + T_f.
    The cubic term in the exponent creates a flatter top than a standard Gaussian.

    Expected keys in args:
        omega_p_amp : amplitude coefficient
        T_c         : start time of control pulse
        T_f         : half-width of the target pulse
        sigma       : width parameter controlling the pulse shape
    """
    amp = args["omega_p_amp"]
    T_c = args["T_c"]
    T_f = args["T_f"]
    sigma = args["sigma"]

    if T_c <= t < T_c + 2 * T_f and T_f > 0:
        t_center = T_c + T_f
        exponent = -((t - t_center) ** 3 / sigma) ** 2
        return (amp /2) * np.exp(exponent)
    else:
        return 0.0


def omega_cc(t: float, args: dict) -> float:
    """
    Piecewise square control pulse on |0> <-> |r> for the CCX gate.

    Segments:
        [0, T_cc)                          : +omega_cc_amp / 2
        [T_cc, T_cc + 3*T_t)              : 0  (gap for target pulses)
        [T_cc + 3*T_t, 2*T_cc + 3*T_t)   : -omega_cc_amp / 2
        otherwise                          : 0
    """
    amp = args["omega_cc_amp"]
    T_cc = args["T_cc"]
    T_t = args["T_t"]

    if 0 <= t < T_cc:
        return amp / 2
    elif T_cc + 3 * T_t <= t < 2 * T_cc + 3 * T_t:
        return -amp / 2
    else:
        return 0.0


def omega_t1(t: float, args: dict) -> float:
    """
    CCX target pulse driving |B> <-> |R> (sub-pulses 1 and 3).

    Segments:
        [T_cc, T_cc + T_t)              : omega_t_amp / 2  (sub-pulse 1)
        [T_cc + 2*T_t, T_cc + 3*T_t)   : omega_t_amp / 2  (sub-pulse 3)
        otherwise                        : 0
    """
    amp = args["omega_t_amp"]
    T_cc = args["T_cc"]
    T_t = args["T_t"]

    if T_cc <= t < T_cc + T_t:
        return amp / 2
    elif T_cc + 2 * T_t <= t < T_cc + 3 * T_t:
        return amp / 2
    else:
        return 0.0


def omega_t2(t: float, args: dict) -> float:
    """
    CCX target pulse driving |A> <-> |R> (sub-pulse 2).

    Segments:
        [T_cc + T_t, T_cc + 2*T_t)  : omega_t_amp / 2
        otherwise                     : 0
    """
    amp = args["omega_t_amp"]
    T_cc = args["T_cc"]
    T_t = args["T_t"]

    if T_cc + T_t <= t < T_cc + 2 * T_t:
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


def compute_two_photon_pulse_area(
    pulse_func_p: Callable,
    pulse_func_R: Callable,
    delta: float,
    t_start: float,
    t_end: float,
    args: dict,
) -> float:
    """
    Numerically integrate sqrt(2) * Omega_p(t) * Omega_R(t) / (2 * delta)
    over [t_start, t_end].

    This computes the effective two-photon Rabi angle for an off-resonant
    Raman transition driven by two fields with detuning ``delta``.

    Parameters
    ----------
    pulse_func_p : callable
        Pulse function for Omega_p with signature ``f(t, args) -> float``.
    pulse_func_R : callable
        Pulse function for Omega_R with signature ``f(t, args) -> float``.
    delta : float
        Detuning (same units as pulse amplitudes).
    t_start, t_end : float
        Integration bounds.
    args : dict
        Arguments forwarded to both pulse functions.

    Returns
    -------
    float
        The integrated two-photon pulse area.
    """

    def integrand(t):
        val_p = pulse_func_p(t, args)
        val_R = pulse_func_R(t, args)
        return np.sqrt(2) * val_p * val_R / (2 * delta)

    result, _ = integrate.quad(integrand, t_start, t_end)
    return result

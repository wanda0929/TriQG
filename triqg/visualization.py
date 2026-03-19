"""
Plotting functions for pulse envelopes and population dynamics.

All functions accept an optional ``ax`` argument and return the axes object.
"""

from __future__ import annotations

from typing import Callable, Dict, List, Optional

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.axes import Axes


def plot_pulses(
    tlist: np.ndarray,
    pulse_funcs: List[Callable],
    labels: List[str],
    args: dict,
    ax: Optional[Axes] = None,
) -> Axes:
    """
    Plot pulse envelopes versus time.

    Parameters
    ----------
    tlist : np.ndarray
        Time array.
    pulse_funcs : list of callable
        Pulse functions with signature ``f(t, args) -> float``.
    labels : list of str
        Label for each pulse.
    args : dict
        Arguments forwarded to each pulse function.
    ax : matplotlib Axes, optional
        Axes to draw on. Created if None.

    Returns
    -------
    matplotlib.axes.Axes
    """
    if ax is None:
        _, ax = plt.subplots()

    for func, label in zip(pulse_funcs, labels):
        values = np.array([func(t, args) for t in tlist])
        ax.plot(tlist, values, label=label)

    ax.set_xlabel("Time")
    ax.set_ylabel("Amplitude")
    ax.set_title("Pulse Envelopes")
    ax.legend()
    return ax


def plot_populations(
    result,
    labels: List[str],
    ax: Optional[Axes] = None,
) -> Axes:
    """
    Plot population dynamics from a mesolve result.

    Parameters
    ----------
    result : object
        Result with ``.times`` and ``.expect`` attributes.
        ``result.expect[i]`` is the population of level i over time.
    labels : list of str
        Label for each population curve.
    ax : matplotlib Axes, optional
        Axes to draw on. Created if None.

    Returns
    -------
    matplotlib.axes.Axes
    """
    if ax is None:
        _, ax = plt.subplots()

    times = np.asarray(result.times)
    for i, label in enumerate(labels):
        ax.plot(times, np.real(result.expect[i]), label=label)

    ax.set_xlabel("Time")
    ax.set_ylabel("Population")
    ax.set_title("Population Dynamics")
    ax.legend()
    return ax


def plot_populations_mc(
    result,
    labels: List[str],
    ax: Optional[Axes] = None,
) -> Axes:
    """
    Plot averaged population dynamics from an mcsolve result with ±1σ bands.

    Parameters
    ----------
    result : object
        Result with ``.times``, ``.expect``, and ``.std_expect`` attributes.
    labels : list of str
        Label for each population curve.
    ax : matplotlib Axes, optional
        Axes to draw on. Created if None.

    Returns
    -------
    matplotlib.axes.Axes
    """
    if ax is None:
        _, ax = plt.subplots()

    times = np.asarray(result.times)
    for i, label in enumerate(labels):
        mean = np.real(result.expect[i])
        std = np.real(result.std_expect[i])
        (line,) = ax.plot(times, mean, label=label)
        ax.fill_between(
            times, mean - std, mean + std, alpha=0.2, color=line.get_color()
        )

    ax.set_xlabel("Time")
    ax.set_ylabel("Population")
    ax.set_title("Population Dynamics (MC)")
    ax.legend()
    return ax

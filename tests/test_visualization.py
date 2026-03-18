"""Tests for triqg.visualization -- plotting functions."""

import numpy as np
import pytest
import matplotlib

matplotlib.use("Agg")  # non-interactive backend for CI/testing
import matplotlib.pyplot as plt

from triqg.visualization import plot_pulses, plot_populations, plot_populations_mc


# Simple pulse functions for testing
def _const_pulse(t, args):
    return 1.0


def _zero_pulse(t, args):
    return 0.0


TLIST = np.linspace(0, 1, 20)
ARGS = {}


class TestPlotPulses:
    def test_returns_axes(self):
        ax = plot_pulses(TLIST, [_const_pulse, _zero_pulse], ["p1", "p2"], ARGS)
        assert isinstance(ax, matplotlib.axes.Axes)
        plt.close("all")

    def test_accepts_existing_ax(self):
        fig, ax = plt.subplots()
        returned = plot_pulses(TLIST, [_const_pulse], ["p1"], ARGS, ax=ax)
        assert returned is ax
        plt.close("all")


class _FakeResult:
    """Minimal result-like object for testing population plots."""

    def __init__(self, expect, times, std_expect=None):
        self.expect = expect
        self.times = times
        self.std_expect = std_expect


class TestPlotPopulations:
    def test_returns_axes(self):
        times = np.linspace(0, 1, 10)
        expect = [np.ones(10), np.zeros(10)]
        result = _FakeResult(expect, times)
        ax = plot_populations(result, ["level 0", "level 1"])
        assert isinstance(ax, matplotlib.axes.Axes)
        plt.close("all")

    def test_accepts_existing_ax(self):
        fig, ax = plt.subplots()
        times = np.linspace(0, 1, 10)
        expect = [np.ones(10)]
        result = _FakeResult(expect, times)
        returned = plot_populations(result, ["level 0"], ax=ax)
        assert returned is ax
        plt.close("all")


class TestPlotPopulationsMC:
    def test_returns_axes_with_bands(self):
        times = np.linspace(0, 1, 10)
        expect = [np.linspace(1, 0, 10), np.linspace(0, 1, 10)]
        std = [0.1 * np.ones(10), 0.1 * np.ones(10)]
        result = _FakeResult(expect, times, std_expect=std)
        ax = plot_populations_mc(result, ["P(g)", "P(e)"])
        assert isinstance(ax, matplotlib.axes.Axes)
        plt.close("all")

    def test_accepts_existing_ax(self):
        fig, ax = plt.subplots()
        times = np.linspace(0, 1, 10)
        expect = [np.ones(10)]
        std = [0.05 * np.ones(10)]
        result = _FakeResult(expect, times, std_expect=std)
        returned = plot_populations_mc(result, ["P(g)"], ax=ax)
        assert returned is ax
        plt.close("all")

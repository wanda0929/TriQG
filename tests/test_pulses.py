"""Tests for triqg.pulses -- pulse shapes and area helper."""

import numpy as np
import pytest

from triqg.pulses import omega_c, omega_p, omega_R, compute_pulse_area


# Shared test parameters
ARGS = {
    "omega_c_amp": 2 * np.pi * 50,  # 2pi * 50 MHz
    "omega_p_amp": 2 * np.pi * 50,
    "omega_R_amp": 2.5 * 2 * np.pi * 50,
    "T_c": np.pi / (2 * np.pi * 50),  # pi / Omega_c
    "T_f": 0.5,  # arbitrary for now
    "sigma": 0.05,
}


class TestOmegaC:
    def test_positive_segment(self):
        """omega_c returns +amp/2 in [0, T_c)."""
        amp = ARGS["omega_c_amp"]
        T_c = ARGS["T_c"]
        # sample at the midpoint of the first segment
        t_mid = T_c / 2
        assert omega_c(t_mid, ARGS) == pytest.approx(amp / 2)

    def test_gap_segment_is_zero(self):
        """omega_c returns 0 in [T_c, T_c + T_f)."""
        T_c = ARGS["T_c"]
        T_f = ARGS["T_f"]
        t_gap = T_c + T_f / 2
        assert omega_c(t_gap, ARGS) == pytest.approx(0.0)

    def test_negative_segment(self):
        """omega_c returns -amp/2 in [T_c + T_f, 2*T_c + T_f)."""
        amp = ARGS["omega_c_amp"]
        T_c = ARGS["T_c"]
        T_f = ARGS["T_f"]
        t_neg = T_c + T_f + T_c / 2
        assert omega_c(t_neg, ARGS) == pytest.approx(-amp / 2)

    def test_zero_before_start(self):
        assert omega_c(-0.1, ARGS) == pytest.approx(0.0)

    def test_zero_after_end(self):
        T_c = ARGS["T_c"]
        T_f = ARGS["T_f"]
        assert omega_c(2 * T_c + T_f + 0.1, ARGS) == pytest.approx(0.0)


class TestOmegaP:
    def test_zero_before_window(self):
        """omega_p returns 0 before T_c."""
        assert omega_p(0.0, ARGS) == pytest.approx(0.0)

    def test_zero_after_window(self):
        """omega_p returns 0 after T_c + T_f."""
        T_c = ARGS["T_c"]
        T_f = ARGS["T_f"]
        assert omega_p(T_c + T_f + 0.01, ARGS) == pytest.approx(0.0)

    def test_peak_at_center(self):
        """At center = T_c + T_f/2, the exponent is 0 so value = amp/2."""
        amp = ARGS["omega_p_amp"]
        T_c = ARGS["T_c"]
        T_f = ARGS["T_f"]
        center = T_c + T_f / 2
        assert omega_p(center, ARGS) == pytest.approx(amp / 2)

    def test_symmetric_near_center(self):
        """Pulse is NOT symmetric (cubic exponent), but |value| decreases away from center."""
        T_c = ARGS["T_c"]
        T_f = ARGS["T_f"]
        center = T_c + T_f / 2
        val_center = omega_p(center, ARGS)
        val_offset = omega_p(center + 0.01, ARGS)
        assert abs(val_offset) < abs(val_center)


class TestOmegaR:
    def test_constant_inside_window(self):
        """omega_R returns amp/2 during [T_c, T_c + T_f)."""
        amp = ARGS["omega_R_amp"]
        T_c = ARGS["T_c"]
        T_f = ARGS["T_f"]
        t_mid = T_c + T_f / 2
        assert omega_R(t_mid, ARGS) == pytest.approx(amp / 2)

    def test_zero_before_window(self):
        assert omega_R(0.0, ARGS) == pytest.approx(0.0)

    def test_zero_after_window(self):
        T_c = ARGS["T_c"]
        T_f = ARGS["T_f"]
        assert omega_R(T_c + T_f + 0.01, ARGS) == pytest.approx(0.0)


class TestComputePulseArea:
    def test_constant_pulse_exact_area(self):
        """
        For a constant pulse f(t)=A over [0, T], area = A^2 * T / (2*delta).
        Choose A, T, delta so the analytic result is known.
        """
        A = 10.0
        delta = 5.0
        T = 2.0
        expected_area = A**2 * T / (2 * delta)  # 10^2 * 2 / 10 = 20

        def const_pulse(t, args):
            return A

        area = compute_pulse_area(const_pulse, delta, 0.0, T, args={})
        assert area == pytest.approx(expected_area, rel=1e-6)

    def test_zero_pulse_gives_zero_area(self):
        def zero_pulse(t, args):
            return 0.0

        area = compute_pulse_area(zero_pulse, 1.0, 0.0, 1.0, args={})
        assert area == pytest.approx(0.0)


class TestQuTiPCompatibility:
    def test_omega_c_works_as_qobjevo_coefficient(self):
        """Pulse functions can be used as QuTiP [Qobj, func] coefficients."""
        import qutip

        H0 = qutip.sigmaz()
        H1 = qutip.sigmax()
        # This should not raise -- QobjEvo accepts f(t, args) callables
        H = qutip.QobjEvo([H0, [H1, omega_c]], args=ARGS)
        # Evaluate at a time in the first segment -- should return H0 + (amp/2)*H1
        result = H(ARGS["T_c"] / 2)
        assert result.shape == (2, 2)

    def test_omega_p_works_as_qobjevo_coefficient(self):
        import qutip

        H0 = qutip.sigmaz()
        H1 = qutip.sigmax()
        H = qutip.QobjEvo([H0, [H1, omega_p]], args=ARGS)
        center = ARGS["T_c"] + ARGS["T_f"] / 2
        result = H(center)
        assert result.shape == (2, 2)

"""Tests for triqg.pulses -- pulse shapes and area helper."""

import numpy as np
import pytest

from triqg.pulses import omega_c, omega_p, omega_R, compute_pulse_area
from triqg.pulses import omega_cc, omega_t1, omega_t2


# Shared test parameters
ARGS = {
    "omega_c_amp": 2 * np.pi * 50,  # 2pi * 50 MHz
    "omega_p_amp": 2 * np.pi * 50,
    "omega_R_amp": 2.5 * 2 * np.pi * 50,
    "T_c": np.pi / (2 * np.pi * 50),  # pi / Omega_c
    "T_f": 0.32,  # Hanning T_f for area=pi at delta=2pi*1200
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
        """omega_c returns -amp/2 in [T_c + 2*T_f, 2*T_c + 2*T_f)."""
        amp = ARGS["omega_c_amp"]
        T_c = ARGS["T_c"]
        T_f = ARGS["T_f"]
        t_neg = T_c + 2 * T_f + T_c / 2
        assert omega_c(t_neg, ARGS) == pytest.approx(-amp / 2)

    def test_zero_before_start(self):
        assert omega_c(-0.1, ARGS) == pytest.approx(0.0)

    def test_zero_after_end(self):
        T_c = ARGS["T_c"]
        T_f = ARGS["T_f"]
        assert omega_c(2 * T_c + T_f + 0.1, ARGS) == pytest.approx(0.0)


class TestOmegaP:
    """Tests for the Hanning (sin^2) probe pulse."""

    def test_zero_before_window(self):
        """omega_p returns 0 before T_c."""
        assert omega_p(0.0, ARGS) == pytest.approx(0.0)

    def test_zero_after_window(self):
        """omega_p returns 0 after T_c + 2*T_f."""
        T_c = ARGS["T_c"]
        T_f = ARGS["T_f"]
        assert omega_p(T_c + 2 * T_f + 0.01, ARGS) == pytest.approx(0.0)

    def test_zero_at_start_edge(self):
        """Hanning starts exactly at 0."""
        T_c = ARGS["T_c"]
        assert omega_p(T_c, ARGS) == pytest.approx(0.0)

    def test_zero_at_end_edge(self):
        """Hanning ends exactly at 0."""
        T_c = ARGS["T_c"]
        T_f = ARGS["T_f"]
        assert omega_p(T_c + 2 * T_f, ARGS) == pytest.approx(0.0)

    def test_peak_at_center(self):
        """At center = T_c + T_f, sin^2(pi/2) = 1, so value = amp."""
        amp = ARGS["omega_p_amp"]
        T_c = ARGS["T_c"]
        T_f = ARGS["T_f"]
        center = T_c + T_f
        assert omega_p(center, ARGS) == pytest.approx(amp)

    def test_symmetric(self):
        """Hanning pulse is symmetric about the center."""
        T_c = ARGS["T_c"]
        T_f = ARGS["T_f"]
        center = T_c + T_f
        dt = 0.05
        val_left = omega_p(center - dt, ARGS)
        val_right = omega_p(center + dt, ARGS)
        assert val_left == pytest.approx(val_right)

    def test_quarter_point(self):
        """At t = T_c + T_f/2, sin^2(pi/4) = 0.5, so value = amp/2."""
        amp = ARGS["omega_p_amp"]
        T_c = ARGS["T_c"]
        T_f = ARGS["T_f"]
        t_quarter = T_c + T_f / 2
        assert omega_p(t_quarter, ARGS) == pytest.approx(amp / 2)


class TestOmegaR:
    def test_constant_inside_window(self):
        """omega_R returns amp/2 during [0, 2*(T_c + T_f))."""
        amp = ARGS["omega_R_amp"]
        T_c = ARGS["T_c"]
        T_f = ARGS["T_f"]
        t_mid = T_c + T_f / 2
        assert omega_R(t_mid, ARGS) == pytest.approx(amp / 2)

    def test_nonzero_at_start(self):
        """omega_R is on from t=0."""
        amp = ARGS["omega_R_amp"]
        assert omega_R(0.001, ARGS) == pytest.approx(amp / 2)

    def test_zero_after_window(self):
        T_c = ARGS["T_c"]
        T_f = ARGS["T_f"]
        t_after = 2 * (T_c + T_f) + 0.01
        assert omega_R(t_after, ARGS) == pytest.approx(0.0)


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


# ---------------------------------------------------------------------------
# CCX gate pulse parameters
# ---------------------------------------------------------------------------
CCX_ARGS = {
    "omega_cc_amp": 2 * np.pi * 100,  # 2pi * 100 MHz
    "omega_t_amp": 2 * np.pi * 50,  # 2pi * 50 MHz
    "T_cc": np.pi / (2 * np.pi * 100),  # pi / omega_cc_amp = 5 ns
    "T_t": np.pi / (2 * np.pi * 50),  # pi / omega_t_amp = 10 ns
}


class TestOmegaCc:
    def test_positive_segment(self):
        """omega_cc returns +amp/2 in [0, T_cc)."""
        amp = CCX_ARGS["omega_cc_amp"]
        T_cc = CCX_ARGS["T_cc"]
        t_mid = T_cc / 2
        assert omega_cc(t_mid, CCX_ARGS) == pytest.approx(amp / 2)

    def test_gap_is_zero(self):
        """omega_cc returns 0 during the gap [T_cc, T_cc + 3*T_t)."""
        T_cc = CCX_ARGS["T_cc"]
        T_t = CCX_ARGS["T_t"]
        # Sample at midpoint of gap
        t_gap = T_cc + 1.5 * T_t
        assert omega_cc(t_gap, CCX_ARGS) == pytest.approx(0.0)

    def test_negative_segment(self):
        """omega_cc returns -amp/2 in [T_cc + 3*T_t, 2*T_cc + 3*T_t)."""
        amp = CCX_ARGS["omega_cc_amp"]
        T_cc = CCX_ARGS["T_cc"]
        T_t = CCX_ARGS["T_t"]
        t_neg = T_cc + 3 * T_t + T_cc / 2
        assert omega_cc(t_neg, CCX_ARGS) == pytest.approx(-amp / 2)

    def test_zero_before_start(self):
        assert omega_cc(-0.1, CCX_ARGS) == pytest.approx(0.0)

    def test_zero_after_end(self):
        T_cc = CCX_ARGS["T_cc"]
        T_t = CCX_ARGS["T_t"]
        t_after = 2 * T_cc + 3 * T_t + 0.1
        assert omega_cc(t_after, CCX_ARGS) == pytest.approx(0.0)


class TestOmegaT1:
    def test_first_subpulse(self):
        """omega_t1 returns amp/2 in [T_cc, T_cc + T_t)."""
        amp = CCX_ARGS["omega_t_amp"]
        T_cc = CCX_ARGS["T_cc"]
        T_t = CCX_ARGS["T_t"]
        t_mid = T_cc + T_t / 2
        assert omega_t1(t_mid, CCX_ARGS) == pytest.approx(amp / 2)

    def test_third_subpulse(self):
        """omega_t1 returns amp/2 in [T_cc + 2*T_t, T_cc + 3*T_t)."""
        amp = CCX_ARGS["omega_t_amp"]
        T_cc = CCX_ARGS["T_cc"]
        T_t = CCX_ARGS["T_t"]
        t_mid = T_cc + 2.5 * T_t
        assert omega_t1(t_mid, CCX_ARGS) == pytest.approx(amp / 2)

    def test_zero_in_second_subpulse_gap(self):
        """omega_t1 returns 0 during the second sub-pulse window [T_cc + T_t, T_cc + 2*T_t)."""
        T_cc = CCX_ARGS["T_cc"]
        T_t = CCX_ARGS["T_t"]
        t_gap = T_cc + 1.5 * T_t
        assert omega_t1(t_gap, CCX_ARGS) == pytest.approx(0.0)

    def test_zero_before_start(self):
        assert omega_t1(0.0, CCX_ARGS) == pytest.approx(0.0)

    def test_zero_after_end(self):
        T_cc = CCX_ARGS["T_cc"]
        T_t = CCX_ARGS["T_t"]
        assert omega_t1(T_cc + 3 * T_t + 0.1, CCX_ARGS) == pytest.approx(0.0)


class TestOmegaT2:
    def test_second_subpulse(self):
        """omega_t2 returns amp/2 in [T_cc + T_t, T_cc + 2*T_t)."""
        amp = CCX_ARGS["omega_t_amp"]
        T_cc = CCX_ARGS["T_cc"]
        T_t = CCX_ARGS["T_t"]
        t_mid = T_cc + 1.5 * T_t
        assert omega_t2(t_mid, CCX_ARGS) == pytest.approx(amp / 2)

    def test_zero_in_first_subpulse_window(self):
        """omega_t2 returns 0 during the first sub-pulse window."""
        T_cc = CCX_ARGS["T_cc"]
        T_t = CCX_ARGS["T_t"]
        t_mid = T_cc + T_t / 2
        assert omega_t2(t_mid, CCX_ARGS) == pytest.approx(0.0)

    def test_zero_in_third_subpulse_window(self):
        """omega_t2 returns 0 during the third sub-pulse window."""
        T_cc = CCX_ARGS["T_cc"]
        T_t = CCX_ARGS["T_t"]
        t_mid = T_cc + 2.5 * T_t
        assert omega_t2(t_mid, CCX_ARGS) == pytest.approx(0.0)

    def test_zero_before_start(self):
        assert omega_t2(0.0, CCX_ARGS) == pytest.approx(0.0)

    def test_zero_after_end(self):
        T_cc = CCX_ARGS["T_cc"]
        T_t = CCX_ARGS["T_t"]
        assert omega_t2(T_cc + 2 * T_t + 0.1, CCX_ARGS) == pytest.approx(0.0)


class TestHanningPulseArea:
    """
    Verify compute_pulse_area of the Hanning (sin^2) probe pulse.

    For the Hanning window, the analytical area is:
        area = omega_p_amp^2 * 2*T_f * (3/8) / (2*delta)

    With delta/2pi = 1200 MHz, omega_p_amp/2pi = 50 MHz, T_f = 0.32 us:
        area = pi
    """

    DELTA = 2 * np.pi * 1200

    def _area(self, T_f):
        a = {
            "omega_c_amp": 2 * np.pi * 50,
            "omega_p_amp": 2 * np.pi * 50,
            "omega_R_amp": 0,
            "T_c": np.pi / (2 * np.pi * 50),
            "T_f": T_f,
        }
        T_c = a["T_c"]
        return compute_pulse_area(omega_p, self.DELTA, T_c, T_c + 2 * T_f, a)

    def test_area_increases_with_Tf(self):
        """Larger T_f -> wider pulse -> more area."""
        Tfs = [0.10, 0.20, 0.32, 0.50]
        areas = [self._area(tf) for tf in Tfs]
        for i in range(len(areas) - 1):
            assert areas[i] < areas[i + 1]

    def test_area_matches_analytical(self):
        """Numerical area matches the analytical 3/8 formula."""
        omega_p_amp = 2 * np.pi * 50
        T_f = 0.32
        expected = omega_p_amp**2 * (2 * T_f) * (3 / 8) / (2 * self.DELTA)
        area = self._area(T_f)
        assert area == pytest.approx(expected, rel=1e-4)

    def test_area_equals_pi_at_design_Tf(self):
        """T_f = 0.32 us gives area = pi with delta/2pi = 1200 MHz."""
        area = self._area(0.32)
        assert area == pytest.approx(np.pi, rel=0.01)

    def test_area_scales_linearly_with_Tf(self):
        """Hanning area is proportional to T_f (since shape is self-similar)."""
        a1 = self._area(0.20)
        a2 = self._area(0.40)
        assert a2 == pytest.approx(2 * a1, rel=1e-3)


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

    def test_omega_cc_works_as_qobjevo_coefficient(self):
        import qutip

        H0 = qutip.sigmaz()
        H1 = qutip.sigmax()
        H = qutip.QobjEvo([H0, [H1, omega_cc]], args=CCX_ARGS)
        result = H(CCX_ARGS["T_cc"] / 2)
        assert result.shape == (2, 2)

    def test_omega_t1_works_as_qobjevo_coefficient(self):
        import qutip

        H0 = qutip.sigmaz()
        H1 = qutip.sigmax()
        H = qutip.QobjEvo([H0, [H1, omega_t1]], args=CCX_ARGS)
        T_cc = CCX_ARGS["T_cc"]
        T_t = CCX_ARGS["T_t"]
        result = H(T_cc + T_t / 2)
        assert result.shape == (2, 2)

    def test_omega_t2_works_as_qobjevo_coefficient(self):
        import qutip

        H0 = qutip.sigmaz()
        H1 = qutip.sigmax()
        H = qutip.QobjEvo([H0, [H1, omega_t2]], args=CCX_ARGS)
        T_cc = CCX_ARGS["T_cc"]
        T_t = CCX_ARGS["T_t"]
        result = H(T_cc + 1.5 * T_t)
        assert result.shape == (2, 2)

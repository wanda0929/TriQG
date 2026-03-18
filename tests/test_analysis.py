"""Tests for triqg.analysis -- fidelity and population extraction."""

import numpy as np
import pytest
import qutip

from triqg.analysis import state_fidelity, extract_populations


class TestStateFidelity:
    def test_ket_fidelity_with_itself_is_one(self):
        psi = qutip.basis(2, 0)
        assert state_fidelity(psi, psi) == pytest.approx(1.0)

    def test_ket_fidelity_orthogonal_is_zero(self):
        psi0 = qutip.basis(2, 0)
        psi1 = qutip.basis(2, 1)
        assert state_fidelity(psi0, psi1) == pytest.approx(0.0, abs=1e-10)

    def test_density_matrix_fidelity_with_itself_is_one(self):
        rho = qutip.ket2dm(qutip.basis(3, 1))
        assert state_fidelity(rho, rho) == pytest.approx(1.0)

    def test_density_matrix_fidelity_orthogonal_is_zero(self):
        rho0 = qutip.ket2dm(qutip.basis(3, 0))
        rho1 = qutip.ket2dm(qutip.basis(3, 2))
        assert state_fidelity(rho0, rho1) == pytest.approx(0.0, abs=1e-10)

    def test_mixed_ket_and_dm_inputs(self):
        """Fidelity works when one input is a ket and the other is a dm."""
        psi = qutip.basis(2, 0)
        rho = qutip.ket2dm(psi)
        assert state_fidelity(psi, rho) == pytest.approx(1.0)
        assert state_fidelity(rho, psi) == pytest.approx(1.0)

    def test_superposition_fidelity(self):
        """Fidelity of |0> with (|0>+|1>)/sqrt(2) should be 0.5."""
        psi0 = qutip.basis(2, 0)
        psi_plus = (qutip.basis(2, 0) + qutip.basis(2, 1)).unit()
        assert state_fidelity(psi0, psi_plus) == pytest.approx(0.5, abs=1e-10)


class TestExtractPopulations:
    def test_returns_correct_shape(self):
        """
        Given a result with 3 expect arrays of length 10,
        extract_populations(result, 3) should return shape (3, 10).
        """

        # Simulate a result-like object with .expect attribute
        class FakeResult:
            expect = [np.ones(10), np.zeros(10), 0.5 * np.ones(10)]

        pops = extract_populations(FakeResult(), 3)
        assert pops.shape == (3, 10)

    def test_returns_correct_values(self):
        """Values should match the input expect arrays."""
        arr0 = np.array([1.0, 0.8, 0.5, 0.2, 0.0])
        arr1 = np.array([0.0, 0.2, 0.5, 0.8, 1.0])

        class FakeResult:
            expect = [arr0, arr1]

        pops = extract_populations(FakeResult(), 2)
        np.testing.assert_array_almost_equal(pops[0], arr0)
        np.testing.assert_array_almost_equal(pops[1], arr1)

    def test_takes_real_part(self):
        """Populations should be real even if expect has small imaginary parts."""
        arr = np.array([1.0 + 1e-15j, 0.5 + 1e-14j])

        class FakeResult:
            expect = [arr]

        pops = extract_populations(FakeResult(), 1)
        assert pops.dtype == np.float64

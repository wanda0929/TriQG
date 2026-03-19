"""Tests for triqg.analysis -- fidelity and population extraction."""

import numpy as np
import pytest
import qutip

from triqg.analysis import state_fidelity, extract_populations, average_gate_fidelity


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


class TestAverageGateFidelity:
    def test_all_perfect_pairs_returns_one(self):
        """When every final state matches its target, average fidelity is 1.0."""
        psi0 = qutip.basis(2, 0)
        psi1 = qutip.basis(2, 1)
        pairs = [(psi0, psi0), (psi1, psi1)]
        assert average_gate_fidelity(pairs) == pytest.approx(1.0)

    def test_mixed_fidelity_returns_correct_mean(self):
        """One perfect pair (F=1) and one orthogonal pair (F=0) averages to 0.5."""
        psi0 = qutip.basis(2, 0)
        psi1 = qutip.basis(2, 1)
        pairs = [(psi0, psi0), (psi0, psi1)]
        assert average_gate_fidelity(pairs) == pytest.approx(0.5, abs=1e-10)

    def test_works_with_density_matrices(self):
        """Average fidelity accepts density matrix pairs."""
        rho0 = qutip.ket2dm(qutip.basis(2, 0))
        rho1 = qutip.ket2dm(qutip.basis(2, 1))
        pairs = [(rho0, rho0), (rho1, rho1)]
        assert average_gate_fidelity(pairs) == pytest.approx(1.0)

    def test_single_pair(self):
        """Works correctly with a single (final, expected) pair."""
        psi0 = qutip.basis(2, 0)
        psi_plus = (qutip.basis(2, 0) + qutip.basis(2, 1)).unit()
        # F(|0>, |+>) = 0.5, so average over one pair is 0.5
        pairs = [(psi0, psi_plus)]
        assert average_gate_fidelity(pairs) == pytest.approx(0.5, abs=1e-10)

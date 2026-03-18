"""Tests for triqg.solver -- unified mesolve/mcsolve interface."""

import numpy as np
import pytest
import qutip

from triqg.solver import simulate, SimulationResult


# -- Trivial 2-level decay test system --
# |e> decays to |g> with rate gamma=1.0. No drive.
H_trivial = qutip.sigmaz()  # static H (just energy splitting)
psi0_excited = qutip.basis(2, 1)  # start in |e>
gamma = 1.0
c_ops_trivial = [np.sqrt(gamma) * qutip.destroy(2)]  # |e> -> |g>
e_ops_trivial = [
    qutip.ket2dm(qutip.basis(2, 0)),  # P(|g>)
    qutip.ket2dm(qutip.basis(2, 1)),
]  # P(|e>)
tlist = np.linspace(0, 3, 50)


class TestMesolve:
    def test_returns_simulation_result_with_times(self):
        result = simulate(
            method="mesolve",
            H=H_trivial,
            psi0=psi0_excited,
            tlist=tlist,
            c_ops=c_ops_trivial,
            e_ops=e_ops_trivial,
        )
        assert isinstance(result, SimulationResult)
        np.testing.assert_array_almost_equal(result.times, tlist)

    def test_expect_has_correct_length_and_decays(self):
        """Excited state population should decay; ground should increase."""
        result = simulate(
            method="mesolve",
            H=H_trivial,
            psi0=psi0_excited,
            tlist=tlist,
            c_ops=c_ops_trivial,
            e_ops=e_ops_trivial,
        )
        assert len(result.expect) == 2
        assert len(result.expect[0]) == len(tlist)
        # P(|e>) should decrease from ~1 to ~0
        pe = np.real(result.expect[1])
        assert pe[0] == pytest.approx(1.0, abs=0.05)
        assert pe[-1] < 0.1

    def test_final_state_is_qobj(self):
        result = simulate(
            method="mesolve",
            H=H_trivial,
            psi0=psi0_excited,
            tlist=tlist,
            c_ops=c_ops_trivial,
            e_ops=e_ops_trivial,
            options={"store_final_state": True},
        )
        assert isinstance(result.final_state, qutip.Qobj)

    def test_mc_attributes_are_none(self):
        """mesolve result should return None for MC-specific attributes."""
        result = simulate(
            method="mesolve",
            H=H_trivial,
            psi0=psi0_excited,
            tlist=tlist,
            c_ops=c_ops_trivial,
            e_ops=e_ops_trivial,
        )
        assert result.std_expect is None
        assert result.num_trajectories is None
        assert result.col_times is None
        assert result.col_which is None


class TestMcsolve:
    def test_returns_simulation_result_with_times(self):
        result = simulate(
            method="mcsolve",
            H=H_trivial,
            psi0=psi0_excited,
            tlist=tlist,
            c_ops=c_ops_trivial,
            e_ops=e_ops_trivial,
            ntraj=20,
        )
        assert isinstance(result, SimulationResult)
        np.testing.assert_array_almost_equal(result.times, tlist)

    def test_expect_returns_averaged_values(self):
        """MC expect should be averaged over trajectories and decay like mesolve."""
        result = simulate(
            method="mcsolve",
            H=H_trivial,
            psi0=psi0_excited,
            tlist=tlist,
            c_ops=c_ops_trivial,
            e_ops=e_ops_trivial,
            ntraj=50,
        )
        assert len(result.expect) == 2
        pe = np.real(result.expect[1])  # P(|e>)
        assert pe[0] == pytest.approx(1.0, abs=0.1)
        assert pe[-1] < 0.2  # should have mostly decayed

    def test_mc_specific_attributes(self):
        result = simulate(
            method="mcsolve",
            H=H_trivial,
            psi0=psi0_excited,
            tlist=tlist,
            c_ops=c_ops_trivial,
            e_ops=e_ops_trivial,
            ntraj=10,
        )
        # std_expect
        assert result.std_expect is not None
        assert len(result.std_expect) == 2

        # num_trajectories
        assert result.num_trajectories == 10

        # col_times and col_which are lists (one per trajectory)
        assert result.col_times is not None
        assert result.col_which is not None
        assert len(result.col_times) == 10

    def test_final_state_is_qobj(self):
        result = simulate(
            method="mcsolve",
            H=H_trivial,
            psi0=psi0_excited,
            tlist=tlist,
            c_ops=c_ops_trivial,
            e_ops=e_ops_trivial,
            ntraj=10,
            options={"store_final_state": True},
        )
        fs = result.final_state
        assert isinstance(fs, qutip.Qobj)


class TestInvalidMethod:
    def test_raises_on_unknown_method(self):
        with pytest.raises(ValueError, match="Unknown method"):
            simulate(
                method="unknown",
                H=H_trivial,
                psi0=psi0_excited,
                tlist=tlist,
            )

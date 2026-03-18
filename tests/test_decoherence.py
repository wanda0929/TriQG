"""Tests for triqg.decoherence -- collapse operators for decay channels."""

import numpy as np
import pytest
import qutip

from triqg.decoherence import build_collapse_operators
from triqg.atoms import composite_basis_state, CsAtom, RbAtom


# Physical decay rates (in MHz, consistent with pulse units)
GAMMA_r = 1.0 / 548.0  # Cs Rydberg |r> lifetime 548 us
GAMMA_R = 1.0 / 505.0  # Rb Rydberg |R> lifetime 505 us
GAMMA_P = 1.0 / 0.131  # Rb intermediate |P> lifetime 0.131 us


class TestBuildCollapseOperators:
    def test_returns_list_of_eight_operators(self):
        c_ops = build_collapse_operators(GAMMA_r, GAMMA_R, GAMMA_P)
        assert isinstance(c_ops, list)
        assert len(c_ops) == 8

    def test_all_operators_are_36x36(self):
        c_ops = build_collapse_operators(GAMMA_r, GAMMA_R, GAMMA_P)
        for i, op in enumerate(c_ops):
            assert isinstance(op, qutip.Qobj), f"c_ops[{i}] not Qobj"
            assert op.shape == (36, 36), f"c_ops[{i}] wrong shape"

    def test_default_branching_rate_scaling(self):
        """
        With default 0.5/0.5 branching, applying c_ops[0] to a state where
        control 1 is in |r> should produce a state with norm sqrt(gamma_r * 0.5).

        c_ops[0] |r, 0, A> = sqrt(gamma_r * 0.5) * |0, 0, A>
        """
        c_ops = build_collapse_operators(GAMMA_r, GAMMA_R, GAMMA_P)
        psi_r = composite_basis_state(2, 0, 0)  # |r, 0, A>
        result = c_ops[0] * psi_r
        expected_norm = np.sqrt(GAMMA_r * 0.5)
        assert result.norm() == pytest.approx(expected_norm, rel=1e-6)


class TestCollapseTransitions:
    """Verify each operator maps the correct |initial> to the correct |final>."""

    cs = CsAtom()
    rb = RbAtom()
    IDX_0 = cs.level_index["0"]
    IDX_1 = cs.level_index["1"]
    IDX_r = cs.level_index["r"]
    IDX_A = rb.level_index["A"]
    IDX_B = rb.level_index["B"]
    IDX_P = rb.level_index["P"]
    IDX_R = rb.level_index["R"]

    def _applied_direction(self, c_op, initial_state, expected_final):
        """
        Check that c_op * |initial> is proportional to |final>,
        i.e. the overlap |<final|result>| / ||result|| == 1.
        """
        result = c_op * initial_state
        if result.norm() < 1e-15:
            return False  # operator annihilates this state
        result_normed = result.unit()
        overlap = abs(expected_final.dag() * result_normed)
        return float(overlap) > 0.999

    def test_c1_r_to_0(self):
        """c_ops[0]: control 1 |r> -> |0>."""
        c_ops = build_collapse_operators(GAMMA_r, GAMMA_R, GAMMA_P)
        initial = composite_basis_state(self.IDX_r, self.IDX_0, self.IDX_A)
        final = composite_basis_state(self.IDX_0, self.IDX_0, self.IDX_A)
        assert self._applied_direction(c_ops[0], initial, final)

    def test_c1_r_to_1(self):
        """c_ops[1]: control 1 |r> -> |1>."""
        c_ops = build_collapse_operators(GAMMA_r, GAMMA_R, GAMMA_P)
        initial = composite_basis_state(self.IDX_r, self.IDX_0, self.IDX_A)
        final = composite_basis_state(self.IDX_1, self.IDX_0, self.IDX_A)
        assert self._applied_direction(c_ops[1], initial, final)

    def test_c2_r_to_0(self):
        """c_ops[2]: control 2 |r> -> |0>."""
        c_ops = build_collapse_operators(GAMMA_r, GAMMA_R, GAMMA_P)
        initial = composite_basis_state(self.IDX_0, self.IDX_r, self.IDX_A)
        final = composite_basis_state(self.IDX_0, self.IDX_0, self.IDX_A)
        assert self._applied_direction(c_ops[2], initial, final)

    def test_c2_r_to_1(self):
        """c_ops[3]: control 2 |r> -> |1>."""
        c_ops = build_collapse_operators(GAMMA_r, GAMMA_R, GAMMA_P)
        initial = composite_basis_state(self.IDX_0, self.IDX_r, self.IDX_A)
        final = composite_basis_state(self.IDX_0, self.IDX_1, self.IDX_A)
        assert self._applied_direction(c_ops[3], initial, final)

    def test_target_R_to_A(self):
        """c_ops[4]: target |R> -> |A>."""
        c_ops = build_collapse_operators(GAMMA_r, GAMMA_R, GAMMA_P)
        initial = composite_basis_state(self.IDX_0, self.IDX_0, self.IDX_R)
        final = composite_basis_state(self.IDX_0, self.IDX_0, self.IDX_A)
        assert self._applied_direction(c_ops[4], initial, final)

    def test_target_R_to_B(self):
        """c_ops[5]: target |R> -> |B>."""
        c_ops = build_collapse_operators(GAMMA_r, GAMMA_R, GAMMA_P)
        initial = composite_basis_state(self.IDX_0, self.IDX_0, self.IDX_R)
        final = composite_basis_state(self.IDX_0, self.IDX_0, self.IDX_B)
        assert self._applied_direction(c_ops[5], initial, final)

    def test_target_P_to_A(self):
        """c_ops[6]: target |P> -> |A>."""
        c_ops = build_collapse_operators(GAMMA_r, GAMMA_R, GAMMA_P)
        initial = composite_basis_state(self.IDX_0, self.IDX_0, self.IDX_P)
        final = composite_basis_state(self.IDX_0, self.IDX_0, self.IDX_A)
        assert self._applied_direction(c_ops[6], initial, final)

    def test_target_P_to_B(self):
        """c_ops[7]: target |P> -> |B>."""
        c_ops = build_collapse_operators(GAMMA_r, GAMMA_R, GAMMA_P)
        initial = composite_basis_state(self.IDX_0, self.IDX_0, self.IDX_P)
        final = composite_basis_state(self.IDX_0, self.IDX_0, self.IDX_B)
        assert self._applied_direction(c_ops[7], initial, final)

    def test_operator_annihilates_wrong_initial(self):
        """c_ops[0] (|r>->|0> on c1) applied to |0,0,A> should give zero."""
        c_ops = build_collapse_operators(GAMMA_r, GAMMA_R, GAMMA_P)
        wrong_state = composite_basis_state(self.IDX_0, self.IDX_0, self.IDX_A)
        result = c_ops[0] * wrong_state
        assert result.norm() < 1e-15


class TestCustomBranching:
    def test_asymmetric_branching_changes_rate(self):
        """
        With branching r=(0.8, 0.2), the |r>->|0> channel (c_ops[0])
        should have norm sqrt(gamma_r * 0.8) when applied to |r,0,A>,
        and |r>->|1> (c_ops[1]) should have norm sqrt(gamma_r * 0.2).
        """
        custom = {"r": (0.8, 0.2), "R": (0.5, 0.5), "P": (0.5, 0.5)}
        c_ops = build_collapse_operators(GAMMA_r, GAMMA_R, GAMMA_P, branching=custom)

        psi = composite_basis_state(2, 0, 0)  # |r, 0, A>

        result_0 = c_ops[0] * psi  # |r> -> |0>
        assert result_0.norm() == pytest.approx(np.sqrt(GAMMA_r * 0.8), rel=1e-6)

        result_1 = c_ops[1] * psi  # |r> -> |1>
        assert result_1.norm() == pytest.approx(np.sqrt(GAMMA_r * 0.2), rel=1e-6)


class TestZeroRate:
    def test_zero_gamma_gives_zero_operators(self):
        """When gamma_r=0, control collapse operators should be zero."""
        c_ops = build_collapse_operators(0.0, GAMMA_R, GAMMA_P)
        # First 4 operators are control channels
        for i in range(4):
            assert c_ops[i].norm() < 1e-15, f"c_ops[{i}] should be zero"
        # Target operators should still be nonzero
        for i in range(4, 8):
            assert c_ops[i].norm() > 1e-15, f"c_ops[{i}] should be nonzero"

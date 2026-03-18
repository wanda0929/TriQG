"""Tests for triqg.atoms -- atomic level structures and composite Hilbert space."""

import pytest
from triqg.atoms import CsAtom, RbAtom, composite_basis_state, composite_projector, DIMS


class TestCsAtom:
    def test_has_three_levels(self):
        atom = CsAtom()
        assert atom.num_levels == 3

    def test_level_labels(self):
        atom = CsAtom()
        assert atom.labels == ["0", "1", "r"]

    def test_level_index_lookup(self):
        atom = CsAtom()
        assert atom.level_index["0"] == 0
        assert atom.level_index["1"] == 1
        assert atom.level_index["r"] == 2


class TestRbAtom:
    def test_has_four_levels(self):
        atom = RbAtom()
        assert atom.num_levels == 4

    def test_level_labels(self):
        atom = RbAtom()
        assert atom.labels == ["A", "B", "P", "R"]

    def test_level_index_lookup(self):
        atom = RbAtom()
        assert atom.level_index["A"] == 0
        assert atom.level_index["B"] == 1
        assert atom.level_index["P"] == 2
        assert atom.level_index["R"] == 3


class TestCompositeBasisState:
    def test_total_dimension(self):
        """3 x 3 x 4 = 36 dimensional Hilbert space."""
        assert DIMS == [3, 3, 4]

    def test_basis_state_shape(self):
        psi = composite_basis_state(0, 0, 0)
        assert psi.shape == (36, 1)

    def test_basis_state_is_ket(self):
        psi = composite_basis_state(0, 0, 0)
        assert psi.isket

    def test_basis_state_is_normalized(self):
        psi = composite_basis_state(1, 2, 3)
        assert psi.norm() == pytest.approx(1.0)

    def test_all_36_basis_states_are_orthonormal(self):
        """Every pair of distinct basis states has zero overlap."""
        states = []
        for c1 in range(DIMS[0]):
            for c2 in range(DIMS[1]):
                for t in range(DIMS[2]):
                    states.append(composite_basis_state(c1, c2, t))
        assert len(states) == 36
        for i, si in enumerate(states):
            for j, sj in enumerate(states):
                overlap = abs(si.overlap(sj))
                if i == j:
                    assert overlap == pytest.approx(1.0)
                else:
                    assert overlap == pytest.approx(0.0)


class TestCompositeProjector:
    def test_projector_shape(self):
        """Projector for control1, level 0 should be 36x36."""
        P = composite_projector(subsystem=0, level=0)
        assert P.shape == (36, 36)

    def test_projector_is_operator(self):
        P = composite_projector(subsystem=0, level=0)
        assert P.isoper

    def test_projector_trace_is_one_per_degeneracy(self):
        """Tr(|0><0|_c1 ⊗ I_c2 ⊗ I_t) = 1 * 3 * 4 = 12."""
        P = composite_projector(subsystem=0, level=0)
        assert float(P.tr().real) == pytest.approx(12.0)

    def test_projector_is_idempotent(self):
        P = composite_projector(subsystem=2, level=3)
        diff = P * P - P
        assert diff.norm() < 1e-12

    def test_projector_is_hermitian(self):
        P = composite_projector(subsystem=1, level=2)
        diff = P - P.dag()
        assert diff.norm() < 1e-12

    def test_subsystem_projectors_sum_to_identity(self):
        """Sum of projectors over all levels of one subsystem = I_36."""
        import qutip

        identity_36 = qutip.tensor(
            qutip.qeye(DIMS[0]),
            qutip.qeye(DIMS[1]),
            qutip.qeye(DIMS[2]),
        )
        for sub in range(3):
            total = sum(composite_projector(sub, lev) for lev in range(DIMS[sub]))
            diff = total - identity_36
            assert diff.norm() < 1e-12, f"subsystem {sub} projectors don't sum to I"

"""Tests for triqg.hamiltonian -- Hamiltonian construction for the OR gate."""

import numpy as np
import pytest
import qutip

from triqg.hamiltonian import build_hamiltonian
from triqg.atoms import CsAtom, RbAtom, DIMS, composite_basis_state


# Reference parameters
DELTA = 2 * np.pi * 1200  # detuning
V_CT = 2 * np.pi * 100  # Rydberg interaction strength

# Level indices for readability
cs = CsAtom()
rb = RbAtom()
IDX_0 = cs.level_index["0"]  # 0
IDX_1 = cs.level_index["1"]  # 1
IDX_r = cs.level_index["r"]  # 2
IDX_A = rb.level_index["A"]  # 0
IDX_B = rb.level_index["B"]  # 1
IDX_P = rb.level_index["P"]  # 2
IDX_R = rb.level_index["R"]  # 3


class TestBuildHamiltonianStructure:
    def test_returns_list_with_static_qobj_first(self):
        H = build_hamiltonian(DELTA, V_CT)
        assert isinstance(H, list)
        assert len(H) > 0
        H_static = H[0]
        assert isinstance(H_static, qutip.Qobj)
        assert H_static.shape == (36, 36)


class TestStaticHamiltonian:
    def test_detuning_on_P_levels(self):
        """
        Every composite state with target in |P> should have diagonal
        element = DELTA (plus possible V_ct if control is also in |r>).
        Check a state with no control in |r> so only DELTA contributes.
        """
        H = build_hamiltonian(DELTA, V_CT)
        H_static = H[0]
        # |0, 0, P> -- controls in |0>, target in |P>
        psi = composite_basis_state(IDX_0, IDX_0, IDX_P)
        energy = qutip.expect(H_static, psi)
        assert energy == pytest.approx(DELTA)

    def test_rydberg_shift_control1_in_r_target_in_R(self):
        """
        |r, 0, R> should have diagonal = V_ct (no DELTA since target not in |P>).
        """
        H = build_hamiltonian(DELTA, V_CT)
        H_static = H[0]
        psi = composite_basis_state(IDX_r, IDX_0, IDX_R)
        energy = qutip.expect(H_static, psi)
        assert energy == pytest.approx(V_CT)

    def test_rydberg_shift_control2_in_r_target_in_R(self):
        """
        |0, r, R> should have diagonal = V_ct.
        """
        H = build_hamiltonian(DELTA, V_CT)
        H_static = H[0]
        psi = composite_basis_state(IDX_0, IDX_r, IDX_R)
        energy = qutip.expect(H_static, psi)
        assert energy == pytest.approx(V_CT)

    def test_both_controls_in_r_target_in_R(self):
        """
        |r, r, R> should have diagonal = 2 * V_ct (both controls contribute).
        """
        H = build_hamiltonian(DELTA, V_CT)
        H_static = H[0]
        psi = composite_basis_state(IDX_r, IDX_r, IDX_R)
        energy = qutip.expect(H_static, psi)
        assert energy == pytest.approx(2 * V_CT)

    def test_no_shift_when_target_not_in_R(self):
        """
        |r, 0, A> should have zero energy (no P detuning, no R interaction).
        """
        H = build_hamiltonian(DELTA, V_CT)
        H_static = H[0]
        psi = composite_basis_state(IDX_r, IDX_0, IDX_A)
        energy = qutip.expect(H_static, psi)
        assert energy == pytest.approx(0.0)

    def test_static_hamiltonian_is_hermitian(self):
        H = build_hamiltonian(DELTA, V_CT)
        H_static = H[0]
        assert (H_static - H_static.dag()).norm() < 1e-12


class TestDriveOperators:
    """Verify each drive operator couples exactly the correct transitions."""

    def _get_drive_op(self, index):
        """Return the operator part of the index-th drive term."""
        H = build_hamiltonian(DELTA, V_CT)
        # H[0] is static; H[1]..H[4] are [op, coeff] pairs
        return H[index][0]

    def _matrix_element(self, op, bra_state, ket_state):
        """Compute <bra|op|ket> as a complex number."""
        return complex(bra_state.dag() * op * ket_state)

    def test_control1_drive_couples_1_to_r(self):
        """
        H_drive_c1 should couple |1,*,*> <-> |r,*,*>.
        Check: <r,0,A| H |1,0,A> is nonzero.
        """
        H_c1 = self._get_drive_op(1)
        bra = composite_basis_state(IDX_r, IDX_0, IDX_A)
        ket = composite_basis_state(IDX_1, IDX_0, IDX_A)
        elem = self._matrix_element(H_c1, bra, ket)
        assert abs(elem) > 0.5  # should be 1.0

    def test_control1_drive_does_not_couple_0_to_r(self):
        """H_drive_c1 should NOT couple |0> <-> |r>."""
        H_c1 = self._get_drive_op(1)
        bra = composite_basis_state(IDX_r, IDX_0, IDX_A)
        ket = composite_basis_state(IDX_0, IDX_0, IDX_A)
        elem = self._matrix_element(H_c1, bra, ket)
        assert abs(elem) < 1e-12

    def test_control2_drive_couples_1_to_r(self):
        H_c2 = self._get_drive_op(2)
        bra = composite_basis_state(IDX_0, IDX_r, IDX_A)
        ket = composite_basis_state(IDX_0, IDX_1, IDX_A)
        elem = self._matrix_element(H_c2, bra, ket)
        assert abs(elem) > 0.5

    def test_probe_drive_couples_A_to_P(self):
        """H_drive_p should couple |*,*,A> <-> |*,*,P>."""
        H_p = self._get_drive_op(3)
        bra = composite_basis_state(IDX_0, IDX_0, IDX_P)
        ket = composite_basis_state(IDX_0, IDX_0, IDX_A)
        elem = self._matrix_element(H_p, bra, ket)
        assert abs(elem) > 0.5

    def test_probe_drive_couples_B_to_P(self):
        """H_drive_p should also couple |*,*,B> <-> |*,*,P>."""
        H_p = self._get_drive_op(3)
        bra = composite_basis_state(IDX_0, IDX_0, IDX_P)
        ket = composite_basis_state(IDX_0, IDX_0, IDX_B)
        elem = self._matrix_element(H_p, bra, ket)
        assert abs(elem) > 0.5

    def test_probe_drive_does_not_couple_A_to_R(self):
        """H_drive_p should NOT couple |A> to |R> directly."""
        H_p = self._get_drive_op(3)
        bra = composite_basis_state(IDX_0, IDX_0, IDX_R)
        ket = composite_basis_state(IDX_0, IDX_0, IDX_A)
        elem = self._matrix_element(H_p, bra, ket)
        assert abs(elem) < 1e-12

    def test_rydberg_coupling_drives_P_to_R(self):
        """H_drive_R should couple |*,*,P> <-> |*,*,R>."""
        H_R = self._get_drive_op(4)
        bra = composite_basis_state(IDX_0, IDX_0, IDX_R)
        ket = composite_basis_state(IDX_0, IDX_0, IDX_P)
        elem = self._matrix_element(H_R, bra, ket)
        assert abs(elem) > 0.5

    def test_rydberg_coupling_does_not_drive_A_to_R(self):
        """H_drive_R should NOT couple |A> to |R> directly."""
        H_R = self._get_drive_op(4)
        bra = composite_basis_state(IDX_0, IDX_0, IDX_R)
        ket = composite_basis_state(IDX_0, IDX_0, IDX_A)
        elem = self._matrix_element(H_R, bra, ket)
        assert abs(elem) < 1e-12

    def test_all_drive_operators_are_hermitian(self):
        H = build_hamiltonian(DELTA, V_CT)
        for i in range(1, len(H)):
            op = H[i][0]
            assert (op - op.dag()).norm() < 1e-12, f"Drive {i} is not Hermitian"

    def test_all_operators_are_36x36(self):
        """Every operator in the Hamiltonian list has shape (36, 36)."""
        H = build_hamiltonian(DELTA, V_CT)
        # H[0] is static Qobj; H[1]..H[4] are [Qobj, func] pairs
        assert H[0].shape == (36, 36)
        for i in range(1, len(H)):
            assert H[i][0].shape == (36, 36), f"Drive {i} wrong shape"

    def test_list_has_five_elements(self):
        """1 static + 4 drives = 5 elements."""
        H = build_hamiltonian(DELTA, V_CT)
        assert len(H) == 5

    def test_drive_coefficients_are_callable(self):
        """Each drive's coefficient is a callable (pulse function)."""
        H = build_hamiltonian(DELTA, V_CT)
        for i in range(1, len(H)):
            assert callable(H[i][1]), f"Drive {i} coeff not callable"

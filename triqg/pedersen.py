"""
Pedersen Haar-averaged gate fidelity for subspace-encoded quantum gates.

Implements the average fidelity formula of

    L. H. Pedersen, N. M. Moller, K. Molmer,
    "Fidelity of quantum operations",
    Phys. Lett. A 367, 47 (2007).

For a completely positive trace-preserving (CPTP) map E acting on the full
simulator Hilbert space (here 36-dim: Cs x Cs x Rb = 3 x 3 x 4), with a
target unitary U0 defined on the d-dimensional computational subspace S
(here d = 2^3 = 8), Pedersen Eqs. (3) and (5) give the Haar-averaged
fidelity over pure input states drawn uniformly from S:

    F_bar = (1 / [d(d+1)]) * [ sum_k Tr(M_k^dag M_k) + sum_k |Tr(M_k)|^2 ],

with Kraus operators M_k = P U0^dag G_k P, {G_k} a Kraus decomposition of
E, and P the projector onto S.  This is equivalent to the Nielsen-Horodecki
form used by QuTiP internally (``qutip.average_gate_fidelity``):

    F_bar = (d * F_pro + 1) / (d + 1),
    F_pro = (1/d^2) * sum_{i,j in S} <pi(i)| E(|i><j|) |pi(j)>,

where U0|i> = |pi(i)> when U0 is a computational-basis permutation (OR,
CCX, ...).

The key difference from the computational-basis average
F_basis = (1/d) sum_i <pi(i)| E(|i><i|) |pi(i)> of Yu et al. / Farouk
et al. is that Pedersen's formula averages over Haar-random pure
superpositions, not just the d classical basis states, which probes the
coherence of the channel off the computational diagonal.

Notation:
    * full_dims : list of subsystem dims of the simulator Hilbert space
                  (e.g. [3, 3, 4] for Cs x Cs x Rb)
    * comp_kets : list of d computational-subspace basis kets, each as a
                  full-Hilbert-space qutip.Qobj (dim = prod(full_dims))
    * U0        : d x d target unitary on the computational subspace,
                  indexed in the same order as ``comp_kets``
"""

from __future__ import annotations

from typing import Callable, List, Sequence, Tuple

import numpy as np
import qutip


def _comp_ket_matrix(comp_kets: Sequence[qutip.Qobj]) -> np.ndarray:
    """Return the (n x d) isometry V whose columns are the computational-basis
    kets embedded in the full n-dimensional Hilbert space.
    """
    V = np.zeros((comp_kets[0].shape[0], len(comp_kets)), dtype=complex)
    for i, ket in enumerate(comp_kets):
        V[:, i] = ket.full().flatten()
    return V


def build_permutation_unitary(
    comp_kets: Sequence[qutip.Qobj],
    permutation: Callable[[int], int],
) -> np.ndarray:
    """Build the d x d target unitary U0 in the computational basis given a
    permutation ``permutation: i -> pi(i)`` with U0|i> = |pi(i)>.

    Parameters
    ----------
    comp_kets : list of qutip.Qobj
        Ordered computational-basis kets.
    permutation : callable (int -> int)
        Map from input basis index to output basis index.

    Returns
    -------
    np.ndarray
        d x d unitary matrix, columns = image of the input basis vectors.
    """
    d = len(comp_kets)
    U0 = np.zeros((d, d), dtype=complex)
    for i in range(d):
        U0[permutation(i), i] = 1.0
    return U0


def channel_choi_on_subspace(
    super_full: qutip.Qobj,
    comp_kets: Sequence[qutip.Qobj],
) -> np.ndarray:
    """Evaluate the CPTP map E (given as a full-space superoperator Qobj) on
    all d^2 basis operators ``|i><j|`` of the computational subspace and
    return the rank-4 tensor

        C[i, j, k, l] = <k| E(|i><j|) |l>,        i, j, k, l = 0, ..., d-1,

    where the inner product is taken between the full-Hilbert-space
    representations of ``|k>``, ``|l>`` (columns of the isometry V whose
    columns are ``comp_kets``).

    Parameters
    ----------
    super_full : qutip.Qobj (type='super')
        Superoperator propagator of the full-space channel (e.g. the
        output of ``qutip.propagator(H, T, c_ops=c_ops, ...)``).
    comp_kets : list of qutip.Qobj
        Computational-subspace basis kets embedded in the full Hilbert
        space.  Length d.

    Returns
    -------
    np.ndarray of shape (d, d, d, d), complex
        The subspace Choi tensor ``C[i, j, k, l]``.
    """
    V = _comp_ket_matrix(comp_kets)
    d = V.shape[1]
    n = V.shape[0]
    full_dims = comp_kets[0].dims[0]  # subsystem dims of the full space

    C = np.zeros((d, d, d, d), dtype=complex)
    for i in range(d):
        v_i = V[:, i]
        for j in range(d):
            v_j = V[:, j]
            rho_in_mat = np.outer(v_i, v_j.conj())  # |i><j| on full space
            rho_in = qutip.Qobj(rho_in_mat, dims=[full_dims, full_dims])
            vec_in = qutip.operator_to_vector(rho_in)
            vec_out = super_full * vec_in
            rho_out = qutip.vector_to_operator(vec_out).full()  # n x n
            # Project into d x d computational block: V^dag rho_out V
            C[i, j, :, :] = V.conj().T @ rho_out @ V
    return C


def pedersen_average_fidelity_from_choi(
    choi_tensor: np.ndarray,
    U0: np.ndarray,
) -> Tuple[float, float, float]:
    """Compute Pedersen's Haar-averaged gate fidelity from the subspace
    Choi tensor and the target unitary.

    Uses the Nielsen-Horodecki reformulation of Pedersen Eqs. (3) and (5):

        F_pro = (1/d^2) * sum_{i,j} conj(U0[k,i]) * C[i,j,k,l] * U0[l,j],
        F_bar = (d * F_pro + 1) / (d + 1).

    In addition returns the trace-preservation indicator

        T_P = (1/d) * sum_{i, k} C[i, i, k, k],

    which equals the fraction of population that remains in the
    computational subspace after the channel (``T_P = 1`` for a perfectly
    trace-preserving channel with no leakage).

    Parameters
    ----------
    choi_tensor : np.ndarray, shape (d, d, d, d)
        Output of :func:`channel_choi_on_subspace`.
    U0 : np.ndarray, shape (d, d)
        Target unitary on the computational subspace.

    Returns
    -------
    F_bar : float
        Pedersen Haar-averaged gate fidelity.
    F_pro : float
        Process (entanglement) fidelity between the projected channel and
        ``U0``.
    survival : float
        Computational-subspace survival probability averaged over the
        d classical basis inputs.
    """
    d = choi_tensor.shape[0]
    assert choi_tensor.shape == (d, d, d, d)
    assert U0.shape == (d, d)

    # F_pro = (1/d^2) * sum_{i,j,k,l} conj(U0[k,i]) * C[i,j,k,l] * U0[l,j]
    F_pro = np.einsum("ki,ijkl,lj->", U0.conj(), choi_tensor, U0) / (d * d)
    F_pro = float(np.real(F_pro))

    F_bar = (d * F_pro + 1.0) / (d + 1.0)

    survival = float(np.real(np.einsum("iikk->", choi_tensor))) / d

    return F_bar, F_pro, survival


def pedersen_average_fidelity(
    H,
    c_ops: List[qutip.Qobj],
    t_total: float,
    comp_kets: Sequence[qutip.Qobj],
    U0: np.ndarray,
    args: dict | None = None,
    options: dict | None = None,
) -> Tuple[float, float, float, qutip.Qobj]:
    """High-level driver: compute the Pedersen Haar-averaged gate fidelity
    for a time-dependent Hamiltonian ``H`` with collapse operators
    ``c_ops`` evolved from t = 0 to t = ``t_total``, restricted to the
    computational subspace spanned by ``comp_kets`` and compared against
    the target unitary ``U0``.

    This calls :func:`qutip.propagator` to get the full-space superoperator
    propagator in one pass, then projects into the d-dim computational
    subspace and evaluates Pedersen's formula.

    Returns
    -------
    F_bar, F_pro, survival, super_full :
        Pedersen average fidelity, process fidelity, computational-subspace
        survival probability, and the full-space superoperator propagator
        (returned so the caller can inspect/reuse it).
    """
    super_full = qutip.propagator(
        H, t_total, c_ops=c_ops, args=args or {}, options=options or {}
    )
    choi = channel_choi_on_subspace(super_full, comp_kets)
    F_bar, F_pro, survival = pedersen_average_fidelity_from_choi(choi, U0)
    return F_bar, F_pro, survival, super_full


def choi_on_subspace_via_mesolve(
    H,
    c_ops: List[qutip.Qobj],
    t_total: float,
    comp_kets: Sequence[qutip.Qobj],
    args: dict | None = None,
    options: dict | None = None,
    verbose: bool = True,
) -> np.ndarray:
    """Compute the d^4 subspace Choi tensor C[i,j,k,l] = <k| E(|i><j|) |l>
    by running d^2 = 64 mesolve calls with (possibly non-Hermitian) initial
    operators ``|i><j|`` embedded in the full Hilbert space.

    This is a propagator-free, drop-in alternative to
    :func:`channel_choi_on_subspace` applied to ``qutip.propagator(...)``.
    For long gate times where the full 36^2-dim superoperator ODE is
    prohibitively expensive, this 64-run approach is typically one to two
    orders of magnitude faster because each mesolve evolves a single
    36-dim density-matrix ODE instead of the 1296-dim superoperator ODE.

    The Lindblad master equation is linear in the initial operator, so
    passing ``|i><j|`` (with ``i != j``) as the mesolve initial "state"
    is mathematically equivalent to applying the Lindblad propagator to
    that operator, even though ``|i><j|`` is not a valid density matrix.
    QuTiP's ``mesolve`` handles this correctly via vectorization.

    Parameters
    ----------
    H, c_ops, t_total, args, options :
        Standard Lindblad master-equation inputs (see ``qutip.mesolve``).
    comp_kets : list of qutip.Qobj
        Computational-subspace basis kets embedded in the full Hilbert
        space.  Length d.
    verbose : bool
        If True, print progress ticks every 8 mesolve calls.

    Returns
    -------
    np.ndarray of shape (d, d, d, d), complex
        The subspace Choi tensor ``C[i, j, k, l] = <k| E(|i><j|) |l>``.
    """
    V = _comp_ket_matrix(comp_kets)
    n, d = V.shape
    full_dims = comp_kets[0].dims[0]

    C = np.zeros((d, d, d, d), dtype=complex)
    total = d * d
    counter = 0
    for i in range(d):
        v_i = V[:, i]
        for j in range(d):
            v_j = V[:, j]
            rho_in_mat = np.outer(v_i, v_j.conj())
            rho_in = qutip.Qobj(rho_in_mat, dims=[full_dims, full_dims])
            result = qutip.mesolve(
                H,
                rho_in,
                [0.0, t_total],
                c_ops=c_ops,
                e_ops=[],
                args=args or {},
                options={**(options or {}), "store_final_state": True},
            )
            rho_out = result.final_state.full()
            C[i, j, :, :] = V.conj().T @ rho_out @ V
            counter += 1
            if verbose and counter % 8 == 0:
                print(f"    ... {counter}/{total} mesolve runs done.",
                      flush=True)
    return C

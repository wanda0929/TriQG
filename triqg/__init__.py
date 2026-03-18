"""
TriQG - Pulse-level simulation of atomic quantum systems.
"""

__version__ = "0.1.0"

from .atoms import (
    Atom,
    CsAtom,
    RbAtom,
    DIMS,
    composite_basis_state,
    composite_projector,
)

from .pulses import (
    omega_c,
    omega_p,
    omega_R,
    compute_pulse_area,
)

from .hamiltonian import build_hamiltonian

from .decoherence import build_collapse_operators

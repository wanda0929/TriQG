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
    omega_cc,
    omega_t1,
    omega_t2,
)

from .hamiltonian import build_hamiltonian

from .decoherence import build_collapse_operators

from .analysis import state_fidelity, extract_populations, average_gate_fidelity

from .solver import simulate, SimulationResult

from .visualization import plot_pulses, plot_populations, plot_populations_mc

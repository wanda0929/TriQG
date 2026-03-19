"""Smoke tests for example scripts -- verify the full pipeline runs without error."""

import numpy as np
import pytest
import matplotlib

matplotlib.use("Agg")

from triqg.atoms import CsAtom, RbAtom, composite_basis_state, composite_projector, DIMS
from triqg.pulses import omega_c, omega_p, omega_R, compute_pulse_area
from triqg.hamiltonian import build_hamiltonian
from triqg.decoherence import build_collapse_operators
from triqg.solver import simulate
from triqg.analysis import state_fidelity, extract_populations
from triqg.visualization import plot_pulses, plot_populations


class TestORGatePipelineSmoke:
    """Run the full OR gate simulation pipeline with minimal resolution."""

    # Reference parameters (Farouk et al.)
    omega_c_amp = 2 * np.pi * 50  # MHz
    omega_p_amp = 2 * np.pi * 50  # MHz
    omega_R_amp = 2.5 * omega_p_amp
    delta = 2 * np.pi * 1200  # MHz
    V_ct = 2 * np.pi * 200  # MHz (example value)

    T_c = np.pi / omega_c_amp
    T_f = 0.5  # arbitrary short value for smoke test
    sigma = 0.05

    gamma_r = 1.0 / 548.0
    gamma_R = 1.0 / 505.0
    gamma_P = 1.0 / 0.131

    def test_mesolve_pipeline_runs(self):
        """Full pipeline: build H, c_ops, run mesolve, extract populations, compute fidelity."""
        args = {
            "omega_c_amp": self.omega_c_amp,
            "omega_p_amp": self.omega_p_amp,
            "omega_R_amp": self.omega_R_amp,
            "T_c": self.T_c,
            "T_f": self.T_f,
            "sigma": self.sigma,
        }

        H = build_hamiltonian(self.delta, self.V_ct)
        c_ops = build_collapse_operators(self.gamma_r, self.gamma_R, self.gamma_P)

        cs = CsAtom()
        rb = RbAtom()
        psi0 = composite_basis_state(
            cs.level_index["1"], cs.level_index["1"], rb.level_index["A"]
        )

        # Very short tlist for speed
        t_total = 2 * self.T_c + self.T_f
        tlist = np.linspace(0, t_total, 5)

        # Build e_ops: population projectors for computational basis states
        e_ops = [
            composite_projector(2, rb.level_index["A"]),
            composite_projector(2, rb.level_index["B"]),
        ]

        result = simulate(
            method="mesolve",
            H=H,
            psi0=psi0,
            tlist=tlist,
            c_ops=c_ops,
            e_ops=e_ops,
            options={"store_final_state": True, "nsteps": 10000},
            args=args,
        )

        # Basic checks
        assert len(result.times) == 5
        assert len(result.expect) == 2

        pops = extract_populations(result, 2)
        assert pops.shape == (2, 5)

        fid = state_fidelity(result.final_state, psi0)
        assert 0.0 <= fid <= 1.0

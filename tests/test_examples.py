"""Smoke tests for example scripts -- verify the full pipeline runs without error."""

import numpy as np
import pytest
import matplotlib

matplotlib.use("Agg")

from triqg.atoms import CsAtom, RbAtom, composite_basis_state, composite_projector, DIMS
from triqg.pulses import omega_c, omega_p, omega_R, compute_pulse_area
from triqg.hamiltonian import build_hamiltonian, build_ccx_hamiltonian
from triqg.decoherence import build_collapse_operators
from triqg.solver import simulate
from triqg.analysis import state_fidelity, extract_populations, average_gate_fidelity
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


class TestCCXGatePipelineSmoke:
    """Run the full CCX gate simulation pipeline with minimal resolution."""

    omega_cc_amp = 2 * np.pi * 100  # MHz
    omega_t_amp = 2 * np.pi * 50  # MHz
    T_cc = np.pi / omega_cc_amp
    T_t = np.pi / omega_t_amp

    V_ct = 2 * np.pi * 200  # MHz

    gamma_r = 1.0 / 548.0
    gamma_R = 1.0 / 505.0
    gamma_P = 1.0 / 0.131

    def test_mesolve_pipeline_runs(self):
        """Full CCX pipeline: build H, c_ops, run mesolve, extract populations, compute fidelity."""
        args = {
            "omega_cc_amp": self.omega_cc_amp,
            "omega_t_amp": self.omega_t_amp,
            "T_cc": self.T_cc,
            "T_t": self.T_t,
        }

        H = build_ccx_hamiltonian(self.V_ct)
        c_ops = build_collapse_operators(self.gamma_r, self.gamma_R, self.gamma_P)

        cs = CsAtom()
        rb = RbAtom()
        # |1, 1, A> -- both controls in |1>, target in |A>
        # CCX should flip target to |B> (no blockade when controls in |1>)
        psi0 = composite_basis_state(
            cs.level_index["1"], cs.level_index["1"], rb.level_index["A"]
        )
        psi_target = composite_basis_state(
            cs.level_index["1"], cs.level_index["1"], rb.level_index["B"]
        )

        t_total = 2 * self.T_cc + 3 * self.T_t
        tlist = np.linspace(0, t_total, 5)

        e_ops = [
            composite_projector(2, rb.level_index["A"]),
            composite_projector(2, rb.level_index["B"]),
        ]

        # max_step is required so the adaptive ODE solver resolves the
        # piecewise pulse windows (QuTiP 5 can skip short active regions).
        max_freq = max(self.omega_cc_amp, self.omega_t_amp) / (2 * np.pi)
        max_step = 1.0 / (20 * max_freq)

        result = simulate(
            method="mesolve",
            H=H,
            psi0=psi0,
            tlist=tlist,
            c_ops=c_ops,
            e_ops=e_ops,
            options={"store_final_state": True, "nsteps": 100000, "max_step": max_step},
            args=args,
        )

        # Basic structure checks
        assert len(result.times) == 5
        assert len(result.expect) == 2

        pops = extract_populations(result, 2)
        assert pops.shape == (2, 5)

        fid = state_fidelity(result.final_state, psi_target)
        assert 0.0 <= fid <= 1.0


class TestCCXTruthTable:
    """Verify the CCX gate truth table across all 8 computational basis inputs."""

    omega_cc_amp = 2 * np.pi * 100
    omega_t_amp = 2 * np.pi * 50
    T_cc = np.pi / omega_cc_amp
    T_t = np.pi / omega_t_amp
    V_ct = 2 * np.pi * 200

    @pytest.fixture(autouse=True)
    def _setup(self):
        """Build shared Hamiltonian, args, and solver options once."""
        self.cs = CsAtom()
        self.rb = RbAtom()
        self.args = {
            "omega_cc_amp": self.omega_cc_amp,
            "omega_t_amp": self.omega_t_amp,
            "T_cc": self.T_cc,
            "T_t": self.T_t,
        }
        self.H = build_ccx_hamiltonian(self.V_ct)
        self.t_total = 2 * self.T_cc + 3 * self.T_t
        self.tlist = np.linspace(0, self.t_total, 5)
        max_freq = max(self.omega_cc_amp, self.omega_t_amp) / (2 * np.pi)
        self.max_step = 1.0 / (20 * max_freq)

    def _simulate_state(self, c1_label, c2_label, t_label):
        """Run CCX simulation for a given input state (no decoherence)."""
        psi0 = composite_basis_state(
            self.cs.level_index[c1_label],
            self.cs.level_index[c2_label],
            self.rb.level_index[t_label],
        )
        result = simulate(
            method="mesolve",
            H=self.H,
            psi0=psi0,
            tlist=self.tlist,
            c_ops=[],
            e_ops=[],
            options={
                "store_final_state": True,
                "nsteps": 100000,
                "max_step": self.max_step,
            },
            args=self.args,
        )
        return result.final_state

    # CCX truth table: (c1_in, c2_in, t_in) -> (c1_out, c2_out, t_out)
    # Only |1,1> controls flip the target; all other control combos leave it unchanged.
    CCX_TRUTH_TABLE = [
        ("0", "0", "A", "0", "0", "A"),
        ("0", "0", "B", "0", "0", "B"),
        ("0", "1", "A", "0", "1", "A"),
        ("0", "1", "B", "0", "1", "B"),
        ("1", "0", "A", "1", "0", "A"),
        ("1", "0", "B", "1", "0", "B"),
        ("1", "1", "A", "1", "1", "B"),  # flip
        ("1", "1", "B", "1", "1", "A"),  # flip
    ]

    @pytest.mark.parametrize(
        "c1_in,c2_in,t_in,c1_out,c2_out,t_out",
        CCX_TRUTH_TABLE,
        ids=[f"|{c1},{c2},{t}>" for c1, c2, t, *_ in CCX_TRUTH_TABLE],
    )
    def test_truth_table_entry(self, c1_in, c2_in, t_in, c1_out, c2_out, t_out):
        """Each computational basis input maps to the expected CCX output."""
        final = self._simulate_state(c1_in, c2_in, t_in)
        expected = composite_basis_state(
            self.cs.level_index[c1_out],
            self.cs.level_index[c2_out],
            self.rb.level_index[t_out],
        )
        fid = state_fidelity(final, expected)
        assert fid > 0.99, (
            f"Fidelity {fid:.4f} too low for "
            f"|{c1_in},{c2_in},{t_in}> -> |{c1_out},{c2_out},{t_out}>"
        )

    def test_average_gate_fidelity(self):
        """Average fidelity across all 8 basis states is high."""
        pairs = []
        for c1_in, c2_in, t_in, c1_out, c2_out, t_out in self.CCX_TRUTH_TABLE:
            final = self._simulate_state(c1_in, c2_in, t_in)
            expected = composite_basis_state(
                self.cs.level_index[c1_out],
                self.cs.level_index[c2_out],
                self.rb.level_index[t_out],
            )
            pairs.append((final, expected))
        avg_fid = average_gate_fidelity(pairs)
        assert avg_fid > 0.99, f"Average gate fidelity {avg_fid:.4f} too low"

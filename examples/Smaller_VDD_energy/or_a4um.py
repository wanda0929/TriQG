"""
OR-gate simulation at LATTICE SPACING a = 4.0 um  (vs FINAL's 5.0 um).

Levels unchanged from FINAL:  Rb 54 D_3/2 + Cs 62 D_5/2.

Only geometry changes:
    r_DA = 2.8284 um   (vs 3.5355)
    r_DD = 4.0000 um   (vs 5.0000)
    r_AA = 5.6569 um   (vs 7.0711)

C_3 and C_6 unchanged (atomic, not geometric):
    tilde C_3^{Rb-Cs}  = 10.0   GHz um^3   (calibrated, level_parameters)
    |C_6^{Rb-Rb}|      =  9.07  GHz um^6   (zero-crossing-suppressed)
    C_6^{Cs-Cs}        = -91.0  GHz um^6   (signed)

Derived interactions:
    V_ct/(2 pi) = +441.94 MHz    (vs +226.27, 1.95x stronger)
    V_DD/(2 pi) = +2.2144 MHz    (vs +0.5805, 3.81x larger)
    V_cc/(2 pi) = -2.7771 MHz    (vs -0.7280, 3.81x larger)
    R_DD = 199.6                  (vs 389.8, still >> 100)
    R_AA = 159.1                  (vs 310.8, still >> 100)

Re-tuning logic
---------------
At a = 5 um the FINAL config has M_2 = V_ct * 2 Delta / (Omega_p Omega_R) = 26
at Delta = 500 MHz.  At a = 4 um the same M_2 is achieved at
    Delta_new = Delta_old * V_ct_old / V_ct_new  ~  500 * 226/442  ~  256 MHz.
T_f scales as Delta/Omega_p^2, so
    T_f_new   ~  T_f_old * 0.51  ~  93 ns
    Gate_new  ~  2 T_c + 2 T_f_new  =  20 + 186  =  206 ns  (vs 385 ns).

Decay-floor estimate (n^3-scaled tau_r = 77 us, matching level_parameters):
    a = 5 um:  2*385/77000 = 1.00 %    (FINAL had |11,*> infidelity ~ 1.3%)
    a = 4 um:  2*206/77000 = 0.54 %    -> projected F_bar ~ 0.996

Scan plan
---------
Pass A:  Delta in {150, 200, 250, 300, 400, 500} MHz at K = 0.95
Pass B:  K in {0.90, 0.93, 0.95, 0.97, 1.00, 1.03} at best Delta of Pass A
Pass C:  side-by-side summary  (a = 5 FINAL  vs  a = 4 best)

Run with:
    cd TriQG/examples/Smaller_VDD_energy
    source ../../.venv/bin/activate
    python or_a4um.py
"""
import math
import numpy as np

from triqg.atoms import CsAtom, RbAtom, composite_basis_state
from triqg.pulses import omega_gaussian, compute_pulse_area
from triqg.hamiltonian import build_hamiltonian
from triqg.decoherence import build_collapse_operators
from triqg.solver import simulate
from triqg.analysis import state_fidelity, average_gate_fidelity

# Import the atomic constants but recompute the geometric V's at a = 4 um.
from level_parameters import (
    C3_tilde, C6_RbRb, C6_CsCs,
    gamma_r, gamma_R, gamma_P,
    tau_r_us, tau_R_us, tau_P_us,
)

# =====================================================================
# Geometry & interactions at a = 4 um  (this script's only change)
# =====================================================================
a_um = 4.00
r_DA = a_um / math.sqrt(2)             # 2.8284
r_DD = a_um                            # 4.0000
r_AA = a_um * math.sqrt(2)             # 5.6569

V_ct_MHz = 1000 * C3_tilde / r_DA**3   # +441.94
V_DD_MHz = 1000 * C6_RbRb  / r_DD**6   # +2.2144
V_cc_MHz = 1000 * C6_CsCs  / r_AA**6   # -2.7771

V_ct = 2 * np.pi * V_ct_MHz
V_DD = 2 * np.pi * V_DD_MHz
V_cc = 2 * np.pi * V_cc_MHz

R_DD = V_ct_MHz / V_DD_MHz
R_AA = V_ct_MHz / abs(V_cc_MHz)

print("=" * 78)
print(f"OR-gate at a = {a_um:.2f} um  (Rb 54 D_3/2 + Cs 62 D_5/2)")
print("=" * 78)
print(f"r_DA = {r_DA:.4f} um,  r_DD = {r_DD:.4f} um,  r_AA = {r_AA:.4f} um")
print(f"V_ct/(2 pi) = {V_ct_MHz:+9.4f} MHz")
print(f"V_DD/(2 pi) = {V_DD_MHz:+9.4f} MHz   (R_DD = {R_DD:.1f})")
print(f"V_cc/(2 pi) = {V_cc_MHz:+9.4f} MHz   (R_AA = {R_AA:.1f})")
print(f"Lifetimes (level_parameters n^3-scaled): tau_r = {tau_r_us} us, "
      f"tau_R = {tau_R_us} us, tau_P = {tau_P_us} us")
print()

# =====================================================================
# Fixed pulse parameters (same as FINAL except Delta will be scanned)
# =====================================================================
omega_p_MHz = 50.0
omega_R_MHz = 3.5 * omega_p_MHz
omega_c_MHz = 50.0

omega_p_amp = 2 * np.pi * omega_p_MHz
omega_R_amp = 2 * np.pi * omega_R_MHz
omega_c_amp = 2 * np.pi * omega_c_MHz
T_c = np.pi / omega_c_amp

ALPHA = 4.0
I_inf = 2 * 0.92770 * 2**(-1.0 / 6.0)

c_ops = build_collapse_operators(gamma_r, gamma_R, gamma_P)

cs = CsAtom(); rb = RbAtom()
ctrl_levels = [cs.level_index["0"], cs.level_index["1"]]
tgt_levels  = [rb.level_index["A"], rb.level_index["B"]]
basis_inputs = []
for c1 in range(2):
    for c2 in range(2):
        for t in range(2):
            psi_in = composite_basis_state(ctrl_levels[c1], ctrl_levels[c2], tgt_levels[t])
            psi_id = composite_basis_state(
                ctrl_levels[c1], ctrl_levels[c2],
                tgt_levels[1 - t if (c1 or c2) else t]
            )
            label = f"|{c1}{c2}{['A','B'][t]}>"
            basis_inputs.append((label, psi_in, psi_id))


def run_one(delta_MHz, K, ntime=350):
    delta = 2 * np.pi * delta_MHz
    sigma_pi4 = ((2 * np.pi * delta) / (omega_p_amp**2 * I_inf)) ** 3
    sigma = (K ** 3) * sigma_pi4
    T_f = (ALPHA * sigma) ** (1.0 / 3.0)
    t_total = 2 * T_c + 2 * T_f
    args = {"omega_c_amp": omega_c_amp, "omega_p_amp": omega_p_amp,
            "omega_R_amp": omega_R_amp, "T_c": T_c, "T_f": T_f, "sigma": sigma}
    H = build_hamiltonian(delta, V_ct, pulse_p=omega_gaussian, V_cc=V_cc)
    tlist = np.linspace(0, t_total, ntime)
    opts = {"store_final_state": True, "nsteps": 200000}
    pairs = []; per_F = []
    for lbl, psi_in, psi_id in basis_inputs:
        res = simulate(method="mesolve", H=H, psi0=psi_in, tlist=tlist,
                       c_ops=c_ops, e_ops=[], options=opts, args=args)
        pairs.append((res.final_state, psi_id))
        per_F.append(state_fidelity(res.final_state, psi_id))
    F_bar = average_gate_fidelity(pairs)
    M1 = V_ct_MHz / (omega_p_MHz**2 / (2 * delta_MHz))
    M2 = V_ct_MHz / (omega_p_MHz * omega_R_MHz / (2 * delta_MHz))
    return F_bar, per_F, T_f * 1e3, t_total * 1e3, M1, M2


# =====================================================================
# Pass A: Delta scan at K = 0.95
# =====================================================================
print("=" * 78)
print("Pass A:  Delta scan at K = 0.95,  Omega_p = 50, Omega_R = 175 MHz")
print("=" * 78)
deltas = [150.0, 200.0, 250.0, 300.0, 400.0, 500.0]
K_A = 0.95
print(f"{'Delta':>7s} {'T_f':>8s} {'gate':>8s} {'M1':>7s} {'M2':>7s} "
      f"{'F_bar':>10s} {'1-F':>10s}")
print("-" * 70)
best_A_F = 0.0; best_A_cell = None
for d in deltas:
    F, perF, Tfn, gn, M1, M2 = run_one(d, K_A)
    print(f"{d:>7.0f} {Tfn:>8.2f} {gn:>8.2f} {M1:>7.1f} {M2:>7.2f} "
          f"{F:>10.6f} {1-F:>10.3e}")
    if F > best_A_F:
        best_A_F = F; best_A_cell = (d, perF, gn, M1, M2)
print()
d_b, perF_b, gate_b, M1_b, M2_b = best_A_cell
print(f"Pass A best: F_bar = {best_A_F:.6f} at Delta = {d_b:.0f} MHz "
      f"(gate = {gate_b:.2f} ns,  M1 = {M1_b:.1f}, M2 = {M2_b:.2f})")
for (lbl, _, _), F in zip(basis_inputs, perF_b):
    print(f"   {lbl:<8s}  F = {F:.6f}")

# =====================================================================
# Pass B: K refinement at best Delta
# =====================================================================
print()
print("=" * 78)
print(f"Pass B:  K refinement at Delta = {d_b:.0f} MHz")
print("=" * 78)
Ks = [0.85, 0.90, 0.93, 0.95, 0.97, 1.00, 1.03]
best_KF = best_A_F; best_K = K_A; best_K_perF = perF_b; best_K_gate = gate_b
print(f"{'K':>5s} {'T_f':>8s} {'gate':>8s} {'F_bar':>10s} {'1-F':>10s}")
print("-" * 55)
for K in Ks:
    F, perF, Tfn, gn, M1, M2 = run_one(d_b, K)
    print(f"{K:>5.2f} {Tfn:>8.2f} {gn:>8.2f} {F:>10.6f} {1-F:>10.3e}")
    if F > best_KF:
        best_KF = F; best_K = K; best_K_perF = perF; best_K_gate = gn

print(f"\nPass B best: F_bar = {best_KF:.6f} at K = {best_K:.2f}, "
      f"gate = {best_K_gate:.2f} ns")
for (lbl, _, _), F in zip(basis_inputs, best_K_perF):
    print(f"   {lbl:<8s}  F = {F:.6f}")

# =====================================================================
# Pass C: side-by-side summary
# =====================================================================
print()
print("=" * 78)
print("Pass C:  side-by-side  a = 5.00 um (FINAL)  vs  a = 4.00 um (this run)")
print("=" * 78)
# FINAL numbers (from the report)
print(f"{'metric':<28s} {'a = 5.00 um':>15s} {'a = 4.00 um':>15s} {'gain':>10s}")
print("-" * 72)
print(f"{'V_ct/(2 pi)  (MHz)':<28s} {'226.27':>15s} {V_ct_MHz:>15.2f} "
      f"{V_ct_MHz/226.27:>10.2f}x")
print(f"{'R_DD':<28s} {'389.8':>15s} {R_DD:>15.1f} {R_DD/389.8:>10.2f}x")
print(f"{'R_AA':<28s} {'310.8':>15s} {R_AA:>15.1f} {R_AA/310.8:>10.2f}x")
print(f"{'best Delta/(2 pi)  (MHz)':<28s} {'500':>15s} {d_b:>15.0f} ")
print(f"{'gate time (ns)':<28s} {'385':>15s} {best_K_gate:>15.2f} "
      f"{385.0/best_K_gate:>9.2f}x faster")
print(f"{'F_bar OR (this run)':<28s} {'0.9924':>15s} {best_KF:>15.6f} ")
print(f"{'1 - F_bar':<28s} {'7.58e-03':>15s} {1-best_KF:>15.3e} "
      f"{7.58e-3/(1-best_KF):>9.2f}x")
print()
# Decay-floor diagnostic
floor_a5 = 2 * 385 / (tau_r_us * 1000)  # 2 T_gate / tau_r
floor_a4 = 2 * best_K_gate / (tau_r_us * 1000)
print(f"Decay floor on |11,*>  (2 T_gate / tau_r):")
print(f"   a = 5 um:  {floor_a5*100:.2f} %    (vs simulated 1 - F_(11,*) ~ 1.3 %)")
print(f"   a = 4 um:  {floor_a4*100:.2f} %    (vs simulated 1 - F_(11,*) = "
      f"{1 - best_K_perF[6]:.2%})")

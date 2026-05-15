"""
Recommended Rb-Cs Rydberg level parameters for the X-error sub-cycle
====================================================================
Companion to:  SelfCorrectingRydberg/plan/Rb_Cs_level_redesign.md
Scan script:   SelfCorrectingRydberg/python_scripts/scan_level_design.py
Date:          2026-05-14

Why this file exists
--------------------
The paper's published pair  (Rb 66 D_{5/2} + Cs 76 D_{3/2}, "Option A")
has a *same-species* problem that was not visible in Appendix A:  Rb-Rb
(D-D) blockade at r_DD = a is V_DD/(2 pi) ~ 91 MHz, which is only 1/6
of the interspecies Rb-Cs blockade.  Because the X-round Rb pulse is
*global* and species-selective, every Rb data atom is driven in
parallel, so D-D blockade enters the gate-selectivity budget on equal
footing with the published A-A entry.  The ratio

    R_DD = V_dd / V_DD  ~  517 / 91  ~  5.7        (paper, Option A)

is far below the R >= 100 needed for clean blockade selectivity.

The recommended redesign moves the Rb data state down to **Rb 54 D_{3/2}**,
which sits very close to a Rb same-species Foerster zero crossing.  The
Cs ancilla becomes **Cs 62 D_{5/2}**, chosen to keep the interspecies
Rb-Cs Foerster pair near-resonant.  Result:

    Pair                            R_DD     R_AA
    Rb 66 D_{5/2} + Cs 76 D_{3/2}   ~5.7  X  ~94   v  (paper)
    Rb 54 D_{3/2} + Cs 62 D_{5/2}   ~389 v   ~309 v  (this file)

|C_6^DD| drops by ~157x; |C_6^AA| drops by ~7.6x;  V_ct slows by ~2.3x,
which raises the gate time from ~10 us to ~23 us.  Acceptable -- still
many orders of magnitude faster than a measure-and-feedback cycle.

Foerster channel
----------------
    | Rb 54 D_{3/2} ;  Cs 62 D_{5/2} >
            <->
    | Rb 55 P_{3/2} ;  Cs 60 F_{5/2} >
    Delta_F / (2 pi) = -5.2 MHz   (~2x better resonance than Option A)

Scope
-----
This module ONLY encodes the recommended (top-1) pair from the scan.
Runner-up candidates (Rb 53 / Rb 52 / Rb 51 / Rb 50 D_{3/2}, etc.)
are listed in the design note Sec. 4 but are NOT instantiated here.
K-round (Z sub-cycle) levels are still open (see design note Sec. 6)
and likewise NOT defined here.

Usage
-----
    from level_parameters import (
        a_um, r_DA, r_DD, r_AA,
        C3_tilde, C6_RbRb, C6_CsCs, Delta_F_MHz,
        V_ct, V_DD, V_cc,
        V_ct_MHz, V_DD_MHz, V_cc_MHz,
        R_DD, R_AA,
        gamma_r, gamma_R, gamma_P,
        tau_r_us, tau_R_us, tau_P_us,
        n_R, l_R, j_R, n_r, l_r, j_r, n_P, l_P, j_P,
    )

Run directly to print a one-page summary of the redesign:
    $ python level_parameters.py
"""

import numpy as np

# =====================================================================
# Lattice geometry  (rotated toric-code unit cell)
# =====================================================================
a_um = 5.00                       # rotated-lattice spacing [um]
r_DA = a_um / np.sqrt(2)          # data-ancilla     ~ 3.5355 um
r_DD = a_um                       # data-data        = 5.00 um
r_AA = a_um * np.sqrt(2)          # ancilla-ancilla  ~ 7.0711 um

# =====================================================================
# Rydberg levels  (recommended X-round redesign)
# =====================================================================
# Rb data Rydberg state |R>:  54 D_{3/2}
n_R, l_R, j_R = 54, 2, 1.5
# Cs ancilla Rydberg state |r>:  62 D_{5/2}
n_r, l_r, j_r = 62, 2, 2.5
# Rb intermediate state |P>:  keep the paper's 7 P_{3/2}
#   (low-n, ~ns lifetime, well off-resonant from any Rydberg manifold)
n_P, l_P, j_P = 7, 1, 1.5

# =====================================================================
# Interaction coefficients
# ---------------------------------------------------------------------
# Geometry:  theta = pi/2  (quantization axis perpendicular to the array).
# Source:    python_scripts/scan_level_design.py
#            -> ARC PairStateInteractions.getC6perturbatively(..., 
#                                    degeneratePerturbation=True)
#            taking the largest-magnitude eigenvalue.
#            See plan/Rb_Cs_level_redesign.md Sec. 4.
#
# Calibration:  the raw radial-matrix-element estimate of C_3 is divided
# by 1.6 to absorb the angular projection factor that the paper's full
# treatment includes (cf. design note Sec. 3, "Calibration").  The 22.84
# vs. 36.15 ratio for Option A reproduces the paper's V_ct = 517 MHz at
# a = 5 um, so the same factor is applied here uniformly.
# =====================================================================
C3_tilde    = 10.0    # Rb-Cs Foerster effective C_3 [GHz * um^3], calibrated
C6_RbRb     =  9.07   # Rb-Rb (D-D) vdW C_6 [GHz * um^6], |value| (Foerster
                      # zero-crossing suppressed; sign not tracked here)
C6_CsCs     = -91.0   # Cs-Cs (A-A) vdW C_6 [GHz * um^6], SIGNED
Delta_F_MHz = -5.2    # Foerster defect / (2 pi) [MHz]

# ---------------------------------------------------------------------
# Derived blockade strengths
#   Time unit throughout TriQG is microseconds, so the 2 pi factors
#   convert each MHz into angular MHz = rad / us.
# ---------------------------------------------------------------------
V_ct_MHz = 1000.0 * C3_tilde / r_DA**3   # ~ +226.27 MHz  (Rb-Cs Foerster)
V_DD_MHz = 1000.0 * C6_RbRb  / r_DD**6   # ~ +0.5805 MHz  (Rb-Rb vdW)
V_cc_MHz = 1000.0 * C6_CsCs  / r_AA**6   # ~ -0.7288 MHz  (Cs-Cs vdW, signed)

V_ct = 2 * np.pi * V_ct_MHz
V_DD = 2 * np.pi * V_DD_MHz
V_cc = 2 * np.pi * V_cc_MHz

# ---------------------------------------------------------------------
# Selectivity ratios at a = 5.00 um  (target: both >= 100)
# ---------------------------------------------------------------------
R_DD = V_ct_MHz / V_DD_MHz             # ~ 389  v
R_AA = V_ct_MHz / abs(V_cc_MHz)        # ~ 310  v

# =====================================================================
# Decoherence rates at T = 300 K  (radiative + blackbody radiation)
# ---------------------------------------------------------------------
# Reference (ARC implementation):
#   [1] N. Sibalic et al., "ARC: An open-source library for calculating
#       properties of alkali Rydberg atoms",
#       Comp. Phys. Commun. 220, 319 (2017), arXiv:1612.05529.
#   [2] I. I. Beterov et al., "Quasiclassical calculations of BBR-induced
#       depopulation rates and effective lifetimes of Rydberg nS, nP, nD
#       alkali-metal atoms with n <= 80",
#       Phys. Rev. A 79, 052504 (2009), arXiv:0902.4995.
#
# Canonical call:
#       atom.getStateLifetime(n, l, j, temperature=300,
#                              includeLevelsUpTo=n+30)
#
# ESTIMATES (n^3 scaling from Option A, BBR-dominated regime):
#   Rb 54 D_{3/2}:  tau_R ~ (54/66)^3 * 134.87 us  ~  74 us
#   Cs 62 D_{5/2}:  tau_r ~ (62/76)^3 * 142.73 us  ~  77 us
#
# >>> TODO: re-verify both lifetimes with a direct ARC call before any
# >>> publication-grade simulation.  The n^3 scaling is a fine first
# >>> estimate but the quantum defect changes between n=66/n=54 (Rb)
# >>> and n=76/n=62 (Cs), so the radiative-rate prefactor shifts.
#
# The Rb intermediate state |P> = |7 P_{3/2}> keeps the paper value
# tau_P = 0.131 us  (n=7 is far below the BBR-dominated regime).
# =====================================================================
tau_r_us = 77.0     # Cs |r> = |62 D_{5/2}>,  scaled (re-verify with ARC)
tau_R_us = 74.0     # Rb |R> = |54 D_{3/2}>,  scaled (re-verify with ARC)
tau_P_us = 0.131    # Rb |P> = |7 P_{3/2}>,   paper value

gamma_r = 1.0 / tau_r_us
gamma_R = 1.0 / tau_R_us
gamma_P = 1.0 / tau_P_us


# =====================================================================
# Stand-alone summary printer
# =====================================================================
def _print_summary() -> None:
    print("=" * 72)
    print("Recommended X-round Rydberg redesign  (smaller V_DD)")
    print(f"Rb 54 D_3/2  +  Cs 62 D_5/2     at  a = {a_um:.2f} um")
    print("=" * 72)
    print("Foerster channel:")
    print("    |Rb 54 D_3/2 ; Cs 62 D_5/2>  <->  |Rb 55 P_3/2 ; Cs 60 F_5/2>")
    print(f"    Delta_F / (2 pi) = {Delta_F_MHz:+.2f} MHz")
    print()
    print("Geometry:")
    print(f"    r_DA = {r_DA:.4f} um   r_DD = {r_DD:.4f} um   r_AA = {r_AA:.4f} um")
    print()
    print("Interaction coefficients:")
    print(f"    C3_tilde (Rb-Cs)  = {C3_tilde:>8.3f}  GHz * um^3")
    print(f"    C6_RbRb  (D-D)    = {C6_RbRb:>8.3f}  GHz * um^6   (|val|, zero-crossing suppressed)")
    print(f"    C6_CsCs  (A-A)    = {C6_CsCs:>8.3f}  GHz * um^6   (signed)")
    print()
    print(f"Blockade strengths at a = {a_um:.2f} um (in MHz, /(2 pi)):")
    print(f"    V_ct / (2 pi)  = {V_ct_MHz:+9.4f} MHz   (Rb-Cs dipole-dipole, gate target)")
    print(f"    V_DD / (2 pi)  = {V_DD_MHz:+9.4f} MHz   (Rb-Rb vdW, want SMALL)")
    print(f"    V_cc / (2 pi)  = {V_cc_MHz:+9.4f} MHz   (Cs-Cs vdW, signed, want SMALL)")
    print()
    print("Selectivity ratios  (target >= 100):")
    print(f"    R_DD = V_ct / V_DD       = {R_DD:7.1f}   {'OK' if R_DD >= 100 else 'FAIL'}")
    print(f"    R_AA = V_ct / |V_cc|     = {R_AA:7.1f}   {'OK' if R_AA >= 100 else 'FAIL'}")
    print()
    print("Lifetimes at T = 300 K  (n^3-scaled estimates; re-verify with ARC):")
    print(f"    tau_r (Cs 62 D_5/2)  ~ {tau_r_us:6.2f} us   ->  gamma_r = {gamma_r:.4e} us^-1")
    print(f"    tau_R (Rb 54 D_3/2)  ~ {tau_R_us:6.2f} us   ->  gamma_R = {gamma_R:.4e} us^-1")
    print(f"    tau_P (Rb 7  P_3/2)  = {tau_P_us:6.3f} us   ->  gamma_P = {gamma_P:.4e} us^-1   (paper)")
    print()
    print("Comparison vs. paper Option A (Rb 66 D_5/2 + Cs 76 D_3/2):")
    print(f"    V_ct slowed by  {517.0 / V_ct_MHz:5.2f}x   ->  gate time grows from ~10 us to ~23 us")
    print(f"    V_DD reduced by {91.0  / V_DD_MHz:7.1f}x   ->  R_DD: 5.7  --> {R_DD:5.1f}")
    print(f"    V_AA reduced by {5.5   / abs(V_cc_MHz):7.2f}x   ->  R_AA: 94   --> {R_AA:5.1f}")
    print("=" * 72)


if __name__ == "__main__":
    _print_summary()

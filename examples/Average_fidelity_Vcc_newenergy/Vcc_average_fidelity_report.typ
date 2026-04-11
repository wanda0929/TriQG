// Option A report: Rb 66 D_5/2 + Cs 76 D_3/2 at a = 5.00 um
// with ARC-computed lifetimes and Omega_R = 2.9 * Omega_p.
//
// Compile with:  typst compile Vcc_average_fidelity_report.typ
//
// Scripts:
//   examples/Average_fidelity_Vcc_newenergy/or_average_gate_fid_gaussian.py
//   examples/Average_fidelity_Vcc_newenergy/ccx_average_gate_fidelity.py

#set document(
  title: "Option A: Rb 66 D_5/2 + Cs 76 D_3/2 at a = 5.00 um",
  author: "TriQG",
)
#set page(
  paper: "a4",
  margin: (x: 2.2cm, y: 2.4cm),
  numbering: "1 / 1",
)
#set text(font: "New Computer Modern", size: 10.5pt, lang: "en")
#set par(justify: true, leading: 0.65em)
#set heading(numbering: "1.1")
#show heading.where(level: 1): set text(size: 14pt, weight: "bold")
#show heading.where(level: 2): set text(size: 12pt, weight: "bold")
#show link: underline

#align(center)[
  #text(size: 17pt, weight: "bold")[
    Option A at $a = 5.00$ μm with ARC lifetimes
  ]

  #v(0.3em)
  #text(size: 11pt)[
    Rb $66 D_(5/2)$ + Cs $76 D_(3/2)$ Förster pair,
    $Omega_R = 2.9 thick Omega_p$,
    super-Gaussian target pulse with $sigma = 0.001771$
  ]

  #v(0.4em)
  #text(size: 10pt)[
    TriQG `examples/Average_fidelity_Vcc_newenergy` · compiled 2026-04-11
  ]
]

#v(0.6em)

#block(
  fill: rgb("#f4faff"),
  inset: 10pt,
  radius: 4pt,
  width: 100%,
)[
  *Headline results.*
  - OR gate: $overline(F)_"OR" = bold(0.991858)$,
    infidelity $1 - overline(F)_"OR" = bold(8.14 times 10^(-3))$.
  - CCX gate: $overline(F)_"CCX" = bold(0.996620)$,
    infidelity $1 - overline(F)_"CCX" = bold(3.38 times 10^(-3))$.
]

= Setup

*Rydberg pair and Förster channel.*

- Rb $|R angle.r = |66 D_(5/2) angle.r$,
  Cs $|r angle.r = |76 D_(3/2) angle.r$,
  Rb intermediate $|P angle.r = |7 P_(3/2) angle.r$.
- Förster channel:
  $|66 D_(5/2); 76 D_(3/2) angle.r <-> |67 P_(3/2); 74 F_(5/2) angle.r$.
- Source for $tilde(C)_3$ and $C_6^"CsCs"$: B. J. Ireland, J. D. Pritchard,
  J. P. Shaffer, *"Interspecies Förster resonances of Rb–Cs Rydberg
  d-states for enhanced multi-qubit gate fidelities"*,
  Phys. Rev. Research *6*, 013293 (2024), arXiv:2401.02308, Table I,
  row 1 ($theta = 90$ deg, quantization axis perpendicular to the
  array plane). Values are also tabulated locally in
  `reference/energy_level_inter.md`.

*Lattice geometry.*

- Rotated-lattice spacing $a = 5.00$ μm.
- Data--ancilla distance $r_"DA" = a \/ sqrt(2) approx 3.5355$ μm.
- Ancilla--ancilla distance $r_"AA" = a sqrt(2) approx 7.0711$ μm.

*Pulse and detuning parameters.*

- $Omega_c = Omega_p = 2 pi times 50$ MHz
- $Omega_(c c) = 2 pi times 100$ MHz, $Omega_t = 2 pi times 50$ MHz
- $Omega_R = 2.9 thick Omega_p$
- $Delta = 2 pi times 500$ MHz (two-photon detuning on the target)
- Super-Gaussian target pulse: control $pi$-pulse duration
  $T_c = pi \/ Omega_c = 10$ ns, pulse half-window $T_f = 150$ ns,
  $sigma = 0.001771$. The effective two-photon pulse area
  $integral Omega_p^2 \/ (2 Delta) dif t$ equals $pi/4$.
- Total OR gate time: $2 T_c + 2 T_f = 320$ ns.
- Total CCX gate time: $2 T_(c c) + 3 T_t = 40$ ns, with
  $T_(c c) = 5$ ns and $T_t = 10$ ns.

= ARC lifetime lookup

Rydberg-state lifetimes were computed with ARC 3.10.2 at $T = 300$ K:

#block(
  fill: luma(245),
  inset: 9pt,
  radius: 3pt,
  width: 100%,
)[
  ```python
  from arc import Rubidium, Caesium
  Rb = Rubidium()
  Cs = Caesium()
  tau_Rb_66D52 = Rb.getStateLifetime(66, 2, 2.5, temperature=300,
                                      includeLevelsUpTo=96)
  tau_Cs_76D32 = Cs.getStateLifetime(76, 2, 1.5, temperature=300,
                                      includeLevelsUpTo=106)
  ```
]

`getStateLifetime` combines the radiative Einstein-A sum with the
Beterov et al. 2009 blackbody-radiation-induced depopulation formula.
Sources:

+ N. Šibalić, J. D. Pritchard, C. S. Adams, K. J. Weatherill,
  *"ARC: An open-source library for calculating properties of alkali
  Rydberg atoms"*, Comput. Phys. Commun. *220*, 319 (2017),
  arXiv:1612.05529. Package homepage:
  #link("https://arc-alkali-rydberg-calculator.readthedocs.io").
+ I. I. Beterov, I. I. Ryabtsev, D. B. Tretyakov, V. M. Entin,
  *"Quasiclassical calculations of blackbody-radiation-induced
  depopulation rates and effective lifetimes of Rydberg $n S$, $n P$,
  and $n D$ alkali-metal atoms with $n <= 80$"*,
  Phys. Rev. A *79*, 052504 (2009), arXiv:0902.4995.

#figure(
  caption: [ARC-computed BBR-included lifetimes at $T = 300$ K.
    $tau_"rad"$ is the radiative component (at $T = 0$); $tau_"BBR"$
    is the BBR-induced depopulation time;
    $tau_"total"^(-1) = tau_"rad"^(-1) + tau_"BBR"^(-1)$.
    The intermediate state $|P angle.r = |7 P_(3/2) angle.r$ keeps the
    paper value $tau_P = 0.131$ μs.],
  kind: table,
  table(
    columns: (auto, auto, auto, auto, auto, auto),
    align: (left, center, center, center, center, center),
    stroke: 0.4pt,
    table.header(
      [*State*], [*$delta_(n l j)$*], [*$n_"eff"$*],
      [*$tau_"rad"$ (μs)*], [*$tau_"BBR"$ (μs)*], [*$bold(tau_"total")$ (μs)*],
    ),
    [Rb $|66 D_(5/2) angle.r$], [$1.3463$], [$64.65$], [$294.31$], [$248.96$], [$bold(134.87)$],
    [Cs $|76 D_(3/2) angle.r$], [$2.4755$], [$73.52$], [$260.15$], [$316.22$], [$bold(142.73)$],
  ),
)

Decoherence rates used in the master equation:
$gamma_r = 1\/142.73 thick "μs"^(-1)$ (Cs $|r angle.r$),
$gamma_R = 1\/134.87 thick "μs"^(-1)$ (Rb $|R angle.r$),
$gamma_P = 1\/0.131 thick "μs"^(-1)$ (Rb $|P angle.r$, paper value).

= Coefficients at $a = 5.00$ μm

#figure(
  caption: [Rydberg interaction strengths and dimensionless ratios for
    Option A at $a = 5.00$ μm. The blockade fidelity $P_(1 r)$ is the
    single-control-level value reported by Ireland et al. 2024 for this
    Förster channel.],
  kind: table,
  table(
    columns: (auto, auto),
    align: (left, center),
    stroke: 0.4pt,
    table.header(
      [*Quantity*], [*Value*],
    ),
    [Rydberg pair],
      [Rb $66 D_(5/2)$ + Cs $76 D_(3/2)$],
    [Förster channel],
      [$67 P_(3/2)$--$74 F_(5/2)$],
    [$tilde(C)_3$ (GHz$dot$μm³)],    [$22.84$],
    [$C_6^"CsCs"$ (GHz$dot$μm⁶)],    [$-692.9$],
    [$V_("ct") \/ (2 pi)$ (MHz)],    [$+516.81$],
    [$V_(c c) \/ (2 pi)$ (MHz)],     [$-5.54$],
    [Selectivity $V_("ct") \/ |V_(c c)|$],
                                     [$93.2$],
    [$V_("ct") \/ Delta$ (OR)],      [$1.034$],
    [$V_("ct") \/ Omega_p$ (OR)],    [$10.34$],
    [$V_("ct") \/ Omega_t$ (CCX)],   [$10.34$],
    [$|V_(c c)| \/ Omega_c$],        [$0.111$],
    [Blockade fidelity $P_(1 r)$],   [$0.9997$],
  ),
)

Both $V_("ct")$ and $V_(c c)$ are derived from the $tilde(C)_3$ and
$C_6^"CsCs"$ coefficients at the lattice distances $r_"DA"$ and $r_"AA"$:
$V_("ct") = tilde(C)_3 \/ r_"DA"^3$,
$V_(c c) = C_6^"CsCs" \/ r_"AA"^6$
(with the sign of $V_(c c)$ inherited from $C_6^"CsCs" < 0$).

= OR gate results

The OR gate implements the truth table

$ |c_1, c_2, t angle.r arrow.bar |c_1, c_2, thick t plus.circle (c_1 or c_2) angle.r, $

i.e. the target flips when at least one control is in $|1 angle.r$.
The protocol uses a control $pi$-pulse on $|1 angle.r <-> |r angle.r$,
a super-Gaussian two-photon Raman pulse on
$|A\/B angle.r <-> |R angle.r$ via the intermediate $|P angle.r$, and a
final de-excitation $pi$-pulse on the controls.

#figure(
  caption: [OR gate: per-input state fidelities for all 8 computational
    basis inputs.],
  kind: table,
  table(
    columns: (auto, auto, auto),
    align: (center, left, center),
    stroke: 0.4pt,
    table.header(
      [*\#*], [*Input*], [*$F_k$*],
    ),
    [1], [$|0,0,A angle.r$], [$0.993561$],
    [2], [$|0,0,B angle.r$], [$0.993561$],
    [3], [$|0,1,A angle.r$], [$0.995259$],
    [4], [$|0,1,B angle.r$], [$0.995259$],
    [5], [$|1,0,A angle.r$], [$0.995259$],
    [6], [$|1,0,B angle.r$], [$0.995259$],
    [7], [$|1,1,A angle.r$], [$0.983351$],
    [8], [$|1,1,B angle.r$], [$0.983351$],
  ),
)

#block(
  fill: rgb("#f4faff"),
  inset: 10pt,
  radius: 4pt,
  width: 100%,
)[
  *OR gate average fidelity:*
  $ quad overline(F)_"OR" = bold(0.991858),
    quad 1 - overline(F)_"OR" = bold(8.14 times 10^(-3)). $
]

*Per-branch observations.*

- *No blockade* ($|0,0,* angle.r$). Target populations at the end of
  the pulse read
  $(P_A, P_B, P_P, P_R) = (0.9936, 0.0034, 0.0001, 0.0030)$. The
  $approx 0.003$ residual in $|R angle.r$ is the dominant error on this
  branch; at $Omega_R = 2.9 thick Omega_p$ the EIT dark state
  $|"dark" angle.r prop Omega_R |A angle.r - Omega_p |R angle.r$ has a
  nontrivial $|R angle.r$ admixture, so imperfect adiabaticity at the
  pulse edges leaks population into the Rydberg manifold.
- *Single blockade* ($|0,1,* angle.r$, $|1,0,* angle.r$). These four
  branches sit at $F_k = 0.995259$, limited by the finite margin
  $V_("ct") \/ Delta = 1.034$ (only $approx 3%$ above the EIT-breaking
  wall at $V_("ct") \/ Delta = 1$) plus $T_1$ decay through
  $gamma_R$, $gamma_r$, $gamma_P$ over the 320 ns gate time.
- *Double blockade* ($|1,1,* angle.r$). The two branches drop to
  $F_k = 0.983351$, dominated by ARC-lifetime $T_1$ decay on the
  doubly-excited control state plus the residual $V_(c c)$ phase
  accumulated on the $|r, r angle.r$ pair; the accumulated phase sits
  near $+0.67 pi$ at this lattice spacing, so the double-blockade
  branch is the worst performer.

= CCX gate results

The CCX (Toffoli) gate flips the target only when both controls are
in $|1 angle.r$. The protocol uses two control $pi$-pulses on
$|1 angle.r <-> |r angle.r$ straddling a three-sub-pulse target
sequence on $|A\/B angle.r <-> |R angle.r$ (no intermediate $|P angle.r$
and no $Omega_R$).

#figure(
  caption: [CCX gate: per-input state fidelities for all 8 computational
    basis inputs.],
  kind: table,
  table(
    columns: (auto, auto, auto),
    align: (center, left, center),
    stroke: 0.4pt,
    table.header(
      [*\#*], [*Input*], [*$F_k$*],
    ),
    [1], [$|0,0,A angle.r$], [$0.992175$],
    [2], [$|0,0,B angle.r$], [$0.991581$],
    [3], [$|0,1,A angle.r$], [$0.996745$],
    [4], [$|0,1,B angle.r$], [$0.997920$],
    [5], [$|1,0,A angle.r$], [$0.996745$],
    [6], [$|1,0,B angle.r$], [$0.997920$],
    [7], [$|1,1,A angle.r$], [$0.999937$],
    [8], [$|1,1,B angle.r$], [$0.999937$],
  ),
)

#block(
  fill: rgb("#f4faff"),
  inset: 10pt,
  radius: 4pt,
  width: 100%,
)[
  *CCX gate average fidelity:*
  $ quad overline(F)_"CCX" = bold(0.996620),
    quad 1 - overline(F)_"CCX" = bold(3.38 times 10^(-3)). $
]

*Per-branch observations.*

- *No blockade* ($|0,0,* angle.r$) sits at
  $F_k approx 0.9919$--$0.9922$; the limiting factors are
  $V_("ct") \/ Delta$ margin and $T_1$ decay during the target sub-pulses.
- *Single blockade* ($|0,1,* angle.r$, $|1,0,* angle.r$) sits at
  $F_k approx 0.9967$--$0.9979$; there is a mild $A$/$B$ asymmetry at
  the $10^(-3)$ level from the slightly different detuned dynamics of
  the two target states.
- *Double blockade* ($|1,1,* angle.r$) is essentially perfect at
  $F_k = 0.999937$ -- the blockade is strong enough that the target
  flip is robust against both lifetime decay and $V_(c c)$ phase errors
  over the short 40 ns CCX window.

= Reproducibility

#block(
  fill: luma(246),
  inset: 8pt,
  radius: 3pt,
)[
  *Run commands* (Apple Silicon laptop, Python 3.12,
  `qutip == 5.2.3`, `arc-alkali-rydberg-calculator == 3.10.2`):
  ```text
  source .venv/bin/activate
  python examples/Average_fidelity_Vcc_newenergy/or_average_gate_fid_gaussian.py
  python examples/Average_fidelity_Vcc_newenergy/ccx_average_gate_fidelity.py
  ```
]

Solver options, timing windows, and all other numerical settings are
documented in the two script headers. The Hamiltonian construction,
collapse operators, and basis conventions live in the `triqg/` package
(`hamiltonian.py`, `decoherence.py`, `atoms.py`, `pulses.py`).

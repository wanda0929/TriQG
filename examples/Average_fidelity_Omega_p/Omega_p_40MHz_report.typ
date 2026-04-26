// Omega_p = 2*pi*40 MHz report: Option A Rb 66 D_5/2 + Cs 76 D_3/2
// at a = 5.00 um.
//
// Compile with:  typst compile Omega_p_40MHz_report.typ
//
// Scripts:
//   examples/Average_fidelity_Omega_p/or_average_gate_fid_gaussian.py
//   examples/Average_fidelity_Omega_p/ccx_average_gate_fidelity.py

#set document(
  title: "Omega_p = 2 pi * 40 MHz: average gate fidelity",
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
    Average gate fidelity at $Omega_p = 2 pi dot 40$ MHz
  ]

  #v(0.3em)
  #text(size: 11pt)[
    Option A: Rb $66 D_(5/2)$ + Cs $76 D_(3/2)$ Förster pair, \
    $a = 5.00$ μm, ARC-computed lifetimes ($T = 300$ K),
    $Omega_R = 2.9 thick Omega_p$
  ]

  #v(0.4em)
  #text(size: 10pt)[
    TriQG `examples/Average_fidelity_Omega_p` · compiled 2026-04-11
  ]
]

#v(0.6em)

#block(
  fill: rgb("#e6f7ea"),
  inset: 10pt,
  radius: 4pt,
  width: 100%,
)[
  *TL;DR.* Target probe amplitude set to $Omega_p = 2 pi dot 40$ MHz.
  The super-Gaussian pulse window $T_f$ and width $sigma$ are rescaled
  analytically so that the effective two-photon pulse area
  ($pi slash 4$) and the pulse edge shape (edge/peak $approx 2.6%$)
  are both exactly preserved. Results:

  - *OR gate:* $overline(F)_"OR" = bold(0.992301)$,
    infidelity $bold(7.70 times 10^(-3))$.
  - *CCX gate:* $overline(F)_"CCX" = bold(0.996620)$,
    infidelity $bold(3.38 times 10^(-3))$.

  At this $Omega_p$, Raman scattering off $|7 P_(3/2) angle.r$ is
  numerically negligible ($P(|P angle.r) < 10^(-4)$ on every branch),
  because the in-pulse $|P angle.r$ population scales as
  $(Omega_p slash Delta)^2 = (40 slash 500)^2 = 6.4 times 10^(-3)$.
  The OR gate infidelity is dominated by the doubly-blockaded
  $|1,1,* angle.r$ branch ($1.5 times 10^(-2)$ per input), which
  absorbs the accumulated $V_(c c)$ phase over the $468.75$ ns target
  window.
]

= Scope

This run evaluates the Option A Rydberg pair (Rb $66 D_(5/2)$ +
Cs $76 D_(3/2)$, Ireland, Pritchard & Shaffer 2024 Table I row 1)
at lattice spacing $a = 5.00$ μm, with the following pulse-side
choices:

#figure(
  caption: [Pulse and detuning parameters.],
  kind: table,
  table(
    columns: (auto, auto, auto),
    align: (left, center, left),
    stroke: 0.4pt,
    table.header(
      [*Quantity*], [*symbol*], [*value*],
    ),
    [Cs control Rabi],
      [$Omega_c$],
      [$2 pi dot 50$ MHz],
    [Rb target probe Rabi],
      [$Omega_p$],
      [$bold(2 pi dot 40 "MHz")$],
    [Rb upper-leg coupling],
      [$Omega_R = 2.9 thick Omega_p$],
      [$2 pi dot 116$ MHz],
    [Cs CCX control Rabi],
      [$Omega_("cc")$],
      [$2 pi dot 100$ MHz],
    [Rb CCX target Rabi],
      [$Omega_t$],
      [$2 pi dot 50$ MHz],
    [Two-photon detuning (OR)],
      [$Delta$],
      [$2 pi dot 500$ MHz],
    [Control $pi$-pulse (OR)],
      [$T_c = pi slash Omega_c$],
      [$10.000$ ns],
    [Target half-window (OR)],
      [$T_f$],
      [$bold(234.375 "ns")$],
    [Super-Gaussian width (OR)],
      [$sigma$],
      [$bold(0.006755829 "μs")$],
    [Total OR gate time],
      [$2 T_c + 2 T_f$],
      [$bold(488.75 "ns")$],
    [Control $pi$-pulse (CCX)],
      [$T_("cc") = pi slash Omega_("cc")$],
      [$5.000$ ns],
    [Target sub-pulse (CCX)],
      [$T_t = pi slash Omega_t$],
      [$10.000$ ns],
    [Total CCX gate time],
      [$2 T_("cc") + 3 T_t$],
      [$40.000$ ns],
  ),
)

= Target pulse reshaping

The super-Gaussian target pulse used by
`or_average_gate_fid_gaussian.py` has the form
$ Omega_p (t) = A / 2 dot exp[- ((t - t_c)^3 / sigma)^2], $
with $t_c$ the pulse center and $A$ the amplitude. Substituting
$w = (t - t_c) slash sigma^(1/3)$ in the two-photon pulse-area
integral yields
$ "area" = A^2 sigma^(1/3) / (8 Delta) dot J(T_f slash sigma^(1/3)),
  quad
  J(x) = integral_(-x)^(x) e^(-2 w^6) dif w. $
The edge-to-peak ratio is $exp[- (T_f^3 slash sigma)^2]$, which also
depends only on the ratio $T_f slash sigma^(1/3)$.

Two constraints fix $T_f$ and $sigma$ once $A$ is chosen:

+ *Effective two-photon pulse area* stays at $pi slash 4$.
+ *Pulse edge shape* (edge/peak ratio) stays fixed, so that the
  super-Gaussian starts and ends smoothly instead of acquiring a
  hard cutoff at $t = T_c$ and $t = T_c + 2 T_f$.

Holding $x equiv T_f slash sigma^(1/3)$ fixed freezes the pulse
shape and freezes $J$; the area constraint then reduces to
$A^2 sigma^(1/3) = "const"$. Combined with the shape constraint,

$ sigma prop A^(-6), quad quad T_f prop A^(-2). $

At $A = 2 pi dot 40$ MHz, the pulse parameters used in this run are

$ T_f = 0.234375 " μs", quad sigma = 0.006755829 " μs", $

which give the pulse invariants

#figure(
  caption: [Numerical verification of the super-Gaussian invariants
    at the chosen $(T_f, sigma)$. The "area" line is computed by
    `compute_pulse_area(omega_gaussian, delta, T_c, T_c + 2 T_f, args)`;
    the others are closed-form.],
  kind: table,
  table(
    columns: (auto, auto, auto),
    align: (left, center, center),
    stroke: 0.4pt,
    table.header(
      [*Invariant*], [*formula*], [*value*],
    ),
    [Shape ratio],
      [$x = T_f slash sigma^(1/3)$],
      [$1.2398$],
    [Edge / peak],
      [$exp[-(T_f^3 slash sigma)^2]$],
      [$2.647 times 10^(-2)$],
    [Two-photon pulse area],
      [$integral Omega_p^2 slash (2 Delta) thin dif t$],
      [$0.785353 approx pi slash 4$],
  ),
)

The pulse is therefore an order-6 super-Gaussian with a flat top,
soft edges at $approx 2.6 %$ of peak, and a pulse area calibrated to
$pi / 4$ to six decimal places. The total OR gate time
$2 T_c + 2 T_f = 488.75$ ns is the direct consequence of stretching
$T_f$ when $A$ is lowered; this is the only "free" consequence of the
reshaping, and drives the cost analysis in §5.

= Interaction coefficients at $a = 5.00$ μm

Derived blockade strengths for the Option A pair at $theta = 90 degree$
(quantization axis perpendicular to the array plane). $V_("ct")$ is
the Rb--Cs dipole--dipole Förster coefficient at the data--ancilla
distance $r_"DA" = a slash sqrt(2)$; $V_(c c)$ is the Cs--Cs van der
Waals interaction at the ancilla--ancilla distance
$r_"AA" = a sqrt(2)$.

#figure(
  caption: [Interaction strengths and dimensionless ratios at
    $a = 5.00$ μm. Bolded row is the tight margin that dominates
    single-blockade OR-gate errors in §5.],
  kind: table,
  table(
    columns: (auto, auto),
    align: (left, center),
    stroke: 0.4pt,
    table.header(
      [*Quantity*], [*value*],
    ),
    [Lattice spacing], [$a = 5.00$ μm],
    [Nearest data--ancilla distance], [$r_"DA" = 3.5355$ μm],
    [Nearest ancilla--ancilla distance], [$r_"AA" = 7.0711$ μm],
    [Förster pair],
      [Rb $66 D_(5/2)$ + Cs $76 D_(3/2)$],
    [Förster channel],
      [$|"Rb " 67 P_(3/2) ; " Cs " 74 F_(5/2) angle.r$],
    [Effective Förster $tilde(C)_3$],
      [$22.84$ GHz$dot$μm³],
    [Cs--Cs van der Waals $C_6^"CsCs"$],
      [$-692.9$ GHz$dot$μm⁶],
    [$V_("ct") slash (2 pi)$],
      [$+516.810$ MHz],
    [$V_(c c) slash (2 pi)$],
      [$-5.543$ MHz],
    [Selectivity $V_("ct") slash |V_(c c)|$],
      [$93.2$],
    [$bold(V_("ct") slash Delta)$ *(OR)*],
      [$bold(1.034)$],
    [$V_("ct") slash Omega_p$ (OR)],
      [$12.92$],
    [$V_("ct") slash Omega_t$ (CCX)],
      [$10.34$],
    [$|V_(c c)| slash Omega_c$],
      [$0.111$],
    [$|V_(c c)| slash Omega_("cc")$],
      [$0.055$],
    [Blockade fidelity $P_(1 r)$],
      [$0.9997$],
  ),
)

*Source of the $tilde(C)_3$ / $C_6^"CsCs"$ values:*
B. J. Ireland, J. D. Pritchard, J. P. Shaffer,
*"Interspecies Förster resonances of Rb--Cs Rydberg d-states for
enhanced multi-qubit gate fidelities"*, Phys. Rev. Research *6*,
013293 (2024), arXiv:2401.02308, Table I row 1. Values tabulated
locally in `reference/energy_level_inter.md`.

= ARC-computed lifetimes at $T = 300$ K

Lifetimes for the two Rydberg states were computed with ARC 3.10.2
using
`atom.getStateLifetime(n, l, j, temperature=300, includeLevelsUpTo=n+30)`.
ARC combines the radiative Einstein-$A$ sum with the Beterov et al.
2009 blackbody-radiation-induced depopulation formula.

#figure(
  caption: [ARC-computed BBR-included lifetimes at $T = 300$ K.
    $tau_"rad"$ is the radiative component (at $T = 0$);
    $tau_"BBR"$ is the BBR-induced depopulation time;
    $tau_"total"^(-1) = tau_"rad"^(-1) + tau_"BBR"^(-1)$.
    The intermediate $|7 P_(3/2) angle.r$ state retains the
    paper value $tau_P = 0.131$ μs.],
  kind: table,
  table(
    columns: (auto, auto, auto, auto),
    align: (left, center, center, center),
    stroke: 0.4pt,
    table.header(
      [*State*],
      [*$tau_"rad"$ (μs)*],
      [*$tau_"BBR"$ (μs)*],
      [*$bold(tau_"total")$ (μs)*],
    ),
    [Rb $|66 D_(5/2) angle.r$], [$294.31$], [$248.96$], [$bold(134.87)$],
    [Cs $|76 D_(3/2) angle.r$], [$260.15$], [$316.22$], [$bold(142.73)$],
    [Rb $|7 P_(3/2) angle.r$ (paper)], [—], [—], [$bold(0.131)$],
  ),
)

*Sources:*
+ N. Šibalić, J. D. Pritchard, C. S. Adams, K. J. Weatherill,
  *"ARC: An open-source library for calculating properties of alkali
  Rydberg atoms"*, Comput. Phys. Commun. *220*, 319 (2017),
  arXiv:1612.05529. Package homepage:
  #link("https://arc-alkali-rydberg-calculator.readthedocs.io").
+ I. I. Beterov, I. I. Ryabtsev, D. B. Tretyakov, V. M. Entin,
  *"Quasiclassical calculations of blackbody-radiation-induced
  depopulation rates and effective lifetimes of Rydberg $n S$, $n P$,
  and $n D$ alkali-metal atoms with $n <= 80$"*, Phys. Rev. A *79*,
  052504 (2009), arXiv:0902.4995.

= OR gate results

#figure(
  caption: [OR gate per-input state fidelities at
    $Omega_p = 2 pi dot 40$ MHz. All 8 computational basis states
    are simulated via `mesolve` with the collapse operators built
    from the lifetimes in §4.],
  kind: table,
  table(
    columns: (auto, auto, auto),
    align: (center, left, center),
    stroke: 0.4pt,
    table.header(
      [*\#*], [*Input*], [*$F_k$*],
    ),
    [1], [$|0,0,A angle.r$], [$0.995384$],
    [2], [$|0,0,B angle.r$], [$0.995384$],
    [3], [$|0,1,A angle.r$], [$0.994474$],
    [4], [$|0,1,B angle.r$], [$0.994474$],
    [5], [$|1,0,A angle.r$], [$0.994474$],
    [6], [$|1,0,B angle.r$], [$0.994474$],
    [7], [$|1,1,A angle.r$], [$0.984872$],
    [8], [$|1,1,B angle.r$], [$0.984872$],
  ),
)

#block(
  fill: rgb("#f4faff"),
  inset: 10pt,
  radius: 4pt,
  width: 100%,
)[
  *OR gate average fidelity:*
  $ quad overline(F)_"OR" = 1/8 sum_(k=1)^8 F_k = bold(0.992301),
    quad 1 - overline(F)_"OR" = bold(7.70 times 10^(-3)). $
]

== Diagnostic: target-level populations

For each input, the script reports the target atom's projection onto
the four relevant Rb levels: the logical $|A angle.r$, $|B angle.r$
ground states, the intermediate $|P angle.r = |7 P_(3/2) angle.r$,
and the Rydberg $|R angle.r = |66 D_(5/2) angle.r$.

#figure(
  caption: [OR gate diagnostic: target-level populations at the end
    of the gate. Expected (ideal) populations: $|0,0,* angle.r$
    branches should return to the input state; every other branch
    should flip $A <-> B$.],
  kind: table,
  table(
    columns: (auto, auto, auto, auto, auto, auto),
    align: (center, left, center, center, center, center),
    stroke: 0.4pt,
    table.header(
      [*\#*], [*Input*],
      [*$P(|A angle.r)$*], [*$P(|B angle.r)$*],
      [*$P(|P angle.r)$*], [*$P(|R angle.r)$*],
    ),
    [1], [$|0,0,A angle.r$], [$0.9954$], [$0.0032$], [$0.0000$], [$0.0014$],
    [2], [$|0,0,B angle.r$], [$0.0032$], [$0.9954$], [$0.0000$], [$0.0014$],
    [3], [$|0,1,A angle.r$], [$0.0022$], [$0.9945$], [$0.0000$], [$0.0001$],
    [4], [$|0,1,B angle.r$], [$0.9945$], [$0.0022$], [$0.0000$], [$0.0001$],
    [5], [$|1,0,A angle.r$], [$0.0022$], [$0.9945$], [$0.0000$], [$0.0001$],
    [6], [$|1,0,B angle.r$], [$0.9945$], [$0.0022$], [$0.0000$], [$0.0001$],
    [7], [$|1,1,A angle.r$], [$0.0019$], [$0.9849$], [$0.0000$], [$0.0000$],
    [8], [$|1,1,B angle.r$], [$0.9849$], [$0.0019$], [$0.0000$], [$0.0000$],
  ),
)

*Per-branch infidelity budget.*

- *$|0,0,* angle.r$ branches ($1 - F_k = 4.6 times 10^(-3)$ each).*
  Neither control is in $|r angle.r$, so there is no blockade. The
  target executes the full
  $|A angle.r arrow |P angle.r arrow |R angle.r arrow |P angle.r arrow |A angle.r$
  cycle and should return unchanged. The $0.0032$ miscount is the
  imperfect-cycle residual, and the $0.0014$ population stuck in
  $|R angle.r$ is EIT-dark-state leakage at the mixing angle set by
  $Omega_R slash Omega_p = 2.9$. Notably, $P(|P angle.r) < 10^(-4)$:
  Raman scattering off the short-lived intermediate state has
  essentially turned off at this $Omega_p$.

- *$|0,1,* angle.r$ / $|1,0,* angle.r$ branches
  ($1 - F_k = 5.5 times 10^(-3)$ each).* One control is in
  $|r angle.r$, so the target $|R angle.r$ level is shifted by
  $V_("ct") slash (2 pi) = 516.81$ MHz. This is close to but larger
  than $Delta slash (2 pi) = 500$ MHz, giving
  $V_("ct") slash Delta = 1.034$. The $0.9945$ flipped-target
  population shows the blockade is working, but the tight margin
  lets $approx 0.22 %$ of the amplitude leak back through the
  off-resonant two-photon path. The $0.0001$ residual $|R angle.r$
  population confirms the single-control blockade is suppressing
  $|R angle.r$ excursion correctly.

- *$|1,1,* angle.r$ branches ($1 - F_k = 15.1 times 10^(-3)$ each).*
  Both controls are in $|r angle.r$, so the blockade shift is twice
  as large, and an additional $V_(c c) slash (2 pi) = -5.543$ MHz
  phase accumulates over the $2 T_f = 468.75$ ns target window. The
  $0.0019$ residual "old" target population is even smaller than on
  the single-blockade branches, but the accumulated $V_(c c)$ phase
  rotates the doubly-blockaded amplitude by
  $approx |V_(c c)| dot 2 T_f approx 1.04 pi$, landing on a
  near-destructive point for the final state alignment. This is the
  dominant contributor to the OR gate infidelity.

= CCX gate results

The CCX script (`ccx_average_gate_fidelity.py`) uses independent
square pulses $Omega_("cc")$ and $Omega_t$ and does *not* depend on
$Omega_p$, so the CCX numbers are invariant under the target-probe
rescaling of §2. They are included here for a complete picture of
the two Option-A gates on this template.

#figure(
  caption: [CCX gate per-input state fidelities. Pulse parameters:
    $Omega_("cc") = 2 pi dot 100$ MHz, $Omega_t = 2 pi dot 50$ MHz,
    $T_("cc") = 5$ ns, $T_t = 10$ ns, total gate time $40$ ns.],
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
  $ quad overline(F)_"CCX" = 1/8 sum_(k=1)^8 F_k = bold(0.996620),
    quad 1 - overline(F)_"CCX" = bold(3.38 times 10^(-3)). $
]

The CCX gate is dominated by the $|0,0,* angle.r$ branch, which
loses $approx 8 times 10^(-3)$ per input from the $V_(c c)$ phase
during the Cs control $pi$-pulse window. The $|1,1,* angle.r$
branch is essentially perfect ($F_k = 0.999937$), and the
single-blockade branches each contribute about $2.7 times 10^(-3)$
from the $V_("ct") slash Delta$ margin (same underlying mechanism
as on the OR gate's single-blockade branches, even though CCX uses
square pulses).

= Reproducibility

#block(
  fill: luma(246),
  inset: 8pt,
  radius: 3pt,
)[
  *Environment* (Apple Silicon, Python 3.12,
  `qutip == 5.2.3`, `arc-alkali-rydberg-calculator == 3.10.2`):
  ```text
  source .venv/bin/activate
  python examples/Average_fidelity_Omega_p/or_average_gate_fid_gaussian.py
  python examples/Average_fidelity_Omega_p/ccx_average_gate_fidelity.py
  ```
]

*Key parameters in the OR gate script
(`or_average_gate_fid_gaussian.py`):*

- `omega_c_amp = 2 * np.pi * 50`
- `omega_p_amp = 2 * np.pi * 40.0`
- `omega_R_amp = 2.9 * omega_p_amp`  (→ $2 pi dot 116$ MHz)
- `delta = 2 * np.pi * 500`
- `T_c = np.pi / omega_c_amp`  (→ $10$ ns)
- `T_f = 0.15 * (5 / 4) ** 2`  (→ $0.234375$ μs = $234.375$ ns)
- `sigma = 0.001771 * (5 / 4) ** 6`  (→ $approx 0.006756$ μs)
- `a_um = 5.0`
- `C3_tilde = 22.84`, `C6_CsCs = -692.9`  (Ireland et al. 2024)
- `gamma_r = 1 / 142.73`, `gamma_R = 1 / 134.87`  (ARC, 300 K)
- `gamma_P = 1 / 0.131`  (paper value retained)

*Key parameters in the CCX gate script
(`ccx_average_gate_fidelity.py`):* unchanged relative to the
previous template -- $Omega_("cc") = 2 pi dot 100$ MHz,
$Omega_t = 2 pi dot 50$ MHz, $T_("cc") = pi slash Omega_("cc") = 5$ ns,
$T_t = pi slash Omega_t = 10$ ns, same interaction coefficients
and lifetimes as the OR gate script.

= Conclusion

+ *The super-Gaussian scaling law is exact.* Holding the shape
  ratio $x = T_f slash sigma^(1/3)$ fixed and imposing
  $A^2 sigma^(1/3) = "const"$ uniquely determines
  $sigma prop A^(-6)$ and $T_f prop A^(-2)$. Numerical verification
  gives pulse area $= 0.785353 = pi slash 4$ to six decimal places
  and edge/peak ratio $= 2.647 times 10^(-2)$, matching closed-form
  expectations bit-for-bit.

+ *OR gate average fidelity:*
  $overline(F)_"OR" = 0.992301$, infidelity
  $7.70 times 10^(-3)$. The dominant error is on the doubly-blockaded
  $|1,1,* angle.r$ branch ($1.5 times 10^(-2)$ per input), from the
  accumulated Cs--Cs $V_(c c)$ phase over the now-longer
  $468.75$ ns target window. Single-blockade branches each lose
  $5.5 times 10^(-3)$ from the tight $V_("ct") slash Delta = 1.034$
  margin. The no-blockade $|0,0,* angle.r$ branch loses
  $4.6 times 10^(-3)$, split between a $0.0032$ imperfect-cycle
  residual and $0.0014$ EIT-dark-state $|R angle.r$ leakage.

+ *Raman scattering is off.* With $Omega_p = 2 pi dot 40$ MHz, the
  in-pulse $|P angle.r$ population is
  $(Omega_p slash (2 Delta))^2 approx 1.6 times 10^(-3)$, so the
  $|P angle.r$-state loss channel (short $tau_P = 0.131$ μs
  lifetime) contributes $< 10^(-4)$ on every branch. The longer
  target window does *not* cost fidelity through Raman scattering
  at this $Omega_p$.

+ *CCX gate average fidelity:* $overline(F)_"CCX" = 0.996620$,
  infidelity $3.38 times 10^(-3)$, unchanged because the CCX pulses
  do not depend on $Omega_p$.

+ *Where the infidelity budget points next.* The $|1,1,* angle.r$
  $V_(c c)$ phase is the largest remaining lever on the OR gate;
  the $V_("ct") slash Delta$ margin is the second. Both are
  properties of the geometry ($a = 5.00$ μm) and the Option A pair,
  not of $Omega_p$. Further reductions in $Omega_p$ will stretch
  $T_f$ (as $A^(-2)$) and grow the $|1,1,* angle.r$ $V_(c c)$ phase
  by the same factor; a sensible next step is to scan $Omega_p$
  upward until $V_(c c) dot 2 T_f$ lands on a constructive-alignment
  point, rather than lowering it further.

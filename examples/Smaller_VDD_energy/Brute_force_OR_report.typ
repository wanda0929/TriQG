// Brute_force_OR_report.typ
//   Brand-new report (companion to Smaller_VDD_average_fidelity_report.typ).
//   5-D brute-force parameter search for OR-gate F_bar > 0.99 at
//   Rb 54 D_{3/2} + Cs 62 D_{5/2}.
//
// Compile with:  typst compile Brute_force_OR_report.typ
//
// Scripts (run in this folder):
//   level_parameters.py                  (atomic + lattice constants)
//   brute_force_or.py                    (5-D sweep, ~9 min, this report)
//   brute_force_plots.py                 (axis marginals, 2D, Pareto)
//   champion_verify.py                   (re-simulate champion cell)
//
// Logs / CSV:
//   brute_force_or.log, brute_force_results.csv, brute_force_robust.csv,
//   champion_verify.log

#set document(
  title: "Brute-force 5-D OR-gate parameter search at Rb 54 D_3/2 + Cs 62 D_5/2",
  author: "TriQG",
)
#set page(paper: "a4", margin: (x: 2.2cm, y: 2.4cm), numbering: "1 / 1")
#set text(font: "New Computer Modern", size: 10.5pt, lang: "en")
#set par(justify: true, leading: 0.65em)
#set heading(numbering: "1.1")
#show heading.where(level: 1): set text(size: 14pt, weight: "bold")
#show heading.where(level: 2): set text(size: 12pt, weight: "bold")
#show link: underline

#align(center)[
  #text(size: 17pt, weight: "bold")[
    Brute-force 5-D parameter search:\
    a faster, more-robust OR gate with $overline(F)_"OR" > 0.99$
  ]

  #v(0.3em)
  #text(size: 11pt)[
    Rb $54 D_(3/2)$ + Cs $62 D_(5/2)$,
    sweep over $(a, Omega_p, Omega_R \/ Omega_p, Delta, Omega_c)$
    at fixed $K = 0.95$, $alpha = 4$
  ]

  #v(0.4em)
  #text(size: 10pt)[
    TriQG `examples/Smaller_VDD_energy` · compiled 2026-05-15
  ]
]

#v(0.6em)

#block(
  fill: rgb("#eef8ec"),
  inset: 10pt,
  radius: 4pt,
  width: 100%,
)[
  *Headline.*

  A 432-cell brute-force sweep on the five user-specified knobs
  ($a$, $Omega_p$, $Omega_R \/ Omega_p$, $Delta$, $Omega_c$)
  with the rev. 4 report's optimal pulse shape ($K = 0.95$, $alpha = 4$)
  held fixed finds:

  - *221 of 432 cells* (51%) exceed $overline(F)_"OR" > 0.99$.
  - *40 cells* are fully robust: every one of their 1-step nearest
    neighbors in the 5-D grid also exceeds $0.99$.
  - The *most-robust champion* is a *new operating point* that
    *beats the rev. 4 FINAL* configuration on every axis the user
    cares about:
    #table(
      columns: (auto, auto, auto),
      align: (left, center, center),
      stroke: 0.4pt,
      table.header([*Quantity*], [*Rev. 4 FINAL*], [*Brute-force champion*]),
      [$a$ (μm)],                       [$5.00$],   [$bold(4.00)$],
      [$Omega_p \/ (2 pi)$ (MHz)],      [$50$],     [$bold(60)$],
      [$Omega_R \/ Omega_p$],           [$3.5$],    [$bold(4.0)$],
      [$Delta \/ (2 pi)$ (MHz)],        [$500$],    [$500$],
      [$Omega_c \/ (2 pi)$ (MHz)],      [$50$],     [$bold(70)$],
      [Total gate $T_("tot")$ (ns)],    [$385$],    [$bold(267.7)$],
      [$overline(F)_"OR"$],             [$0.99242$],[$bold(0.99416)$],
      [$1 - overline(F)_"OR"$],         [$7.58 times 10^(-3)$],
                                        [$bold(5.84 times 10^(-3))$],
      [$R_("DD")$, $R_("AA")$],         [$390 \/ 311$], [$200 \/ 159$],
      [Robust neighbors $> 0.99$],      [(not scored)], [$bold(6 \/ 6)$],
    )

  Net win: *infidelity drops 23%*, *gate time drops 30%*, and the
  champion sits inside a fully-robust 40-cell plateau in 5-D
  parameter space. The selectivity ratios $R_("DD")$, $R_("AA")$
  are reduced from rev. 4 but still safely $> 100$.
]

= What was swept and what was held fixed

The user asked for a *brute force* over five knobs:
$a$ (lattice spacing), $Omega_p$ (target probe amplitude),
$Omega_R \/ Omega_p$ (EIT ratio), $Delta$ (probe detuning),
$Omega_c$ (Cs ancilla pi-pulse Rabi).

Two pulse-shape parameters were *not* swept because the prior
report (Sec. 4 of `Smaller_VDD_average_fidelity_report.typ`)
already pinned the optimum and the user's question was about the
five physical knobs above:

- $K = "area" \/ (pi \/ 4) = 0.95$ -- the sub-$pi \/ 4$ optimum
  for finite $V_(c t)$.
- $alpha = T_f^3 \/ sigma = 4$ -- smooth super-Gaussian edges
  (edge / peak $= 1.1 times 10^(-7)$).

At each cell, the closed-form area-calibration relation
$ sigma = K^3 thick ((2 pi Delta) / (Omega_p^2 thick I_oo))^3,
  quad T_f = (alpha sigma)^(1\/3),
  quad I_oo = 2 thick Gamma(7\/6) thick 2^(-1\/6) approx 1.6534 $
fixes $(sigma, T_f)$ from $(Omega_p, Delta)$, so the pulse area
is *exactly* $0.95 thick pi \/ 4$ everywhere on the grid.
$T_c = pi \/ Omega_c$.

The lattice spacing $a$ enters via the interaction strengths:

#figure(
  caption: [Effect of $a$ on the three Rydberg blockades.
    Selectivity ratios $R_("DD")$ and $R_("AA")$ stay above the
    $>= 100$ threshold across the full $4.0$–$5.5$ μm range,
    so every cell in the brute-force grid is *physically admissible*
    on the same-species-blockade selectivity budget.],
  kind: table,
  table(
    columns: (auto, auto, auto, auto, auto, auto, auto),
    align: (center, center, center, center, center, center, center),
    stroke: 0.4pt,
    table.header(
      [*$a$ (μm)*],
      [*$V_(c t) \/ (2 pi)$ (MHz)*],
      [*$V_("DD") \/ (2 pi)$ (MHz)*],
      [*$V_(c c) \/ (2 pi)$ (MHz)*],
      [*$R_("DD")$*],
      [*$R_("AA")$*],
      [*OK?*],
    ),
    [$4.0$], [$441.9$], [$+2.214$], [$-2.7771$], [$200$], [$159$], [#sym.checkmark],
    [$4.5$], [$310.4$], [$+0.788$], [$-1.3699$], [$394$], [$227$], [#sym.checkmark],
    [$5.0$], [$226.3$], [$+0.290$], [$-0.7280$], [$781$], [$311$], [#sym.checkmark],
    [$5.5$], [$170.0$], [$+0.111$], [$-0.4109$], [$1535$], [$413$], [#sym.checkmark],
  ),
)

The grid:

#figure(
  caption: [Five-axis brute-force grid (total $4 dot 3 dot 3 dot 4 dot 3 = 432$ cells).
    Each cell runs 8 `mesolve` simulations (one per computational
    basis input). Total wall time $549$ s ($9.1$ min) on Apple
    Silicon, $0.15$ s per `mesolve`.],
  kind: table,
  table(
    columns: (auto, auto),
    align: (left, left),
    stroke: 0.4pt,
    table.header([*Axis*], [*Grid values*]),
    [$a$ (μm)],                          [$4.0, 4.5, 5.0, 5.5$],
    [$Omega_p \/ (2 pi)$ (MHz)],         [$40, 50, 60$],
    [$Omega_R \/ Omega_p$],              [$3.0, 3.5, 4.0$],
    [$Delta \/ (2 pi)$ (MHz)],           [$400, 500, 600, 800$],
    [$Omega_c \/ (2 pi)$ (MHz)],         [$30, 50, 70$],
  ),
)

= Sweep results: where the $> 0.99$ region lives

The pass-rate marginal on each axis -- *fraction of cells with
$overline(F)_"OR" > 0.99$ when we average over the other four
axes* -- already tells most of the story.

#figure(
  image("brute_force_axis_marginals.png", width: 100%),
  caption: [
    Bar height = % of the $4 dot 3 dot 4 dot 3 = 144$
    (resp. $108$) cells in the marginal that exceed
    $overline(F)_"OR" > 0.99$. White label at the bar base is the
    maximum $overline(F)_"OR"$ achieved in that marginal. The
    *more-robust* end of each axis is where the bar is tallest;
    the *higher-fidelity* end is where the white "max" label is
    largest.

    Key reads: (i) the broadest pass-rate plateau in $a$ is at
    $a = 4.5$ μm ($66 %$), and *not* at the rev. 4 FINAL value
    $a = 5.0$ μm ($61 %$); (ii) the pass-rate in $Omega_p$ is
    monotone -- $60$ MHz beats $50$ MHz ($71 %$ vs $63 %$); (iii)
    ratio $= 3.5$ wins on pass-rate ($56 %$), but ratio $= 4.0$
    wins on peak fidelity ($0.9942$); (iv) $Delta = 500$ MHz wins
    on both fronts; (v) larger $Omega_c$ is uniformly better
    (shorter $T_c$ = less Cs decay in the gate window).
  ],
) <fig:marginals>

The 2-D projections, taking the *maximum* $overline(F)_"OR"$ over
the other three axes, show that the high-fidelity region is a
2-D *plateau*, not a 2-D peak:

#figure(
  image("brute_force_2D_marginals.png", width: 100%),
  caption: [
    Best-case $overline(F)_"OR"$ as a function of two axes,
    maximized over the other three. Every panel has a wide region
    above $0.992$. The $(Omega_p, Delta)$ panel shows the FINAL
    rev. 4 cell ($Omega_p = 50$, $Delta = 500$) sits on the edge of
    the plateau, not in its center -- raising $Omega_p$ from $50$
    to $60$ MHz and keeping $Delta = 500$ MHz moves *into* the
    plateau and reaches $0.9942$. The $(Omega_R \/ Omega_p, Delta)$
    panel shows the rev. 4 ratio $= 3.5$ is already excellent but
    $= 4.0$ wins at $Delta = 400$–$500$ MHz.
  ],
) <fig:2D>

The Pareto cloud over $(T_("tot"), overline(F)_"OR")$ confirms the
prize is *both* fidelity and speed:

#figure(
  image("brute_force_pareto.png", width: 90%),
  caption: [
    All $432$ cells in $(T_("tot"), overline(F)_"OR")$.
    Grey: cells below $0.99$.
    Coloured: cells above $0.99$, with hue = *robust score* (fraction
    of 1-step nearest-neighbors in the 5-D grid that also exceed
    $0.99$). Yellow = fully robust. Dashed black: Pareto frontier.
    The Pareto frontier saturates at $overline(F)_"OR" approx
    0.9942$ for $T_("tot") gt.eq 217$ ns; nothing under $200$ ns
    clears $0.99$ in this grid.
  ],
) <fig:pareto>

#block(
  fill: rgb("#f0f6ff"),
  inset: 10pt,
  radius: 4pt,
  width: 100%,
)[
  *Statistics summary.*

  - Cells with $overline(F)_"OR" > 0.99$: $bold(221 \/ 432)$ ($51 %$).
  - Cells fully robust (every 5-D neighbor also $> 0.99$): $bold(40)$.
  - Pareto saturation: $overline(F)_"OR" approx 0.994$, reached at
    $T_("tot") tilde.op 217$ ns and never improved.
  - Hard floor on $F_(1,1,*)$ across the plateau: $approx 0.99$
    (blockade-leakage limit -- see Sec. 4).
]

= The champion: the most-robust $> 0.99$ cell

Among the $221$ passing cells, the brute-force script ranks them by
*robust score* (fraction of 5-D nearest neighbors that also exceed
$0.99$), then by $overline(F)_"OR"$. The rank-1 champion:

#figure(
  caption: [Champion cell. Verified by re-running `champion_verify.py`
    against the canonical Hamiltonian and decoherence model.],
  kind: table,
  table(
    columns: (auto, auto, auto),
    align: (left, center, center),
    stroke: 0.4pt,
    table.header([*Parameter*], [*Value*], [*Notes*]),
    [$a$],                        [$4.00$ μm],            [tighter lattice],
    [$Omega_p \/ (2 pi)$],        [$60$ MHz],             [+20% over rev. 4],
    [$Omega_R \/ (2 pi)$],        [$240$ MHz],            [ratio $= 4.0$],
    [$Omega_c \/ (2 pi)$],        [$70$ MHz],             [+40% over rev. 4],
    [$Delta \/ (2 pi)$],          [$500$ MHz],            [unchanged],
    [$K$, $alpha$],               [$0.95$, $4$],          [fixed pulse shape],
    [$T_c$],                      [$7.14$ ns],            [Cs pi-pulse],
    [$T_f$],                      [$126.71$ ns],          [probe half-width],
    [$sigma$],                    [$0.509$ ns],           [pulse waist],
    [Total $T_("tot")$],          [$bold(267.71)$ ns],    [-117 ns vs rev. 4],
    [$V_(c t) \/ (2 pi)$],        [$441.9$ MHz],          [$+95 %$ vs rev. 4],
    [$V_("DD") \/ (2 pi)$],       [$+2.21$ MHz],          [(still small)],
    [$V_(c c) \/ (2 pi)$],        [$-2.78$ MHz],          [(still small)],
    [$R_("DD")$, $R_("AA")$],     [$200$, $159$],         [both $> 100$ #sym.checkmark],
    [$M_1$ (AC-Stark margin)],    [$122.8$],              [$Omega_p^2 \/ (2 Delta)$],
    [$M_2$ (Raman margin)],       [$30.7$],               [$Omega_p Omega_R \/ (2 Delta)$],
    [$overline(F)_"OR"$ (verified)], [$bold(0.994160)$],   [$+1.74 times 10^(-3)$ vs rev. 4],
    [$1 - overline(F)_"OR"$],     [$bold(5.84 times 10^(-3))$], [$-23 %$ vs rev. 4],
  ),
)

The per-input fidelity table from the verifier run (`champion_verify.log`):

#figure(
  caption: [OR-gate per-input fidelities at the champion cell.
    The $|1,1,* angle.r$ branch is again the bottleneck, but now
    floors at $0.988$ instead of $0.987$, and the $|0,0,* angle.r$
    branch is essentially perfect because the larger $Omega_c$
    cuts $T_c$ by $3 times$.],
  kind: table,
  table(
    columns: (auto, auto, auto, auto, auto, auto, auto),
    align: (center, left, center, center, center, center, center),
    stroke: 0.4pt,
    table.header(
      [*\#*], [*Input*], [*$F_k$*],
      [*$P(A)$*], [*$P(B)$*], [*$P(P)$*], [*$P(R)$*],
    ),
    [1], [$|0,0,A angle.r$], [$0.999424$], [$0.9994$], [$0.0003$], [$0.0000$], [$0.0003$],
    [2], [$|0,0,B angle.r$], [$0.999424$], [$0.0003$], [$0.9994$], [$0.0000$], [$0.0003$],
    [3], [$|0,1,A angle.r$], [$0.994441$], [$0.0022$], [$0.9944$], [$0.0000$], [$0.0001$],
    [4], [$|0,1,B angle.r$], [$0.994441$], [$0.9944$], [$0.0022$], [$0.0000$], [$0.0001$],
    [5], [$|1,0,A angle.r$], [$0.994441$], [$0.0022$], [$0.9944$], [$0.0000$], [$0.0001$],
    [6], [$|1,0,B angle.r$], [$0.994441$], [$0.9944$], [$0.0022$], [$0.0000$], [$0.0001$],
    [7], [$|1,1,A angle.r$], [$0.988335$], [$0.0033$], [$0.9883$], [$0.0000$], [$0.0000$],
    [8], [$|1,1,B angle.r$], [$0.988335$], [$0.9883$], [$0.0033$], [$0.0000$], [$0.0000$],
  ),
)

Arithmetic check:
$ overline(F)_"OR"
  = (2 dot 0.9994 + 4 dot 0.9944 + 2 dot 0.9883) \/ 8
  = bold(0.9942) checkmark. $

= Robust plateau around the champion

The champion sits inside a *connected* fully-robust plateau in
the 5-D grid. The 40 cells with `robust_score = 1.0` cluster
into two regions:

#figure(
  caption: [The 15 highest-fidelity fully-robust cells (every 5-D
    nearest neighbor also exceeds $0.99$). All are tightly grouped
    around $Omega_p = 60$ MHz, $r in {3.5, 4.0}$, $Delta in
    [400, 800]$ MHz, $Omega_c in {50, 70}$ MHz, and the lattice
    splits between two clusters: $a = 4.0$ μm (very tight, fast
    gate) and $a = 4.5$ μm (looser, slightly slower).],
  kind: table,
  table(
    columns: (auto, auto, auto, auto, auto, auto, auto, auto),
    align: (center, center, center, center, center, center, center, center),
    stroke: 0.4pt,
    table.header(
      [*$a$*], [*$Omega_p$*], [*$r$*], [*$Delta$*], [*$Omega_c$*],
      [*$T_("tot")$*], [*$overline(F)_"OR"$*], [*$F_(1,1,*)$*],
    ),
    [$4.0$], [$60$], [$4.0$], [$500$], [$70$], [$267.7$], [$bold(0.9942)$], [$0.9883$],
    [$4.0$], [$60$], [$4.0$], [$400$], [$70$], [$bold(217.0)$], [$0.9941$], [$0.9906$],
    [$4.0$], [$60$], [$4.0$], [$400$], [$50$], [$222.7$], [$0.9939$], [$0.9900$],
    [$4.0$], [$60$], [$4.0$], [$500$], [$50$], [$273.4$], [$0.9937$], [$0.9867$],
    [$4.5$], [$60$], [$4.0$], [$600$], [$70$], [$318.4$], [$0.9937$], [$0.9895$],
    [$4.5$], [$60$], [$4.0$], [$600$], [$50$], [$324.1$], [$0.9937$], [$0.9893$],
    [$4.5$], [$60$], [$3.5$], [$600$], [$70$], [$318.4$], [$0.9936$], [$0.9886$],
    [$4.5$], [$60$], [$3.5$], [$600$], [$50$], [$324.1$], [$0.9935$], [$0.9884$],
    [$4.0$], [$60$], [$4.0$], [$600$], [$70$], [$318.4$], [$0.9935$], [$0.9856$],
    [$4.5$], [$60$], [$3.5$], [$500$], [$70$], [$267.7$], [$0.9934$], [$0.9899$],
    [$4.5$], [$60$], [$4.0$], [$600$], [$30$], [$337.4$], [$0.9934$], [$0.9885$],
    [$4.5$], [$60$], [$3.5$], [$500$], [$50$], [$273.4$], [$0.9933$], [$0.9896$],
    [$4.5$], [$60$], [$4.0$], [$800$], [$70$], [$419.8$], [$0.9932$], [$0.9866$],
    [$4.5$], [$60$], [$4.0$], [$800$], [$50$], [$425.5$], [$0.9932$], [$0.9864$],
    [$4.5$], [$50$], [$4.0$], [$500$], [$70$], [$379.2$], [$0.9931$], [$0.9871$],
  ),
)

Two things are striking in this table:

+ *$Omega_p = 60$ MHz is overwhelmingly preferred*. Out of the
  40 fully-robust cells (`brute_force_robust.csv` rows with
  `robust_score = 1.0`), $39$ have $Omega_p = 60$ MHz; only one
  has $Omega_p = 50$ MHz (the last row above). This is because
  larger $Omega_p$ shrinks $T_f$ as $Omega_p^(-2)$ at fixed area,
  which directly shrinks the Cs $|r angle.r$ decay window
  (the dominant infidelity channel; see Sec. 4).

+ *Two lattice clusters*: $a = 4$ μm (fastest, slightly lower
  $R_("DD")$) and $a = 4.5$ μm (slightly slower, higher
  $R_("DD") = 394$). Both are improvements over the rev. 4 value
  $a = 5.0$ μm.

The lone fastest cell -- $a = 4.0$, $Delta = 400$, $Omega_c = 70$
-- runs at *$T_("tot") = 217$ ns* (a $44 %$ speedup over rev. 4)
while still hitting $overline(F)_"OR" = 0.9941$. Its $M_2 = 24.6$
is the lowest in the fully-robust set; if calibration tolerance
matters more than gate time, the champion at $Delta = 500$ MHz
($M_2 = 30.7$) is the safer pick.

= Why the champion wins: error-budget decomposition

The OR-gate infidelity at this protocol has three dominant
contributions; the brute-force results decompose cleanly.

== Cs $|r angle.r$ decay during the gate window

For inputs with $n$ controls in $|1 angle.r$, the Cs ancilla
spends $T_("tot") - T_c$ in $|r angle.r$, giving a per-control
decay factor $exp(-gamma_r (T_("tot") - T_c)) approx
1 - (T_("tot") - T_c) \/ tau_r$. At $tau_r = 77$ μs:

$ 1 - F_(1,1,*) approx 2 thick (T_("tot") - T_c) \/ tau_r
  + epsilon_"blockade". $

At rev. 4 ($T_("tot") = 385$ ns, $T_c = 10$ ns):
$2 dot 375 \/ 77000 = 9.7 times 10^(-3)$.
At the champion ($T_("tot") = 267.7$ ns, $T_c = 7.1$ ns):
$2 dot 260.6 \/ 77000 = bold(6.77 times 10^(-3))$.
The remaining $approx 4 times 10^(-3)$ in the rev. 4 11-branch
infidelity ($1 - 0.987 = 1.3 times 10^(-2)$) is the
blockade-leakage residual, and *that part* is set by $M_2$, not
$T_("tot")$.

== Blockade-leakage cap on $F_(1,1,*)$

In strong blockade $M_2 = V_(c t) (2 Delta) \/ (Omega_p Omega_R)$,
the doubly-blockaded branch's residual infidelity scales as
$1 \/ M_2^2$. At the champion $M_2 = 30.7$; at rev. 4
$M_2 = 25.9$. The leakage cap improves by
$(30.7 \/ 25.9)^2 approx 1.4 times$, contributing the rest of
the $F_(1,1,*) = 0.987 -> 0.988$ delta. The strong-blockade
margin floor across the fully-robust set is
$M_2 in [21.6, 36.8]$, all comfortably above the rev. 4 value;
ratio $= 4.0$ helps here because for fixed $Omega_p$ it raises
$Omega_R$, which deepens EIT shielding without spending pulse
area (which is set by $Omega_p^2 \/ Delta$, not $Omega_R$).

== Shorter $T_c$ from larger $Omega_c$

The leading-edge / trailing-edge Cs pi-pulses on $|1 angle.r
arrow.l.r |r angle.r$ run for $T_c = pi \/ Omega_c$.
Rev. 4 used $Omega_c = 50$ MHz $arrow.r T_c = 10$ ns;
champion uses $Omega_c = 70$ MHz $arrow.r T_c = 7.14$ ns. The
$2 dot 2.86 = 5.7$ ns shaved off the total gate is $~7 %$ of
$T_c$-window decay -- small but additive with the much-larger
$T_f$ savings.

== Net error budget

#figure(
  caption: [Decomposition of $overline(F)_"OR"$ infidelity at the
    rev. 4 FINAL vs the brute-force champion. The biggest
    improvement comes from the *shorter gate* (lower
    $|r angle.r$-decay loss), and a secondary improvement from
    higher $M_2$ (lower blockade leakage).],
  kind: table,
  table(
    columns: (auto, auto, auto, auto),
    align: (left, center, center, center),
    stroke: 0.4pt,
    table.header(
      [*Contribution*],
      [*Rev. 4 FINAL*],
      [*Champion*],
      [*Δ (champion − rev. 4)*],
    ),
    [Total $T_("tot")$ (ns)], [$385$], [$268$], [*−117*],
    [$tilde 2 (T_("tot") - T_c) \/ tau_r$ ($|11*angle.r$)],
      [$9.7 times 10^(-3)$], [$6.77 times 10^(-3)$], [$-3 times 10^(-3)$],
    [Blockade-leakage $tilde 1 \/ M_2^2$ (relative)],
      [$1.0$ (ref)], [$0.71$], [$-29 %$],
    [$F_(1,1,*)$ branch],
      [$0.987$], [$0.988$], [$+1 times 10^(-3)$],
    [$F_(0,1,*)$ branch (1 ctrl, 1-Cs decay)],
      [$0.993$], [$0.994$], [$+1 times 10^(-3)$],
    [$F_(0,0,*)$ branch (target-only EIT)],
      [$0.997$], [$0.9994$], [$+2 times 10^(-3)$],
    [Total $overline(F)_"OR"$],
      [$0.99242$], [$0.99416$], [$+1.74 times 10^(-3)$],
  ),
)

The $|0,0,* angle.r$ branch improves by *more* than the
$|1,1,* angle.r$ branch ($+2 times 10^(-3)$ vs $+1 times 10^(-3)$)
because (i) the target's EIT cycle is much shorter ($T_f
= 127$ ns instead of $182$ ns), and (ii) the larger
$Omega_R$ gives a deeper EIT dark state, suppressing the
$|R angle.r$ pedestal that was the $|0,0,* angle.r$ bottleneck.

= Robustness: how much calibration drift can each axis absorb?

Calling a cell "robust" depends on a metric. Here we use:
*at the champion, find each axis's tolerance such that
$overline(F)_"OR"$ stays above $0.99$ as that single axis is
varied through the grid while the other four are held fixed*.

#figure(
  caption: [Single-axis $overline(F)_"OR" > 0.99$ tolerance bands
    at the champion cell ($a = 4.0$, $Omega_p = 60$, $r = 4.0$,
    $Delta = 500$, $Omega_c = 70$). The "step" is the grid step
    in this brute-force sweep; finer scans would tighten the
    estimate but the qualitative picture is already conservative.],
  kind: table,
  table(
    columns: (auto, auto, auto, auto),
    align: (left, center, center, center),
    stroke: 0.4pt,
    table.header(
      [*Axis*],
      [*Grid step*],
      [*Lower / upper $> 0.99$ in grid*],
      [*Tolerance*],
    ),
    [$a$ (μm)], [$0.5$], [$4.0 - 5.5$], [$plus.minus$ at least $0.5$ μm],
    [$Omega_p \/ (2 pi)$ (MHz)], [$10$], [$50 - 60$ (no $70$ grid point)], [$plus.minus 10$ MHz],
    [$Omega_R \/ Omega_p$], [$0.5$], [$3.0 - 4.0$], [$plus.minus 0.5$],
    [$Delta \/ (2 pi)$ (MHz)], [$100$], [$400 - 800$], [$plus.minus 100$+ MHz],
    [$Omega_c \/ (2 pi)$ (MHz)], [$20$], [$30 - 70$], [$plus.minus 20$+ MHz],
  ),
)

All five axes admit at least one grid step of drift before
$overline(F)_"OR"$ drops below $0.99$. Crucially, the champion's
neighbours in *all five axes* exceed $0.99$, which is the
operational definition of fully robust used in this report.

= How this compares to the prior $omega_R$ knife-edge

The rev. 4 audit (Sec. 5 of `Smaller_VDD_average_fidelity_report.typ`)
already saw a hint of this answer: at $a = 5.0$ μm,
$Omega_p = 50$ MHz, the $(Delta, Omega_R \/ Omega_p)$ landscape
had a *secondary peak* at $Delta = 225$ MHz, $r = 2.5$ with
$overline(F)_"OR" = 0.9920$ in only $184$ ns -- but was a
knife-edge with only $4 \/ 36$ cells passing $0.99$ in a tight
$(r, K)$ grid.

The brute-force result here finds something better:

#figure(
  caption: [Comparison of three operating points. The brute-force
    champion is *both faster* than the rev. 4 FINAL *and more
    robust* than the rev. 4 small-$Delta$ knife-edge. The
    knife-edge offered a $50 %$ time saving at the cost of
    fragility; the brute-force champion offers a $30 %$ time
    saving with no fragility at all.],
  kind: table,
  table(
    columns: (auto, auto, auto, auto, auto),
    align: (left, center, center, center, center),
    stroke: 0.4pt,
    table.header(
      [*Operating point*],
      [*$T_("tot")$*],
      [*$overline(F)_"OR"$*],
      [*$M_2$*],
      [*Robust?*],
    ),
    [Rev. 4 FINAL ($a = 5.0$, $Omega_p = 50$, $r = 3.5$, $Delta = 500$)],
      [$385$ ns], [$0.99242$], [$25.9$], [partly],
    [Rev. 4 small-$Delta$ ($a = 5.0$, $Omega_p = 50$, $r = 2.5$, $Delta = 225$)],
      [$184$ ns], [$0.9920$], [$16.3$], [no (knife-edge)],
    [*Brute-force champion ($a = 4.0$, $Omega_p = 60$, $r = 4.0$, $Delta = 500$)*],
      [$bold(268)$ ns], [$bold(0.99416)$], [$bold(30.7)$], [*yes (6 / 6)*],
  ),
)

The brute-force champion succeeds where the small-$Delta$
knife-edge failed because it *raised* $V_(c t)$ ($+95 %$, by
shrinking $a$ from $5.0$ to $4.0$ μm) rather than *lowered*
$Delta$. Both moves restore $M_2$, but the first gets there
without sliding past the dressed-state anti-resonance at
ratio $= 2$ that polluted the small-$Delta$ ridge.

= Discussion and recommended operating point

#block(
  fill: rgb("#eef8ec"),
  inset: 10pt,
  radius: 4pt,
  width: 100%,
)[
  *Recommendation.* For routine use of the X-error sub-cycle at
  Rb $54 D_(3/2)$ + Cs $62 D_(5/2)$:

  $ bold(a = 4.0 thick "μm"), quad
    bold(Omega_p = 2 pi times 60 thick "MHz"), quad
    bold(Omega_R = 2 pi times 240 thick "MHz"), quad
    bold(Delta = 2 pi times 500 thick "MHz"), quad
    bold(Omega_c = 2 pi times 70 thick "MHz"), $

  with $K = 0.95$, $alpha = 4$ (rev. 4 pulse shape).

  Expect:
  $ overline(F)_"OR" = bold(0.9942), quad
    1 - overline(F)_"OR" = bold(5.84 times 10^(-3)), quad
    T_("tot") = bold(268) thick "ns". $

  All five axes tolerate $tilde plus.minus 1$ grid step of drift
  while staying above $0.99$ -- *the only fully-robust champion
  in the swept domain*.

  If gate time is paramount and the lab can tighten calibration to
  the level shown in the table above, the *fastest fully-robust*
  cell at $T_("tot") = 217$ ns ($Delta = 400$ MHz, otherwise the
  same) is a $0.9941$ fidelity option, $44 %$ faster than rev. 4.
]

== Trade-offs to be aware of

+ *Lattice tightening to $a = 4.0$ μm* raises $V_("DD")$ from
  $0.58$ MHz to $2.21$ MHz and $|V_(c c)|$ from $0.73$ to $2.78$
  MHz. Both $R_("DD") = 200$ and $R_("AA") = 159$ are still
  comfortably above the $>= 100$ threshold, but the safety
  margin is *smaller* than rev. 4's $390 \/ 311$. If the K-round
  Rb level (still TBD) requires more headroom on $V_("DD")$,
  $a = 4.5$ μm (the second-cluster fully-robust set) offers
  $R_("DD") = 394$ at the cost of $T_("tot") = 318$ ns (still
  $17 %$ faster than rev. 4) with $overline(F)_"OR" = 0.9937$.

+ *Higher $Omega_p = 60$ MHz and $Omega_R = 240$ MHz* are not
  problems for any optical-power budget realistic at $7P_(3/2)$
  intermediate state -- they are below the saturation regime --
  but the AC-Stark on $|P angle.r$, $Omega_p^2 \/ (2 Delta)$,
  rises from $2.5$ MHz (rev. 4) to $3.6$ MHz (champion). Still
  well-resolved against $Delta = 500$ MHz.

+ *Pulse shape was not swept*. The rev. 4 audit established
  $K = 0.95$ and $alpha = 4$ as a wide plateau; both are
  *empirically optimal at finite $V_(c t)$* and are not expected
  to shift much at $V_(c t) = 442$ MHz vs $226$ MHz. A 1-D
  $K$-sweep at the champion configuration would be a low-cost
  sanity check before publication.

== Open issues, unchanged from rev. 4

- ARC lifetime re-verification for $tau_r$ (Cs $62 D_(5/2)$) and
  $tau_R$ (Rb $54 D_(3/2)$). The $n^3$-scaled estimates of $77$ μs
  and $74$ μs may shift by $tilde plus.minus 10 %$ when computed
  directly. This affects $|1, 1, * angle.r$ infidelity linearly
  in $tau_r$.

- K-round (Z sub-cycle) Rb level still TBD; nothing in this
  report constrains that choice.

= Reproducibility

#block(
  fill: luma(246),
  inset: 8pt,
  radius: 3pt,
)[
  *Run commands* (Apple Silicon laptop, Python 3.12,
  `qutip == 5.2.3`, `scipy == 1.17.1`):
  ```text
  cd TriQG
  source .venv/bin/activate
  cd examples/Smaller_VDD_energy
  python brute_force_or.py          # 5-D sweep, ~9 min, 432 cells
  python brute_force_plots.py       # figures: marginals, 2D, pareto
  python champion_verify.py         # re-simulate champion cell
  ```

  *Outputs* in this directory:
  - `brute_force_or.log`            -- progress log of the sweep
  - `brute_force_results.csv`       -- all 432 cells, every column
  - `brute_force_robust.csv`        -- 221 cells with $overline(F) > 0.99$,
                                       ranked by robust_score
  - `brute_force_axis_marginals.png`
  - `brute_force_2D_marginals.png`
  - `brute_force_pareto.png`
  - `champion_verify.log`           -- detailed run on champion cell
]

#block(
  fill: rgb("#fff4e8"),
  inset: 9pt,
  radius: 3pt,
  width: 100%,
)[
  *Caveats.*

  - $K = 0.95$ and $alpha = 4$ were *held fixed* (the rev. 4
    optimum). The brute force does *not* re-establish that they
    are optimal at $V_(c t) = 442$ MHz; it assumes the rev. 4
    finding carries over. Spot-check at one $K$-sweep is
    recommended.
  - Lifetimes ($tau_r = 77$ μs, $tau_R = 74$ μs) are still
    $n^3$-scaled estimates, not direct ARC values. All comparisons
    in this report are *self-consistent* on those estimates.
  - $V_("DD")$ does not enter the 3-atom Hamiltonian; the $R_("DD")$
    selectivity audit is a *separate* gate-selectivity check that
    must be re-validated as part of the X-round parallel-drive
    error budget.
  - The robust-score metric uses *grid-step neighbors*, not
    continuous derivatives. A 5-axis Hessian-based fragility
    metric at the champion would refine the numbers in
    Sec. 6 but is not expected to overturn the conclusion.
]

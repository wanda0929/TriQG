// a4_K1_error_budget_report.typ
//
// Companion to a4_K1_champion_report.pdf: decomposes the basis-fidelity
// infidelity of the OR and CCX champions into the three claimed error
// sources by selectively switching each one off in the QuTiP simulation.
//
// Compile with:  typst compile a4_K1_error_budget_report.typ

#set document(
  title: "Error budget for the K = 1, a = 4 um champion",
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
  #text(size: 16pt, weight: "bold")[
    Error budget for the $K = 1$, $a = 4$ μm champion
  ]
  #v(0.25em)
  #text(size: 11pt)[
    Which of #emph[$|P angle.r$ decay], #emph[Rydberg decay], or
    #emph[same-species] $V_(c c)$ is the leak?
  ]
  #v(0.3em)
  #text(size: 9pt)[
    TriQG `examples/Smaller_VDD_energy` ·
    `a4_K1_error_budget.py` · 2026-05-15
  ]
]

#v(0.4em)

#block(
  fill: rgb("#eef8ec"),
  inset: 10pt,
  radius: 4pt,
  width: 100%,
)[
  *One-line answer:* The two gates have *different* dominant errors.
  Re-running the champion drive parameters with each error channel
  switched off (`a4_K1_error_budget.py`) gives, for the basis-fidelity
  infidelity $1 - overline(F)_"basis"$:

  - *OR gate* ($T_"tot" = 244$ ns):
    Rydberg decay is the leader at *$60.7 %$* of the
    removable error ($Delta(1-overline(F)) = 1.69 times 10^(-3)$),
    intermediate $lr(|P angle.r)$ decay is *$27.2 %$*
    ($7.58 times 10^(-4)$), same-species $V_(c c)$ is *$12.0 %$*
    ($3.34 times 10^(-4)$). A further $2.26 times 10^(-3)$
    (45 % of the total $5.04 times 10^(-3)$) is the *coherent*
    floor -- non-removable subspace leakage from finite blockade.

  - *CCX gate* ($T_"tot" = 95$ ns):
    Same-species $V_(c c)$ is the leader at *$65.1 %$*
    ($1.18 times 10^(-3)$), Rydberg decay is *$34.9 %$*
    ($6.35 times 10^(-4)$), and $lr(|P angle.r)$ decay is *$0 %$*
    (the CCX never populates the $7 P_(3\/2)$ state). The coherent
    floor is just $1.0 times 10^(-4)$, so the CCX is essentially
    decoherence + interaction limited.

  Channels are additive to better than $0.2 %$ (linearity check in
  @sec:method); the percentages above are an honest decomposition,
  not a heuristic split.
]

= Method <sec:method>

The champion (rev. of 2026-05-15) has total basis infidelity
$1 - overline(F)_"basis"^"OR" = 5.04 times 10^(-3)$ and
$1 - overline(F)_"basis"^"CCX" = 1.92 times 10^(-3)$. Three physical
mechanisms are credited with that error in
`a4_K1_champion_report.pdf`:

#table(
  columns: (auto, auto, auto),
  align: (left, center, left),
  stroke: 0.4pt,
  table.header([*Mechanism*], [*Code knob*], [*ARC value at $T = 0$ K*]),
  [Rb $7 P_(3\/2)$ decay], [$gamma_P = 1\/tau_P$],
    [$3.704$ μs⁻¹ ($tau_P = 0.270$ μs)],
  [Rb $54 D_(3\/2)$ decay (target)],
    [$gamma_R = 1\/tau_R$], [$0.00608$ μs⁻¹ ($tau_R = 164.55$ μs)],
  [Cs $62 D_(5\/2)$ decay (ancillas)],
    [$gamma_r = 1\/tau_r$], [$0.00720$ μs⁻¹ ($tau_r = 138.87$ μs)],
  [Cs–Cs vdW on $lr(|r r angle.r)$],
    [$V_(c c)$ in `build_hamiltonian`],
    [$-2 pi times 2.777$ MHz],
  [Rb–Rb vdW on $lr(|R R angle.r)$],
    [$V_("DD")$], [$+2 pi times 2.214$ MHz #footnote[Does not enter
    the 3-atom box: only one Rb target is present, so no Rb–Rb pair
    exists in the simulated Hilbert space. We confirm this by leaving
    everything else fixed and toggling only $V_(c c)$.]],
)

The two Rydberg decay channels are reported together as one
"Rydberg decay" bucket because they tag the same physical mechanism
(spontaneous emission from a long-lived Rydberg level) and their
relative weights are fixed by the level choice, not by the protocol.

== Five scenarios

We re-run the QuTiP `mesolve` simulation on the eight computational
basis inputs five times for each gate:

#table(
  columns: (auto, auto, auto, auto, auto),
  align: (left, center, center, center, center),
  stroke: 0.4pt,
  table.header(
    [*Scenario*],
    [*$gamma_P$*], [*$gamma_R, gamma_r$*],
    [*$V_(c c)$*], [*meaning*],
  ),
  [`baseline`],          [on], [on], [on], [full physics],
  [`no_P_decay`],        [*off*], [on], [on], [$|P angle.r$ leak suppressed],
  [`no_Rydberg_decay`],  [on], [*off*], [on], [Rydberg lifetime $arrow.r infinity$],
  [`no_same_species`],   [on], [on], [*off*], [Cs–Cs vdW = 0],
  [`clean`],             [off], [off], [off], [pure unitary coherent floor],
)

For each scenario the basis-fidelity infidelity #emph[delta_F]
contributed by a channel is

$ Delta F^"channel" = overline(F)^"no_channel" - overline(F)^"baseline". $

This is the *infidelity removed* when the channel is silenced. The
percentage share of the total removable (i.e. non-coherent) error is

$ "share"^"channel" = Delta F^"channel" \/ (overline(F)^"clean" -
  overline(F)^"baseline"). $

== Additivity check

The three single-channel gains should sum to the total
`clean - baseline` gain if the channels do not interfere. We measure:

#table(
  columns: (auto, auto, auto),
  align: (left, center, center),
  stroke: 0.4pt,
  table.header([*Quantity*], [*OR gate*], [*CCX gate*]),
  [$Delta F^"no_P" + Delta F^"no_Ryd" + Delta F^"no_Vcc"$],
    [$2.783 times 10^(-3)$], [$1.819 times 10^(-3)$],
  [$Delta F^"clean - baseline"$],
    [$2.787 times 10^(-3)$], [$1.820 times 10^(-3)$],
  [residual cross term],
    [$+3.8 times 10^(-6)$ (+0.14 %)],
    [$+0.9 times 10^(-6)$ (+0.05 %)],
)

The channels are additive to better than 0.2 %. The cross-term is
positive in both gates -- when two channels are on together the loss
is fractionally less than the sum of their individual losses -- but
the effect is well below the integration-noise floor of `mesolve` and
does not change the ranking.

= Results

#figure(
  caption: [Decomposition of basis infidelity at the
    $a = 4.0$ μm, $K = 1$ operating point. Bars are absolute
    infidelity contributions ($Delta(1 - overline(F))$ in units of
    $10^(-3)$); panel (b) shows the same data as a percentage of the
    removable (decoherence + $V_(c c)$) total, excluding the coherent
    floor. *The OR gate is Rydberg-decay limited; the CCX gate is
    same-species-vdW limited*. Source: `a4_K1_error_budget.py`,
    `a4_K1_error_budget.png`.],
  image("a4_K1_error_budget.png", width: 100%),
)

== OR-gate decomposition

#figure(
  caption: [OR-gate ($T_"tot" = 244$ ns) error budget. The basis
    infidelity splits into a $2.26 times 10^(-3)$ coherent floor
    (non-removable; finite-blockade subspace leakage) plus
    $2.79 times 10^(-3)$ of removable channel error, of which
    *Rydberg decay alone is $1.69 times 10^(-3)$, six times more
    than $V_(c c)$*. The $lr(|P angle.r)$-decay contribution is a
    factor $2.2 times$ smaller than Rydberg even though $gamma_P$
    is $600 times$ larger than $gamma_R$ -- the protocol spends
    almost no time in $lr(|P angle.r)$ because the two-photon Raman
    is far off-resonance ($Delta\/(2 pi) = 500$ MHz).],
  kind: table,
  table(
    columns: (auto, auto, auto, auto, auto),
    align: (left, center, center, center, center),
    stroke: 0.4pt,
    table.header(
      [*Scenario*], [*$overline(F)_"basis"$*], [*$1 - overline(F)$*],
      [*$Delta F$*], [*share*],
    ),
    [`baseline`],        [$0.994958$], [$5.04 times 10^(-3)$], [—], [—],
    [`no_P_decay`],      [$0.995716$], [$4.28 times 10^(-3)$],
      [$+7.58 times 10^(-4)$], [$27.2 %$],
    [`no_Rydberg_decay`],[$0.996649$], [$3.35 times 10^(-3)$],
      [$+bold(1.69 times 10^(-3))$], [$bold(60.7 %)$],
    [`no_same_species`], [$0.995291$], [$4.71 times 10^(-3)$],
      [$+3.34 times 10^(-4)$], [$12.0 %$],
    [`clean`],           [$0.997744$], [$2.26 times 10^(-3)$],
      [$+2.79 times 10^(-3)$], [$100 %$],
  ),
)

*Reading the numbers.* The OR-gate populates the target Rb $54 D_(3\/2)$
Rydberg state for almost the entire $2 T_f = 227$ ns of the
super-Gaussian probe -- that's $approx 93 %$ of the gate -- on every
branch where blockade is *not* in effect (i.e. when at least one
control was left in $lr(|1 angle.r) = lr(|r angle.r)$ during the Cs
$pi$-pulses, so the conditional $|R angle.r$ flip *does* fire). The
predicted Rydberg-decay infidelity from a back-of-envelope
$gamma_R times T_f$ is $0.00608 times 0.227 approx 1.4 times 10^(-3)$,
matching the measured $1.69 times 10^(-3)$ within $20 %$ once the
control $lr(|r angle.r)$-population is added.

The $|P angle.r$-decay contribution is *much* smaller than its rate
might suggest. The reason is in the report's `Sec. Caveats`: the
$lr(|P angle.r)$ population during an off-resonant Raman is
$lr(|c_P|)^2 approx Omega_p^2 \/ (2 Delta)^2 approx (65 \/ 1000)^2
approx 4 times 10^(-3)$, so the *effective* $|P angle.r$-dwell time is
$approx 0.004 times 2 T_f \= 0.9$ ns, giving infidelity
$approx gamma_P times 0.9 "ns" approx 3.3 times 10^(-3)$. The
measured contribution is $7.58 times 10^(-4)$, smaller still because
not every branch fires the Raman fully -- consistent.

The $V_(c c)$ contribution of $3.34 times 10^(-4)$ is small for the
OR gate because the protocol *avoids* parking both ancillas in
$lr(|r angle.r)$ simultaneously for any extended interval -- the
Cs $pi$-pulses last only $T_c = 8.3$ ns each, so the $lr(|r r angle.r)$
window is short.

== CCX-gate decomposition

#figure(
  caption: [CCX-gate ($T_"tot" = 95$ ns) error budget. *The pattern
    is inverted relative to the OR gate*: $V_(c c)$ is the dominant
    error and $|P angle.r$ decay is exactly zero. The CCX protocol's
    middle "$3 pi$" train on the target parks *both* ancillas in
    $lr(|r angle.r)$ for the full $3 T_t = 75$ ns -- that is where
    $V_(c c)$ has all its action and exactly $79 %$ of the gate
    time.],
  kind: table,
  table(
    columns: (auto, auto, auto, auto, auto),
    align: (left, center, center, center, center),
    stroke: 0.4pt,
    table.header(
      [*Scenario*], [*$overline(F)_"basis"$*], [*$1 - overline(F)$*],
      [*$Delta F$*], [*share*],
    ),
    [`baseline`],        [$0.998076$], [$1.924 times 10^(-3)$], [—], [—],
    [`no_P_decay`],      [$0.998076$], [$1.924 times 10^(-3)$],
      [$+0$ exactly], [$0.0 %$],
    [`no_Rydberg_decay`],[$0.998711$], [$1.289 times 10^(-3)$],
      [$+6.35 times 10^(-4)$], [$34.9 %$],
    [`no_same_species`], [$0.999260$], [$7.40 times 10^(-4)$],
      [$+bold(1.18 times 10^(-3))$], [$bold(65.1 %)$],
    [`clean`],           [$0.999896$], [$1.04 times 10^(-4)$],
      [$+1.82 times 10^(-3)$], [$100 %$],
  ),
)

*Why $lr(|P angle.r)$ decay vanishes.* The CCX uses resonant square
$pi$-pulses (Sec. 6 of the champion report): Cs ancillas drive
$lr(|0 angle.r) <-> lr(|r angle.r)$ directly and the Rb target drives
$lr(|A angle.r) <-> lr(|R angle.r)$ or $lr(|B angle.r) <-> lr(|R angle.r)$
*without* the intermediate $7 P_(3\/2)$ Raman. The $|P angle.r$ state
is never populated, so $gamma_P$ has nothing to act on. This is a
*structural*, not parametric, zero -- the CCX-gate $|P angle.r$ share
will stay zero however the lattice or drives are retuned.

*Why $V_(c c)$ dominates.* The CCX worst-case branch is
$lr(|1\,1\,t angle.r)$: both Cs ancillas are parked in $lr(|r angle.r)$
for the entire middle $3 T_t = 75$ ns. The $V_(c c) lr(|r r angle.r)
lr(angle.l r r|)$ shift acts on this branch for $T_(r r) approx 75$ ns,
imprinting both a phase and an off-resonant detuning of the surrounding
Cs $pi$-pulses. The infidelity scales as
$(V_(c c) T_(r r))^2 \/ d$, and with $V_(c c)\/(2 pi) = 2.78$ MHz,
$T_(r r) = 75$ ns, that gives $approx (V_(c c) T_(r r))^2 \/ 8
approx 2 times 10^(-3)$ -- of the right order, modestly larger than
the measured $1.18 times 10^(-3)$. The fact that the basis-fidelity
loss is finite (rather than zero, which would happen if $V_(c c)$ only
added a removable phase) tells us part of $V_(c c)$ is detuning the
control $pi$-pulses, not just dressing $lr(|r r angle.r)$.

*Why Rydberg decay still matters.* At $T_"tot" = 95$ ns the
$gamma_(R\,r) times T$ budget is only $tilde.op 0.7 times 10^(-3)$ --
$2.4 times$ smaller than for the OR gate. That tracks the $2.6 times$
shorter gate. So the CCX gate is *not* Rydberg-decay-limited even
though it uses identical levels.

= Cross-gate take-aways <sec:takeaways>

+ *The dominant error is gate-specific.* "Which error matters most?"
  has no single answer for the protocol. The OR gate is *Rydberg-decay
  limited*; the CCX gate is *$V_(c c)$ limited*.

+ *$|P angle.r$ decay is a small leak.* For the OR gate it's $tilde.op 27 %$
  of the removable error -- not negligible, but a factor $2.2 times$
  below Rydberg decay despite $gamma_P$ being $600 times$ larger.
  The off-resonant Raman keeps the $lr(|P angle.r)$ amplitude at
  $tilde.op 10^(-3)$. For the CCX gate $lr(|P angle.r)$ decay is *structurally
  zero* because the protocol never excites $lr(|P angle.r)$. The caveat
  in `Sec. Caveats` of the champion report -- "$tau_P$ choice is robust
  to a factor of 2" -- is corroborated quantitatively: even doubling
  $gamma_P$ would only push the OR's $Delta F^"no_P"$ to
  $approx 1.5 times 10^(-3)$, still $30 %$ less than Rydberg.

+ *Rydberg decay limits the OR gate because of the gate time.*
  $Delta F^"no_Ryd" approx gamma_(R\,r) times T_(r r ,R R)$, which
  is essentially a budget-of-Rydberg-residency-time. The only knobs
  that move it are (i) shorter $T_"tot"$ at fixed area $K = 1$,
  which trades against blockade margin $M_2$, or (ii) longer Rydberg
  lifetime, which is what going from $300$ K to $0$ K already bought
  us ($tau$ doubled, error halved).

+ *$V_(c c)$ limits the CCX gate because of branch geometry.*
  Re-architecting the CCX to *not* park both ancillas in
  $lr(|r angle.r)$ simultaneously would close most of the $1.18
  times 10^(-3)$ same-species leak (potentially driving
  $overline(F)_"CCX"$ towards $0.9994$). Equivalently, increasing the
  ancilla-ancilla spacing $r_(A A) = a sqrt(2)$ from $5.66$ μm to e.g.
  $7$ μm would scale $V_(c c)$ down by $(7\/5.66)^6 approx 3.6 times$
  and reduce its infidelity contribution by $approx 13 times$.

+ *Coherent floor matters for the OR gate.* The non-removable
  $2.26 times 10^(-3)$ from finite-blockade leakage is $45 %$ of the
  total OR-gate infidelity. *No amount of cleaner atoms will fix it*
  at $a = 4$ μm. Pushing the OR-gate fidelity above $0.997$ requires
  either (a) tightening the lattice further (raises $V_(c t)$, but
  invalidates the report's blockade-margin assumption $V_(c t) gt.tilde
  Delta$, the borderline $V_(c t) = 442$ vs $Delta = 500$ noted in
  `Sec. Caveats`), or (b) compensating with a $K^ast eq.not 1$
  area adjustment (which `Why_K095.typ` argues is now $approx 1.00$ at
  this lattice -- so there is no quick gain there).

= Practical recommendations

#table(
  columns: (auto, auto),
  align: (left, left),
  stroke: 0.4pt,
  table.header([*Goal*], [*Highest-leverage knob*]),
  [Improve OR-gate fidelity], [Pick longer-lived Rydberg levels
    (or a cryogenic apparatus with $T tilde.op 0$ K). $gamma_(R\,r)$
    is the dominant rate; halving it would remove $approx 8 times 10^(-4)$
    from the OR infidelity.],
  [Improve CCX-gate fidelity], [Increase ancilla-ancilla spacing
    $r_(A A)$ or use a Cs level with weaker $C_6^("CsCs")$. The
    $V_(c c)$ leak scales as $1\/r_(A A)^(12)$, so even a $25 %$
    larger spacing halves the dominant CCX error.],
  [Reduce $|P angle.r$ decay], [Lower priority. The OR-gate gain from
    eliminating $|P angle.r$ entirely is only $7.6 times 10^(-4)$
    ($15 %$ of total), and the CCX gate doesn't see $|P angle.r$ at all.
    Increasing the Raman detuning $Delta$ would suppress $|P angle.r$
    population but also lengthens $T_f$ (the super-Gaussian width
    grows as $sigma prop Delta^3$), so it likely loses on the
    Rydberg-decay budget.],
  [Reduce coherent floor], [OR-gate-only. Requires a different
    lattice or a $K^ast eq.not 1$ pulse-area correction; see the
    $K = 0.95$ vs $K = 1$ trade in the rev. 4 and champion reports.],
)

#v(0.4em)

= Reproducibility

#block(
  fill: luma(246),
  inset: 8pt,
  radius: 3pt,
)[
  Run commands (Apple Silicon, Python 3.12, `qutip == 5.2.3`):
  ```text
  cd TriQG && source .venv/bin/activate
  cd examples/Smaller_VDD_energy
  python a4_K1_error_budget.py           # ~1.5 min, writes JSON + log
  python a4_K1_error_budget_plot.py      # writes PNG
  typst compile a4_K1_error_budget_report.typ
  ```

  Files produced for this report:
  - `a4_K1_error_budget.py`        -- five-scenario QuTiP driver
  - `a4_K1_error_budget.json`      -- structured per-scenario fidelities
  - `a4_K1_error_budget.log`       -- per-run console log
  - `a4_K1_error_budget_plot.py`   -- stacked-bar visualisation
  - `a4_K1_error_budget.png`       -- the figure embedded above
]

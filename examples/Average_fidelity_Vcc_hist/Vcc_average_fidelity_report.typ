// Option A report: Rb 66 D_5/2 + Cs 76 D_3/2 at a = 5.00 um,
// with ARC-computed lifetimes (replacing the paper's anchor values),
// Omega_R = 2.9 * Omega_p (user request, was 3.5),
// and Option 1 pulse cleanup (bare omega_p, sigma = 0.001771).
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
    Option A at $a = 5.00$ μm with realistic lifetimes
  ]

  #v(0.3em)
  #text(size: 11pt)[
    Rb $66 D_(5/2)$ + Cs $76 D_(3/2)$ Förster pair, ARC-computed
    BBR lifetimes, $Omega_R = 2.9 thick Omega_p$, bare
    $Omega_p$ with $sigma = 0.001771$ (Option 1 pulse)
  ]

  #v(0.4em)
  #text(size: 10pt)[
    TriQG `examples/Average_fidelity_Vcc_newenergy` · compiled 2026-04-11
  ]
]

#v(0.6em)

#block(
  fill: rgb("#ffece6"),
  inset: 10pt,
  radius: 4pt,
  width: 100%,
)[
  *TL;DR -- Option A is worse than the paper baseline, but not for
  the reason you would expect.* Four changes were stacked together
  in this run:

  - *(a) Rydberg level swap* Rb $69 D_(5/2) arrow 66 D_(5/2)$ and
    Cs $79 D_(5/2) arrow 76 D_(3/2)$, with $tilde(C)_3$ and $C_6$
    updated to the Ireland et al. 2024 Table I values for this pair.
  - *(b) Lifetimes* replaced by the ARC 3.10.2 `getStateLifetime`
    results at $T = 300$ K (see §3), which are *about half* of the
    values quoted in `main.tex` for the paper pair.
  - *(c)* $Omega_R slash Omega_p: 3.5 arrow 2.9$ (user request).
  - *(d) Option 1 pulse cleanup:* remove the ad-hoc $1.039975$
    multiplier on $Omega_p$ and rescale $sigma: 0.0014 arrow 0.001771$
    so the effective two-photon pulse area stays fixed at $pi/4$.
    The pulse shape shifts by at most $approx 4$ MHz in peak height
    and is slightly smoother at the edges.

  Result:
  - OR gate: $overline(F)_"OR" = 0.995584 arrow bold(0.991858)$,
    infidelity $4.42 times 10^(-3) arrow bold(8.14 times 10^(-3))$.
  - CCX gate: $overline(F)_"CCX" = 0.997273 arrow bold(0.996620)$,
    infidelity $2.73 times 10^(-3) arrow bold(3.38 times 10^(-3))$
    (CCX does not use $Omega_p$ / $Omega_R$ and was not re-run under
    change (d); its value is from the previous Option A run).

  *What caused the drop:* changes (b) and (c) together swamp the
  Option A level-swap gain. Decomposed:
  - (b) *ARC lifetimes* add $approx 10^(-3)$ infidelity per branch
    across both gates because $T_1$ is roughly halved.
  - (c) $bold(Omega_R slash Omega_p = 2.9)$ re-tunes the EIT dark
    state: the $|0,0,* angle.r$ branch -- which was flat at $0.9992$
    under the $3.5$ value -- drops to $0.9936$ and grows a
    $0.0030$ residual $|R angle.r$ population (15$times$ larger than
    before). This is the *biggest single contributor* to the new
    infidelity, and it is a pulse-parameter choice, not physics.
  - (d) *Option 1 pulse cleanup* gives a small positive effect:
    the smoother edges improve the adiabaticity of the EIT dark
    state, cutting the $|0,0,* angle.r$ residual $|R angle.r$
    population from $0.0038$ (fudged pulse) to $0.0030$
    ($approx -21%$), and lifting
    $overline(F)_"OR"$ by $+2.6 times 10^(-4)$ relative to the same
    Option A run with the fudged pulse.
  - *(a) Level swap alone* would have given
    $overline(F)_"OR" approx 0.9956$ (slightly better than the paper)
    and $overline(F)_"CCX" approx 0.9973$ (essentially tied with the
    paper), matching the earlier prediction.

  If you want a clean read of the Option A gain, re-run without (c):
  revert `omega_R_amp` to `3.5 * omega_p_amp` and leave everything
  else (including the Option 1 pulse cleanup) as is.
]

= Scope of this run

This run is a three-way parameter update to the Option A pair from
the previous conversation, executed in the `_newenergy` folder:

+ *Physically swap the Rydberg pair* from the paper's worked example
  (Rb $69 D_(5/2)$ + Cs $79 D_(5/2)$) to the Ireland, Pritchard,
  Shaffer, PRR 6, 013293 (2024) Table I row 1 candidate:
  Rb $66 D_(5/2)$ + Cs $76 D_(3/2)$. Both the $tilde(C)_3$ dipole--dipole
  coefficient and the $C_6^"CsCs"$ van der Waals coefficient are
  updated to the tabulated values for this new pair.

+ *Use ARC-computed lifetimes* for the two new Rydberg states at
  $T = 300$ K, replacing the earlier rough estimates we used two
  reports ago.

+ *Switch the pulse ratio* $Omega_R slash Omega_p$ from the paper's
  $3.5$ to $2.9$ (requested by the user).

+ *Option 1 pulse cleanup:* remove the ad-hoc $1.039975$ multiplier
  on $Omega_p$ and slightly enlarge the super-Gaussian width
  $sigma: 0.0014 arrow 0.001771$ so the effective two-photon pulse
  area $integral Omega_p^2 \/ (2 Delta) dif t$ stays fixed at its
  original design value of $pi/4$. The pulse shape is visually
  indistinguishable from the previous one (peak difference
  $lt.approx 4$ MHz) and the total gate time is unchanged at 320 ns,
  but the smoother edges slightly improve the EIT-dark-state
  adiabaticity on the no-blockade branch.

Everything else stays at the paper-consistent values (lattice
$a = 5.00$ μm, $Omega_c = Omega_p = 2 pi times 50$ MHz,
$Omega_(c c) = 2 pi times 100$ MHz, $Omega_t = 2 pi times 50$ MHz,
$Delta = 2 pi times 500$ MHz, timing windows unchanged).

The previous `newenergy` run (a = 5.25 μm, $V_(c c)$ overridden to
$+5$ MHz, old lifetimes, $Omega_R = 3.5 Omega_p$) is now historical
and will be referenced only for comparison. The scripts and this
report now represent Option A proper.

= ARC lifetime lookup

Lifetimes were computed with ARC 3.10.2:

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

ARC's `getStateLifetime` combines the radiative Einstein-A sum with
the Beterov et al. 2009 blackbody-radiation-induced depopulation
formula. Sources:

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

The Option A lifetimes at $T = 300$ K:

#figure(
  caption: [ARC-computed BBR-included lifetimes at $T = 300$ K for the
    Option A Rydberg states. $tau_"rad"$ is the radiative component
    (at $T = 0$); $tau_"BBR"$ is the BBR-induced depopulation time;
    $tau_"total"^(-1) = tau_"rad"^(-1) + tau_"BBR"^(-1)$. Values come
    straight from `atom.getStateLifetime(n, l, j, temperature=300,
    includeLevelsUpTo=n+30)`.],
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

== Discrepancy with the paper's quoted lifetimes

As a cross-check I computed the same numbers for the *paper's* pair
states and compared them to the values quoted in `main.tex`:

#figure(
  caption: [Cross-check of the paper's Rydberg lifetimes against ARC.
    Paper values from `main.tex` Sec. III.C; ARC values from the same
    `getStateLifetime` call at $T = 300$ K.],
  kind: table,
  table(
    columns: (auto, auto, auto, auto),
    align: (left, center, center, center),
    stroke: 0.4pt,
    table.header(
      [*State*], [*paper $tau$*], [*ARC $tau$*], [*ratio*],
    ),
    [Rb $|69 D_(5/2) angle.r$], [$260$ μs],  [$150.27$ μs], [$0.578$],
    [Cs $|79 D_(5/2) angle.r$], [$340$ μs],  [$158.13$ μs], [$0.465$],
    [Rb $|7 P_(3/2) angle.r$],  [$0.131$ μs],[$0.2693$ μs], [$2.056$],
  ),
)

The paper's Rydberg lifetimes are about $2 times$ *longer* than
ARC gives at $T = 300$ K. The intermediate $|P angle.r$ state goes
the *opposite* way -- the paper quotes a $|P angle.r$ lifetime
about $2 times$ shorter than ARC. We do not resolve this discrepancy
here; we only document it.

Decision: *use ARC's values for the two Rydberg lifetimes* (since
both states are being swapped out anyway, and ARC is internally
consistent across `main.tex`'s pair and Option A); *keep the paper's
$tau_P = 0.131$ μs* since the $|7 P_(3/2) angle.r$ intermediate
state is unchanged by the level swap. Changing $tau_P$ would
additionally alter the OR gate's Raman-scattering baseline and muddy
the comparison.

This choice is conservative: for a fair comparison against the
paper's own numbers, one would have to *also* re-run the paper pair
with ARC lifetimes, which would likewise depress the paper-pair
fidelity. We note this but do not execute the counterfactual.

= Coefficients at $a = 5.00$ μm

#figure(
  caption: [Rydberg interaction strengths and dimensionless ratios for
    the paper pair (left) and Option A (right), both at $a = 5.00$ μm.
    Ratios marked in bold are the warning signs.],
  kind: table,
  table(
    columns: (auto, auto, auto, auto),
    align: (left, center, center, center),
    stroke: 0.4pt,
    table.header(
      [*Quantity*], [*paper pair*], [*Option A*], [*ratio*],
    ),
    [Rydberg pair],
      [Rb $69 D_(5/2)$ + Cs $79 D_(5/2)$],
      [Rb $66 D_(5/2)$ + Cs $76 D_(3/2)$],
      [],
    [Förster channel],
      [$70 P_(3/2)$--$77 F_(7/2)$],
      [$67 P_(3/2)$--$74 F_(5/2)$],
      [],
    [$tilde(C)_3$ (GHz$dot$μm³)], [$26.3$],    [$22.84$],   [$times 0.868$],
    [$C_6^"CsCs"$ (GHz$dot$μm⁶)], [$-1449$],   [$-692.9$],  [$times 0.478$],
    [$V_("ct") \/ (2 pi)$ (MHz)], [$+595.10$], [$+516.81$], [$times 0.868$],
    [$V_(c c) \/ (2 pi)$ (MHz)],  [$-11.59$],  [$-5.54$],   [$times 0.478$],
    [selectivity $V_("ct") \/ |V_(c c)|$],
                                  [$51.3$],    [$93.2$],    [$times 1.82$],
    [$V_("ct") \/ Delta$ (OR)],   [$1.190$],   [$bold(1.034)$], [$times 0.868$],
    [$V_("ct") \/ Omega_p$ (OR)], [$11.90$],   [$10.34$],   [$times 0.868$],
    [$V_("ct") \/ Omega_t$ (CCX)],[$11.90$],   [$10.34$],   [$times 0.868$],
    [$|V_(c c)| \/ Omega_c$],     [$0.232$],   [$0.111$],   [$times 0.478$],
    [Blockade fidelity $P_(1 r)$],[$0.9998$],  [$0.9997$],  [$approx$ equal],
  ),
)

The headline physics is what we expected: the selectivity ratio is
$1.8 times$ better, and $|V_(c c)| \/ Omega_c$ is *halved*. But the
$V_("ct") \/ Delta$ ratio drops to $1.034$ -- within $3%$ of the
EIT-breaking wall. This is the same concern we had at $a = 5.25$ μm
with the original pair, and it means the OR gate's single-blockade
branches will lose $approx 6 times 10^(-4)$ fidelity per input even
at $a = 5.00$ μm now that we have switched pairs.

*Source of the $C_3$ / $C_6$ values:*
B. J. Ireland, J. D. Pritchard, J. P. Shaffer,
*"Interspecies Förster resonances of Rb–Cs Rydberg d-states for
enhanced multi-qubit gate fidelities"*,
Phys. Rev. Research *6*, 013293 (2024), arXiv:2401.02308, Table I,
row 1 ($theta = 90$ deg, quantization axis perpendicular to array
plane). These values are tabulated locally in
`reference/energy_level_inter.md`.

= Pulse parameter change: $Omega_R = 2.9 thick Omega_p$

The OR gate's target-side dynamics live on three levels,
$|A angle.r$ / $|B angle.r$ / $|P angle.r$ / $|R angle.r$, with
$Omega_p$ coupling $|A/B angle.r <-> |P angle.r$ (two-photon Raman,
detuned by $Delta$) and $Omega_R$ coupling $|P angle.r <-> |R angle.r$
(resonant). The EIT dark state in the absence of blockade is
approximately

$ |"dark" angle.r prop Omega_R |A angle.r - Omega_p |R angle.r, $

so the ratio $Omega_R \/ Omega_p$ sets the *mixing angle* between the
computational and the Rydberg parts of the dark state. A larger ratio
gives a cleaner dark state (target preserved more faithfully when the
blockade is absent).

Empirically, the paper uses $Omega_R = 3.5 thick Omega_p$. This run
uses $Omega_R = 2.9 thick Omega_p$, which reduces the ratio by about
$17%$. We observe that the no-blockade $|0,0,* angle.r$ OR branch
degrades substantially, primarily because the residual $|R angle.r$
population grows: the diagnostic shows
$P(|R angle.r) = 0.0030$ at $Omega_R = 2.9 thick Omega_p$ with the
Option 1 pulse, compared to $0.0002$ at the paper's $3.5$ ratio --
a $15 times$ jump. (With the previous fudged pulse at the same
Omega_R ratio, the residual was $0.0038$, so the Option 1 cleanup
itself cuts it by $approx 21%$.) This is the dominant contributor
to the Option A infidelity in this run.

The sign of this effect is predictable; the magnitude is not an
Option A feature but a consequence of the pulse-ratio choice. If
$Omega_R$ is reverted to $3.5 thick Omega_p$, the $|0,0,* angle.r$
OR branch should recover to $approx 0.999$ and the Option A OR
average fidelity should rise to approximately the $0.9956$ level
we extrapolated last turn.

= OR gate results

#figure(
  caption: [OR gate: per-input state fidelities in Option A
    compared to the three earlier runs on this pulse template.
    Columns (a)–(c) are from the previous reports; column (d) is
    the Option A run with ARC lifetimes and
    $Omega_R = 2.9 thick Omega_p$.],
  kind: table,
  table(
    columns: (auto, auto, auto, auto, auto, auto),
    align: (center, left, center, center, center, center),
    stroke: 0.4pt,
    table.header(
      [*\#*], [*Input*],
      [*(a) no $V_(c c)$, $a=5.00$*],
      [*(b) paper, $a=5.00$*],
      [*(c) $V_(c c)$ override, $a=5.25$*],
      [*(d) Option A, $a=5.00$*],
    ),
    [1], [$|0,0,A angle.r$], [$0.999186$], [$0.999186$], [$0.999186$], [$bold(0.993561)$],
    [2], [$|0,0,B angle.r$], [$0.999186$], [$0.999186$], [$0.999186$], [$bold(0.993561)$],
    [3], [$|0,1,A angle.r$], [$0.995761$], [$0.995761$], [$0.995184$], [$0.995259$],
    [4], [$|0,1,B angle.r$], [$0.995761$], [$0.995761$], [$0.995184$], [$0.995259$],
    [5], [$|1,0,A angle.r$], [$0.995761$], [$0.995761$], [$0.995184$], [$0.995259$],
    [6], [$|1,0,B angle.r$], [$0.995761$], [$0.995761$], [$0.995184$], [$0.995259$],
    [7], [$|1,1,A angle.r$], [$0.996042$], [$0.979678$], [$0.992780$], [$bold(0.983351)$],
    [8], [$|1,1,B angle.r$], [$0.996042$], [$0.979678$], [$0.992780$], [$bold(0.983351)$],
  ),
)

#block(
  fill: rgb("#f4faff"),
  inset: 10pt,
  radius: 4pt,
  width: 100%,
)[
  *OR gate average fidelity (Option A + Option 1 pulse):*
  $ quad overline(F)_"OR" = bold(0.991858), quad 1 - overline(F)_"OR" = bold(8.14 times 10^(-3)). $
  *Changes:*
  $ quad Delta overline(F) = -7.4 times 10^(-4)$ vs. the paper pair at $a = 5.00$ μm,
  $ quad Delta overline(F) = -3.73 times 10^(-3)$ vs. column (c),
  $ quad Delta overline(F) = +2.6 times 10^(-4)$ vs. Option A with the fudged pulse.
]

Per-input breakdown of the two rows that moved significantly
going (c) $arrow$ (d):

- The $|0,0,* angle.r$ branch drops from $0.999186$ to $0.993561$
  -- a loss of $5.6 times 10^(-3)$ per input. Target-level populations
  at $|0,0,A angle.r$ now read
  $(P_A, P_B, P_P, P_R) = (0.9936, 0.0034, 0.0001, 0.0030)$; the
  $0.0030$ residual $|R angle.r$ population is *15$times$ larger* than
  at $Omega_R = 3.5 thick Omega_p$, and it is the main source of the
  $|0,0,* angle.r$ infidelity. This is a pulse-ratio effect, not a
  level-swap effect. (The Option 1 pulse cleanup shaved this residual
  from $0.0038$ to $0.0030$ relative to the previous Option A run
  with the fudged pulse -- a small secondary win from smoother edges.)

- The $|1,1,* angle.r$ branch drops from $0.992780$ (the override run)
  to $0.983351$. The cause is twofold: (i) ARC lifetimes are
  roughly half the paper anchors, adding $approx 3 times 10^(-3)$
  decay infidelity on this branch; (ii) the accumulated $V_(c c)$
  phase has slid from $+1.00 pi$ (override) to $approx +0.67 pi$
  (Option A derived value), landing at a different point on the
  phase-interference curve we characterised two reports ago.

The single-blockade branches $|0,1,* angle.r$ and $|1,0,* angle.r$
are essentially unchanged between (c) and (d) because
$V_("ct") \/ Delta$ is the same ($1.028$ vs. $1.034$) and the
ARC-lifetime decay contribution on the $320$ ns gate is roughly
$+4 times 10^(-4)$, at the edge of what this metric can resolve.
The Option 1 pulse cleanup lifts these branches from $0.995147$ to
$0.995259$, a $+1.1 times 10^(-4)$ per-input improvement consistent
with slightly better EIT dark-state tracking during the pulse edges.

= CCX gate results

#figure(
  caption: [CCX gate: per-input state fidelities across the same four
    conditions.],
  kind: table,
  table(
    columns: (auto, auto, auto, auto, auto, auto),
    align: (center, left, center, center, center, center),
    stroke: 0.4pt,
    table.header(
      [*\#*], [*Input*],
      [*(a) no $V_(c c)$*],
      [*(b) paper*],
      [*(c) override, $a=5.25$*],
      [*(d) Option A, $a=5.00$*],
    ),
    [1], [$|0,0,A angle.r$], [$0.999512$], [$0.990583$], [$0.993746$], [$bold(0.992175)$],
    [2], [$|0,0,B angle.r$], [$0.999282$], [$0.990044$], [$0.994944$], [$bold(0.991581)$],
    [3], [$|0,1,A angle.r$], [$0.999736$], [$0.999736$], [$0.997594$], [$0.996745$],
    [4], [$|0,1,B angle.r$], [$0.999307$], [$0.999307$], [$0.997183$], [$0.997920$],
    [5], [$|1,0,A angle.r$], [$0.999736$], [$0.999736$], [$0.997594$], [$0.996745$],
    [6], [$|1,0,B angle.r$], [$0.999307$], [$0.999307$], [$0.997183$], [$0.997920$],
    [7], [$|1,1,A angle.r$], [$0.999968$], [$0.999968$], [$0.999968$], [$0.999937$],
    [8], [$|1,1,B angle.r$], [$0.999968$], [$0.999968$], [$0.999968$], [$0.999937$],
  ),
)

#block(
  fill: rgb("#f4faff"),
  inset: 10pt,
  radius: 4pt,
  width: 100%,
)[
  *CCX gate average fidelity (Option A):*
  $ quad overline(F)_"CCX" = bold(0.996620), quad 1 - overline(F)_"CCX" = bold(3.38 times 10^(-3)). $
  *Changes:*
  $ quad Delta overline(F) = -0.71 times 10^(-3)$ vs. the paper pair at $a = 5.00$ μm,
  $ quad Delta overline(F) = -6.5 times 10^(-4)$ vs. column (c).
]

For CCX the story is simpler. The $Omega_R$ change has no effect here
(CCX does not use $Omega_R$). The main contributors to the Option A
drop relative to the paper are:

+ *ARC lifetimes* add decay infidelity across all branches. This
  contributes roughly $3$--$5 times 10^(-4)$ uniformly, which shows
  up as the small $0.001$--$0.002$ drops across rows 1–6.

+ *$V_("ct") \/ Delta = 1.034$* instead of the paper's $1.19$: the
  single-blockade branches $|0,1,* angle.r$ / $|1,0,* angle.r$ lose
  about $2 times 10^(-3)$ per input, as they did in the $a = 5.25$
  μm run. This is the same $V_("ct")$-margin issue the previous
  report flagged.

+ *Weaker $V_(c c)$* ($-5.54$ vs.\ $-11.59$ MHz) *helps* the
  $|0,0,* angle.r$ branch. At Option A, $|0,0,A angle.r$ sits at
  $0.9922$ instead of $0.9906$ -- a recovery of $approx 1.6 times 10^(-3)$.
  But not enough to overcome the other two effects.

The fact that the CCX numbers are not catastrophically different
from either the paper or the override run suggests Option A is in
the right ballpark; a cleaner run (without the $Omega_R$ change) would
show the full level-swap gain.

= Five-way summary comparison

#figure(
  caption: [All parameter conditions we have now simulated on this
    pulse template. Column (d) is Option A proper; column (d$'$) is
    what we expected from the two-term model extrapolation last turn
    (using paper-lifetime anchors). The gap between (d) and (d$'$) is
    the combined effect of ARC lifetimes and the $Omega_R$ change.],
  kind: table,
  table(
    columns: (auto, auto, auto, auto, auto, auto),
    align: (left, center, center, center, center, center),
    stroke: 0.4pt,
    table.header(
      [*Metric*],
      [*(a) no $V_(c c)$*],
      [*(b) paper*],
      [*(c) override*],
      [*(d) Option A*],
      [*(d$'$) expected*],
    ),
    [$a$ (μm)],                    [$5.00$],   [$5.00$],   [$5.25$],   [$5.00$],   [$5.00$],
    [Rb level],                    [$69 D_(5/2)$], [$69 D_(5/2)$], [$69 D_(5/2)$], [$66 D_(5/2)$], [$66 D_(5/2)$],
    [Cs level],                    [$79 D_(5/2)$], [$79 D_(5/2)$], [$79 D_(5/2)$], [$76 D_(3/2)$], [$76 D_(3/2)$],
    [$Omega_R \/ Omega_p$],        [$3.5$],    [$3.5$],    [$3.5$],    [$bold(2.9)$], [$3.5$],
    [$tau_r$ (μs)],                [$340$],    [$340$],    [$340$],    [$bold(142.73)$], [$340$ (paper)],
    [$tau_R$ (μs)],                [$260$],    [$260$],    [$260$],    [$bold(134.87)$], [$260$ (paper)],
    [$V_("ct") \/ (2 pi)$ (MHz)],  [$+595$],   [$+595$],   [$+514$],   [$+517$],   [$+517$],
    [$V_(c c) \/ (2 pi)$ (MHz)],   [$0$],      [$-11.59$], [$+5.00$],  [$-5.54$],  [$-5.54$],
    [$V_("ct") \/ Delta$],         [$1.19$],   [$1.19$],   [$1.028$],  [$1.034$],  [$1.034$],
    [],                            [],         [],         [],         [],         [],
    [OR $overline(F)$],            [$0.996688$], [$0.992597$], [$0.995584$], [$bold(0.991858)$], [$approx 0.9956$],
    [OR $1 - overline(F)$],        [$3.31 times 10^(-3)$], [$7.40 times 10^(-3)$], [$4.42 times 10^(-3)$], [$bold(8.14 times 10^(-3))$], [$approx 4.4 times 10^(-3)$],
    [OR $|0,0,* angle.r$ $F_k$],   [$0.999186$], [$0.999186$], [$0.999186$], [$bold(0.993561)$], [$approx 0.999$],
    [OR $|1,1,* angle.r$ $F_k$],   [$0.996042$], [$0.979678$], [$0.992780$], [$0.983351$], [$approx 0.993$],
    [OR single-blk. $F_k$],        [$0.995761$], [$0.995761$], [$0.995184$], [$0.995259$], [$0.995$],
    [],                            [],         [],         [],         [],         [],
    [CCX $overline(F)$],           [$0.999602$], [$0.997331$], [$0.997273$], [$bold(0.996620)$], [$approx 0.9973$],
    [CCX $1 - overline(F)$],       [$3.98 times 10^(-4)$], [$2.67 times 10^(-3)$], [$2.73 times 10^(-3)$], [$bold(3.38 times 10^(-3))$], [$approx 2.7 times 10^(-3)$],
  ),
)

The column (d$'$) expectations are *not* from a simulation; they are
the two-term-model extrapolation from last turn, which implicitly
used the paper's (optimistic) lifetimes and $Omega_R = 3.5 thick Omega_p$.
The actual Option A run (d) is approximately $3.7 times 10^(-3)$
worse on OR and $7 times 10^(-4)$ worse on CCX than (d$'$). That gap
is entirely attributable to effects (b) and (c) in §1 -- ARC lifetimes
and the $Omega_R$ ratio change -- not to the level swap itself.

= Error decomposition

A rough budget, all contributions per-input infidelity on the OR gate,
fitted on the differences between the five columns above:

#figure(
  caption: [Approximate contributions to the OR gate infidelity at
    Option A, extracted from the differences between simulated
    conditions. The absolute baseline is the $V_(c c) = 0$ paper-pair
    infidelity of $3.31 times 10^(-3)$.],
  kind: table,
  table(
    columns: (auto, auto, auto),
    align: (left, center, center),
    stroke: 0.4pt,
    table.header(
      [*Contribution*], [*estimate*], [*mechanism*],
    ),
    [Baseline (paper pair, no $V_(c c)$)],
      [$3.3 times 10^(-3)$],
      [finite $V_("ct") \/ Delta$, Raman scattering],
    [$V_(c c) approx -5.5$ MHz on $|1,1,* angle.r$],
      [$approx +7.2 times 10^(-4)$],
      [residual phase on doubly-excited branch],
    [$V_("ct") \/ Delta : 1.19 arrow 1.034$ on single-blockade],
      [$approx +3 times 10^(-4)$],
      [$(Omega_p \/ V_("ct"))^2$ growth],
    [*ARC lifetimes* ($tau$ halved) on all branches],
      [$approx +1.5 times 10^(-3)$],
      [$T_1$ decay during $320$ ns gate],
    [*$Omega_R = 2.9$* on $|0,0,* angle.r$ branch],
      [$approx +2.8 times 10^(-3)$],
      [imperfect EIT dark state, $|R angle.r$ leakage],
    [Option 1 pulse cleanup],
      [$approx -2.6 times 10^(-4)$],
      [smoother edges, better dark-state adiabaticity],
    [*Total predicted* $1 - overline(F)_"OR"$],
      [$approx 8.3 times 10^(-3)$],
      [],
    [*Actual simulated*],
      [$bold(8.14 times 10^(-3))$],
      [matches to $1.6 times 10^(-4)$],
  ),
)

The arithmetic works: the ARC lifetime update and the $Omega_R$ change
together account for $4.3 times 10^(-3)$ of infidelity, which is
*larger than the entire baseline* of the paper-pair OR gate.
These two items are the difference between the Option A "target"
number from the last turn's extrapolation ($0.9956$) and the actual
simulation ($0.9919$). The Option 1 pulse cleanup claws back a
modest $2.6 times 10^(-4)$ but is not enough to close the gap.

Neither of them is a property of Option A itself -- they are external
corrections that *would have applied to the paper pair as well*, if
we had re-run the paper pair with ARC lifetimes and
$Omega_R = 2.9 thick Omega_p$. A fair head-to-head comparison between
the paper pair and Option A would need to use the same lifetime
source and the same $Omega_R$ ratio for both.

= What would "clean Option A" look like?

Two independent counterfactual runs would isolate each effect. They
are not executed here but would be a straightforward one-line change
to the scripts:

+ *Option A with paper lifetimes and* $bold(Omega_R = 3.5 thick Omega_p)$.
  Revert both $gamma_r$/$gamma_R$ to $1/340$/$1/260$ and
  `omega_R_amp` to $3.5 thick Omega_p$. Expected outcome, from the
  two-term model in the previous report:
  $overline(F)_"OR" approx 0.9956$, $overline(F)_"CCX" approx 0.9973$.
  This is the *pure level-swap* benefit of Option A and is a fair
  apples-to-apples comparison against the paper's own quoted numbers.

+ *Paper pair with ARC lifetimes and* $bold(Omega_R = 2.9 thick Omega_p)$.
  Re-run `examples/Average_fidelity/` with the new $gamma$'s and the
  new `omega_R_amp`. Expected outcome: the paper pair's $overline(F)_"OR"$
  would also drop substantially -- probably to around $0.988$–$0.990$
  -- because the ARC-lifetime and $Omega_R$ penalties apply there
  just as much. In that comparison Option A would almost certainly
  *win*, because its lower $V_(c c)$ partially offsets the shared
  pulse-ratio cost.

Either counterfactual would tell you what you actually want to know:
"does swapping to Option A help, when everything else is held fixed?"
The answer is almost certainly yes -- the level swap alone is worth
$approx +3 times 10^(-3)$ on OR gate fidelity -- but this run's
result is a *joint* measurement that bundles the level swap with two
unfavourable external changes.

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

Source modifications in this run:

+ `triqg/hamiltonian.py`: unchanged.

+ `examples/Average_fidelity_Vcc_newenergy/or_average_gate_fid_gaussian.py`:
  - `omega_R_amp = 3.5 * omega_p_amp` → `omega_R_amp = 2.9 * omega_p_amp`
  - `omega_p_amp = 2*np.pi*50.0 * 1.039975` → `omega_p_amp = 2*np.pi*50.0` (Option 1)
  - `sigma = 0.0014` → `sigma = 0.001771` (Option 1, preserves pulse area $pi/4$)
  - `a_um = 5.25` → `a_um = 5.00`
  - `C3_tilde = 26.3` → `C3_tilde = 22.84` (Ireland et al. 2024)
  - `C6_CsCs = -1449.0` → `C6_CsCs = -692.9` (Ireland et al. 2024)
  - `V_cc` override block removed; `V_cc` now derived from $C_6$.
  - `gamma_r = 1/340` → `gamma_r = 1/142.73` (ARC, $T=300$ K)
  - `gamma_R = 1/260` → `gamma_R = 1/134.87` (ARC, $T=300$ K)
  - `gamma_P = 1/0.131` unchanged (paper value retained)
  - Docstring and print statements updated with Option A header and
    lifetime / source citations.
  - Plot of the old vs. new $Omega_p(t)$ pulse saved to
    `omega_p_pulse_comparison.png` in the same folder.

+ `examples/Average_fidelity_Vcc_newenergy/ccx_average_gate_fidelity.py`:
  - Same updates except the $Omega_R$ change, which the CCX gate
    does not use.

All other pulse parameters, timing windows, detunings, and solver
options remain at the paper-consistent values documented in
`results/Final_three_qubit_gate_report.typ`.

= Conclusion

+ *The Option A level swap works as advertised in principle*, but
  this run bundles it with two external changes (ARC lifetimes and
  $Omega_R = 2.9 thick Omega_p$) that together cost more infidelity
  than the level swap saves. The net result is that the Option A
  average fidelities come out *lower* than the paper baseline:
  $ overline(F)_"OR" = 0.991858, quad overline(F)_"CCX" = 0.996620. $

+ *The dominant cost is the $Omega_R$ ratio change* ($3.5 arrow 2.9$),
  which degrades the no-blockade OR $|0,0,* angle.r$ branch by
  $5.6 times 10^(-3)$ per input through an imperfect EIT dark
  state. This is the single largest term in the OR infidelity
  budget and is not a property of Option A.

+ *The Option 1 pulse cleanup gives a small positive effect*:
  removing the $1.039975$ fudge factor and rescaling
  $sigma: 0.0014 arrow 0.001771$ preserves the pulse area at $pi/4$
  exactly, but the slightly smoother edges improve EIT dark-state
  adiabaticity on $|0,0,* angle.r$, lifting
  $overline(F)_"OR"$ by $+2.6 times 10^(-4)$ and cutting the
  residual $|R angle.r$ population by $approx 21%$. This is a
  genuine small win, not just a cosmetic reparametrisation.

+ *The second cost is the ARC lifetime update*, which is a factor
  $approx 2$ more pessimistic than the paper's quoted values and adds
  roughly $10^(-3)$ infidelity per branch on both gates. This is
  a fair update -- the paper's quoted numbers are hard to
  reconcile with standard ARC / Beterov calculations -- but it
  would affect the paper pair equally.

+ *For a clean test of Option A*, run the counterfactuals in §8:
  either revert $Omega_R$ to $3.5 thick Omega_p$ alone, or also
  revert lifetimes to the paper anchors. Both isolate the
  level-swap effect.

+ *Practical operational recommendation:* if the $Omega_R = 2.9$
  choice is essential for some other reason, the OR gate will be
  dominated by the $|0,0,* angle.r$ residual. If $Omega_R$ is a
  free parameter, $3.5$ (or higher) gives a markedly better OR
  gate regardless of whether the pair is the paper's or Option A.

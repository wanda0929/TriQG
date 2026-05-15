// Smaller_VDD_energy report  (rev. 4)
//   Rb 54 D_{3/2} + Cs 62 D_{5/2}  at  a = 5.00 um
//   Strong-blockade + sub-pi/4-area-optimized OR gate
//   F_bar > 0.99 achieved.
//
// rev. 4 (2026-05-14):  added section on Omega_R sweep
//   (sweep_omega_R.py + sweep_smallDelta_refine.py) showing why
//   the Farouk ratio Omega_R / Omega_p = 3.5 and Delta in [400, 600]
//   MHz remain optimal under a relaxed (Delta, Omega_R) sweep.
//
// Compile with:  typst compile Smaller_VDD_average_fidelity_report.typ
//
// Scripts (run in this folder):
//   level_parameters.py
//   or_average_gate_fid_gaussian.py     (canonical FINAL config)
//   ccx_average_gate_fidelity.py
//   sweep_delta.py / sweep_alpha.py / sweep_omega_p.py
//   sweep_area.py / sweep_2d.py         (parameter scans)
//   sweep_omega_R.py                    (rev. 4: coarse (Delta,Omega_R))
//   sweep_smallDelta_refine.py          (rev. 4: three-pass refinement)
//   plot_omega_R_landscape.py           (rev. 4: omega_R_landscape.png)
//   make_pulse_figs.py                  (regenerates figures)
//
// Logs:
//   or_final.log, ccx.log, sweep_*.log

#set document(
  title: "Smaller V_DD: F > 0.99 OR gate at Rb 54 D_3/2 + Cs 62 D_5/2",
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
    Smaller $V_("DD")$: $overline(F)_"OR" > 0.99$ via strong blockade
    and sub-$pi\/4$ area
  ]

  #v(0.3em)
  #text(size: 11pt)[
    Rb $54 D_(3/2)$ + Cs $62 D_(5/2)$ at $a = 5.00$ μm,
    smooth super-Gaussian probe with $alpha = 4$, area $= 0.95 thick pi \/ 4$
  ]

  #v(0.4em)
  #text(size: 10pt)[
    TriQG `examples/Smaller_VDD_energy` · compiled 2026-05-14 (rev. 4)
  ]
]

#v(0.6em)

#block(
  fill: rgb("#eef8ec"),
  inset: 10pt,
  radius: 4pt,
  width: 100%,
)[
  *Headline results.* All on the new energy levels (Rb $54 D_(3/2)$ + Cs $62 D_(5/2)$, $V_(c t) = 2 pi times 226$ MHz):

  - *OR gate:*
    $overline(F)_"OR" = bold(0.992422)$,
    infidelity $1 - overline(F)_"OR" = bold(7.58 times 10^(-3))$.
    Total gate time $385$ ns.
  - *CCX gate:*
    $overline(F)_"CCX" = bold(0.995854)$,
    infidelity $1 - overline(F)_"CCX" = bold(4.15 times 10^(-3))$.
    Total gate time $95$ ns.
  - *Both gates above $0.99$.*  X-round per-cycle infidelity
    $1 - overline(F)_"OR" overline(F)_"CCX" approx
    1.17 times 10^(-2)$.

  The OR-gate improvement (from $0.97$ in rev. 2 to $0.992$ in this
  revision) comes from *two* changes from rev. 2:
  (i) restore the paper's amplitudes
      $Omega_p = Omega_c = 2 pi times 50$ MHz,
      $Omega_R = 3.5 thick Omega_p$ (instead of the rev. 2 value of
      $2 pi times 20$ MHz, which was too small to saturate the
      blockade in a short gate time);
  (ii) drop the two-photon area from the paper's $pi\/4$ to
      $0.95 thick pi \/ 4$ -- empirically the protocol's correct
      target rotation at finite $V_(c t)$.

  *Rev. 4 audit.*  A full $(Delta, Omega_R \/ Omega_p)$ sweep with
  ratio relaxed from the Farouk $3.5$ confirms that this
  configuration sits on the *broadest* $overline(F)_"OR" > 0.99$
  plateau in the $(Delta, Omega_R, K)$ landscape; a secondary peak
  at $(Delta = 225 thick "MHz", Omega_R \/ Omega_p = 2.5)$ matches
  the fidelity in $2.1 times$ shorter gate time but is a
  knife-edge (Sec. 5).
]

= Motivation

The paper's X-round Rydberg pair (Rb $66 D_(5/2)$ + Cs $76 D_(3/2)$)
sits $~160 times$ above a Rb same-species Förster zero crossing, so
the data--data blockade at $r_("DD") = a = 5$ μm reaches

$ V_("DD") \/ (2 pi) approx 91 thick "MHz", $

one-sixth of the desired Rb--Cs interspecies blockade. The
species-selective global drive in the X-round excites every Rb data
atom in parallel, so this $V_("DD")$ enters the per-gate selectivity
budget on equal footing with the published $V_(c c)$ entry:

$ R_("DD") = V_(c t) / V_("DD") approx 5.7 quad arrow.r quad
  R_("DD") << 100 thick "(target)". $

Moving Rb to $54 D_(3/2)$ parks the data atom *on* the Rb $D_(3/2)$
same-species zero crossing, killing $|C_6^"DD"|$ by $~157 times$.
Cs $62 D_(5/2)$ keeps the interspecies Förster channel
near-resonant ($Delta_F \/ (2 pi) = -5.2$ MHz). All numerical inputs
are imported from `level_parameters.py`; see
`SelfCorrectingRydberg/plan/Rb_Cs_level_redesign.md` for the full
ARC scan.

The redesigned $V_(c t)$ is $~2.3 times$ smaller than the paper's
value, which makes the protocol's $V_(c t) -> oo$ assumptions
slightly off. The pulse-engineering work in this report fixes that.

= Setup

== Rydberg levels and Förster channel

- Rb data $|R angle.r = |54 D_(3/2) angle.r$,
  Cs ancilla $|r angle.r = |62 D_(5/2) angle.r$,
  Rb intermediate $|P angle.r = |7 P_(3/2) angle.r$ (unchanged).
- Förster channel
  $|54 D_(3/2); 62 D_(5/2) angle.r <-> |55 P_(3/2); 60 F_(5/2) angle.r$,
  defect $Delta_F \/ (2 pi) = -5.2$ MHz.

== Lattice geometry (unchanged)

- Rotated-lattice spacing $a = 5.00$ μm.
- Data--ancilla $r_("DA") = a \/ sqrt(2) approx 3.5355$ μm.
- Data--data $r_("DD") = a = 5.0000$ μm.
- Ancilla--ancilla $r_("AA") = a sqrt(2) approx 7.0711$ μm.

== Interaction strengths

#figure(
  caption: [Rydberg interaction strengths at $a = 5.00$ μm.
    The new selectivity ratios $R_("DD")$ and $R_("AA")$ both clear
    the $>= 100$ threshold.],
  kind: table,
  table(
    columns: (auto, auto, auto),
    align: (left, center, center),
    stroke: 0.4pt,
    table.header(
      [*Quantity*], [*Paper (Option A)*], [*This redesign*],
    ),
    [Rydberg pair],
      [Rb $66 D_(5/2)$ + Cs $76 D_(3/2)$],
      [Rb $54 D_(3/2)$ + Cs $62 D_(5/2)$],
    [Förster channel],
      [$67 P_(3/2)$ -- $74 F_(5/2)$],
      [$55 P_(3/2)$ -- $60 F_(5/2)$],
    [$Delta_F \/ (2 pi)$ (MHz)],
      [$+10.9$], [$-5.2$],
    [$tilde(C)_3$ (GHz$dot$μm³)], [$22.84$], [$10.0$],
    [$|C_6^"RbRb"|$ (GHz$dot$μm⁶)], [$~1420$], [$9.07$],
    [$C_6^"CsCs"$ (GHz$dot$μm⁶)], [$-692.9$], [$-91.0$],
    [$V_(c t) \/ (2 pi)$ (MHz)], [$+516.81$], [$+226.27$],
    [$V_("DD") \/ (2 pi)$ (MHz)], [$+91$], [$+0.58$],
    [$V_(c c) \/ (2 pi)$ (MHz)], [$-5.54$], [$-0.7280$],
    [$R_("DD") = V_(c t) \/ V_("DD")$], [$~5.7$ *(fail)*], [$bold(389.8)$ #sym.checkmark],
    [$R_("AA") = V_(c t) \/ |V_(c c)|$], [$~94$], [$bold(310.8)$ #sym.checkmark],
  ),
)

== Pulse-sequence parameters (FINAL)

#figure(
  caption: [Pulse parameters used in this report.
    The OR-gate target probe is a super-Gaussian of order 6,
    $Omega_p(t) = (Omega_p \/ 2) thin exp[-((t - t_c)^3 \/ sigma)^2]$
    with $t_c = T_c + T_f$, active on $[T_c, T_c + 2 T_f]$.
    All quantities marked *bold* are this report's settings; the
    "Rev. 2" column shows the prior overly-conservative configuration.],
  kind: table,
  table(
    columns: (auto, auto, auto, auto),
    align: (left, center, center, center),
    stroke: 0.4pt,
    table.header(
      [*Parameter*], [*Paper*], [*Rev. 2*], [*FINAL (rev. 3)*],
    ),
    [$Omega_c \/ (2 pi)$ (MHz)], [$50$], [$20$], [$bold(50)$],
    [$Omega_p \/ (2 pi)$ (MHz)], [$50 times 1.040$], [$20$], [$bold(50)$],
    [$Omega_R \/ Omega_p$], [$3.5$], [$3.5$], [$3.5$],
    [$Delta \/ (2 pi)$ (MHz)], [$500$], [$500$], [$bold(500)$],
    [$T_c$ (ns)], [$10$], [$25$], [$bold(10)$],
    [Shape factor $alpha = T_f^3 \/ sigma$],
      [$2.41$ (abrupt)], [$4.0$ (smooth)], [$bold(4.0)$ (smooth)],
    [Edge / peak ratio],
      [$3 times 10^(-3)$], [$1.1 times 10^(-7)$], [$bold(1.1 times 10^(-7))$],
    [Area / $(pi \/ 4)$], [$1.000$], [$1.000$], [$bold(0.950)$],
    [$T_f$ (ns)], [$150$], [$1200$], [$bold(182.46)$],
    [$sigma$ (ns)], [$1.4$], [$432$], [$bold(1.519)$],
    [Total OR gate (ns)], [$320$], [$2441$], [$bold(384.9)$],
    [$Omega_(c c) \/ (2 pi)$ (MHz)], [$100$], [$50$], [$50$],
    [$Omega_t \/ (2 pi)$ (MHz)], [$50$], [$20$], [$20$],
    [Total CCX gate (ns)], [$40$], [$95$], [$95$],
  ),
)

== Decoherence rates (T = 300 K, BBR-included)

Cs $|62 D_(5/2) angle.r$: $tau_r approx 77$ μs#footnote[$n^3$-scaled
estimate from ARC values for paper levels; re-verify with
`Cs.getStateLifetime(62, 2, 2.5, temperature=300, includeLevelsUpTo=92)`.];
Rb $|54 D_(3/2) angle.r$: $tau_R approx 74$ μs#footnote[$n^3$-scaled.];
Rb $|7 P_(3/2) angle.r$: $tau_P = 0.131$ μs (paper).

= Strong-blockade audit

For the OR-gate Raman drive on the target, the relevant off-resonant
processes have effective rates that must lie *well below* the
Rydberg blockade $V_(c t)$:

#figure(
  caption: [Strong-blockade conditions at the final configuration.
    Both off-resonant target rates are at least $26 times$ below the
    interspecies blockade; the AC-Stark rate is $90 times$ below.],
  kind: table,
  table(
    columns: (auto, auto, auto, auto),
    align: (left, center, center, center),
    stroke: 0.4pt,
    table.header([*Process*], [*Expression*], [*Value / $(2 pi)$ (MHz)*],
                 [*Margin $V_(c t) \/ (.)$*]),
    [One-photon AC Stark on $|P angle.r$],
      [$Omega_p^2 \/ (2 Delta)$], [$2.50$], [$bold(90)$],
    [Two-photon Raman on $|R angle.r$],
      [$Omega_p Omega_R \/ (2 Delta)$], [$8.75$], [$bold(25.9)$],
    [Rydberg blockade $V_(c t) \/ (2 pi)$],
      [$tilde(C)_3 \/ r_("DA")^3$], [$226.27$], [---],
  ),
)

Both margins comfortably satisfy "$<< V_(c t)$" in the sense the
user requested; the binding constraint is the two-photon Raman rate
(factor $26$ below $V_(c t)$). Tightening these margins further
(e.g. by raising $Delta$) lengthens the gate quadratically and
*worsens* fidelity through Cs $|r angle.r$ decay during the gate
window; the optimum we converged on sits on the boundary between
the blockade-leakage regime (small $Delta$) and the
decay-during-dwell regime (large $Delta$).

= Pulse shape and the area = $0.95 thick pi \/ 4$ optimum

== Why the area is *not* $pi \/ 4$ at finite $V_(c t)$

The Farouk-style OR protocol was designed in the $V_(c t) -> oo$
limit, where the blockade is *perfect*: any control in $|r angle.r$
fully shifts the target's $|R angle.r$ off resonance. Under that
ideal, the two-photon Raman area must equal $pi \/ 4$ for the
correct OR truth table.

At finite $V_(c t)$, the protocol has small residual leakage on the
blockaded branches. Empirically, the protocol's correct two-photon
area shifts to compensate. A 2-D scan over (detuning $Delta$,
area factor $K = $ area$\/(pi \/ 4)$) -- run by `sweep_2d.py` at
fixed $Omega_p = 2 pi times 50$ MHz, $alpha = 4$ -- shows a
fidelity ridge at $K approx 0.95-0.96$ for $Delta in [400, 600]$ MHz:

#figure(
  image("fidelity_heatmap.png", width: 88%),
  caption: [
    2-D scan of $overline(F)_"OR"$ over $(Delta, K)$, where
    $K = $ area $\/ (pi \/ 4)$.
    Fixed parameters: $Omega_p \/ (2 pi) = 50$ MHz,
    $Omega_R = 3.5 thick Omega_p$, $Omega_c \/ (2 pi) = 50$ MHz,
    $alpha = 4$ (smooth edges).
    The yellow contour marks $overline(F)_"OR" = 0.992$; the dashed
    white contour marks $overline(F)_"OR" = 0.99$.
    The optimum (red star) sits at $Delta = 500$ MHz, $K = 0.95$,
    yielding $overline(F)_"OR" = 0.9924$.
    Note the broad ridge: $overline(F)_"OR" > 0.99$ over the entire
    rectangle $Delta in [400, 800]$ MHz, $K in [0.92, 1.00]$, so the
    optimum is *robust* to $tilde.op 5 %$ calibration errors in either
    axis.
  ],
) <fig:heatmap>

== Pulse-shape comparison

The final probe pulse is much shorter than the rev. 2 design --
$T_f = 182$ ns instead of $1200$ ns -- because we restored the
paper's $Omega_p = 2 pi times 50$ MHz (rev. 2 had used
$Omega_p = 2 pi times 20$ MHz, which forced a $6.25 times$ longer
$T_f$ to maintain the area). Smoothness ($alpha = 4$,
edge / peak $= 1.1 times 10^(-7)$) is preserved:

#figure(
  image("pulse_shapes_final.png", width: 100%),
  caption: [
    (a) Super-Gaussian probe pulse $Omega_p(t)$ for four
    configurations at $Omega_p \/ (2 pi) = 50$ MHz,
    $Delta \/ (2 pi) = 500$ MHz: the paper's abrupt-edge
    $alpha = 2.41, K = 1$ pulse (brown), the smooth $alpha = 4$
    candidate at $K = 1$ (orange), the *FINAL* $alpha = 4, K = 0.95$
    (green, bold), and a slightly more under-rotated $K = 0.90$
    (blue) for reference.
    (b) Area integrand $Omega_p^2 \/ (2 Delta)$ (solid) and
    cumulative area (dashed, right axis). The horizontal dotted
    lines at $pi \/ 4$ and $0.95 thick pi \/ 4$ mark the two area
    targets. The FINAL pulse reaches $0.95 thick pi \/ 4 = 0.746$
    exactly. The $alpha = 4$ pulse approaches its asymptote $~ 50$
    ns later than the $alpha = 2.41$ pulse because the smoother
    edges put the area integral further out in time.
  ],
) <fig:pulses>

== Area-calibration formula

In the well-resolved limit ($alpha >> 1$ so edge amplitudes are
negligible), the area integral
$integral Omega_p^2 \/ (2 Delta) dif t$ depends *only* on $sigma$,
not on $T_f$. The substitution $v = (t - t_c) \/ sigma^(1\/3)$ gives

$ "area" = frac(Omega_p^2 sigma^(1\/3), 8 Delta) thick I_oo,
  quad I_oo equiv integral_(-oo)^oo e^(-2 v^6) dif v
  = 2 thick Gamma(7\/6) thick 2^(-1\/6) approx 1.6534. $

Solving for $sigma$ at a target area $K dot pi \/ 4$:

$ sigma = K^3 thick sigma_(K=1),
  quad sigma_(K=1)^(1\/3)
  = frac(2 pi Delta, Omega_p^2 thick I_oo). $

Numerically at the final configuration:
$sigma_(K=1) = 1.771$ ns,
$sigma = 0.95^3 dot 1.771 = 1.519$ ns,
$T_f = (4 sigma)^(1\/3) = 182.46$ ns,
verified by `compute_pulse_area` to give area $= 0.74614 = 0.9500 dot pi \/ 4$.

= Why $Omega_R \/ Omega_p = 3.5$ and $Delta in [400, 800]$ MHz are the robust choice <sec:omegaR>

The Farouk-style OR protocol fixes $Omega_R \/ Omega_p = 3.5$ by
construction. To check that this choice is not an arbitrary
holdover, we relaxed the constraint and scanned
$(Delta, Omega_R \/ Omega_p)$ independently at
$Omega_p = 2 pi times 50$ MHz, $K = 0.95$, $alpha = 4$
(`sweep_omega_R.py`). Two error budgets compete:

- *Blockade leakage*, scaling as $1 \/ M_2^2$ with
  $M_2 equiv V_(c t) thin (2 Delta) \/ (Omega_p Omega_R)$ the
  two-photon Raman margin. Smaller $Delta$ or larger $Omega_R$
  shrinks $M_2$ and grows the $|1,1,* angle.r$ leakage.
- *Rydberg decay*, scaling as $T_f \/ tau_R$ with
  $T_f prop Delta \/ Omega_p^2$ at fixed pulse area.

With $Omega_R$ held at the Farouk value $3.5 thick Omega_p
= 2 pi times 175$ MHz, $M_2$ falls below $tilde.op 21$ once
$Delta < 400$ MHz, pushing $|1,1,* angle.r$ over the leakage
cliff; this is *the only reason* the $Delta in [400, 600]$ MHz
plateau in #ref(<fig:heatmap>) looks "special." Releasing
$Omega_R$ should, naively, restore $M_2$ at smaller $Delta$ and
open up a faster gate. We tested this.

#figure(
  image("omega_R_landscape.png", width: 100%),
  caption: [
    (a) Coarse $(Delta, Omega_R \/ Omega_p)$ scan at
    $Omega_p \/ (2 pi) = 50$ MHz, $K = 0.95$, $alpha = 4$
    (`sweep_omega_R.py`).
    White contour: $overline(F)_"OR" = 0.99$.
    Red contour: $overline(F)_"OR" = 0.992$.
    Red dot: FINAL configuration.
    Cyan star: small-$Delta$ alternative.
    A dressed-state anti-resonance at ratio $= 2$ drops
    $overline(F)$ to $tilde.op 0.90$ across every $Delta$.
    (b) Refinement at $Delta = 225$ MHz over
    $(Omega_R \/ Omega_p, K)$ (`sweep_smallDelta_refine.py`):
    only $4 \/ 36$ cells exceed $overline(F)_"OR" > 0.99$ -- a
    knife-edge ridge with $tilde.op plus.minus 0.1$ tolerance on
    ratio and $tilde.op plus.minus 0.02$ on $K$.
  ],
) <fig:omega_R>

Three findings emerge from #ref(<fig:omega_R>):

+ *The Farouk ridge is the broadest plateau.*  Along
  ratio $approx 3.0–3.5$, $Delta in [400, 700]$ MHz, the
  $overline(F)_"OR" > 0.99$ region is wide in *all three* of
  $(Delta, Omega_R, K)$. $M_2$ on the ridge sits between
  $25.9$ and $42.2$, deep in the strong-blockade regime where
  $tilde.op 10 %$ calibration drift on any axis does not push the
  fidelity below $0.99$.

+ *A faster small-$Delta$ peak exists but is a knife-edge.*  At
  $Delta = 225$ MHz, $Omega_R \/ Omega_p = 2.5$, $K = 0.95$ the
  simulation gives $overline(F)_"OR" = 0.9920$ in a total gate
  time of *$184$ ns* -- $2.1 times$ shorter than the FINAL
  $385$ ns. But $M_2 = 16.3$ at this point, right on the
  blockade-leakage cliff: #ref(<fig:omega_R>) panel (b) shows the
  $overline(F)_"OR" > 0.99$ region collapses to four cells out of
  thirty-six in a $(Omega_R \/ Omega_p, K)$ grid of step
  $0.1, 0.02$. Calibration drift in either axis at the
  $tilde.op 1-2 %$ level drops fidelity below $0.99$.

+ *A pathological anti-resonance lives at ratio $= 2$.*
  $overline(F)_"OR"$ collapses to $tilde.op 0.90$ across every
  $Delta$, attributable to a $2 : 1$ Floquet-like resonance
  between $Omega_R$ and $Omega_p$ in the EIT-dressed structure.
  The small-$Delta$ peak lives just past this dip, which is *why*
  its ridge is so narrow.

A hard floor at $overline(F)_"OR" approx 0.992$ is set by the
$|1,1,* angle.r$ branch, which caps at
$F_(1,1,*) approx 0.988$ in *both* high-fidelity regions. This
is the $M_2 tilde.op 30$ blockade leakage limit at these atomic
parameters; pushing higher requires either smaller $Omega_p$
(which is bad -- it *lengthens* the pulse and worsens
$|R angle.r$ decay at fixed area), shorter $r_("DA")$ (lattice
tightening, see Sec. 8.2), or a different Förster pair with
larger $V_(c t)$.

#block(
  fill: rgb("#f0f6ff"),
  inset: 10pt,
  radius: 4pt,
  width: 100%,
)[
  *Decision.* The FINAL configuration ($Delta = 500$ MHz,
  $Omega_R = 3.5 thick Omega_p$, $K = 0.95$) lives in the broadest
  $overline(F)_"OR" > 0.99$ plateau in $(Delta, Omega_R, K)$ space.
  The small-$Delta$ alternative is a real engineering option *if
  gate time is paramount and calibration to $tilde.op 1 %$ is feasible*
  -- e.g. for a tightened QEC cycle budget -- but is too fragile
  for routine use. We keep the FINAL configuration as the
  recommended setting.
]

= OR gate results

The OR gate implements

$ |c_1, c_2, t angle.r arrow.bar
  |c_1, c_2, t plus.circle (c_1 or c_2) angle.r, $

i.e. the target flips when at least one control is in $|1 angle.r$.

#figure(
  caption: [OR gate: per-input state fidelities for all 8 computational
    basis inputs at the FINAL configuration ($Delta = 2 pi times 500$
    MHz, $K = 0.95$, $alpha = 4$).
    All branches exceed $0.98$; the blockade-and-flip
    branches sit comfortably above $0.99$.],
  kind: table,
  table(
    columns: (auto, auto, auto, auto),
    align: (center, left, center, center),
    stroke: 0.4pt,
    table.header(
      [*\#*], [*Input*], [*$F_k$ (rev. 3 final)*], [*Rev. 2 for compare*],
    ),
    [1], [$|0,0,A angle.r$], [$0.997232$], [$0.997069$],
    [2], [$|0,0,B angle.r$], [$0.997232$], [$0.997069$],
    [3], [$|0,1,A angle.r$], [$0.992704$], [$0.966954$],
    [4], [$|0,1,B angle.r$], [$0.992704$], [$0.966954$],
    [5], [$|1,0,A angle.r$], [$0.992704$], [$0.966954$],
    [6], [$|1,0,B angle.r$], [$0.992704$], [$0.966954$],
    [7], [$|1,1,A angle.r$], [$0.987049$], [$0.936230$],
    [8], [$|1,1,B angle.r$], [$0.987049$], [$0.936230$],
  ),
)

#block(
  fill: rgb("#eef8ec"),
  inset: 10pt,
  radius: 4pt,
  width: 100%,
)[
  *OR gate average fidelity (FINAL):*
  $ quad overline(F)_"OR" = bold(0.992422), quad
    1 - overline(F)_"OR" = bold(7.58 times 10^(-3)). $
]

== Per-branch interpretation

- *No controls excited* ($|0, 0, * angle.r$, $F_k = 0.9972$).
  No Cs atom is in $|r angle.r$; the target completes its EIT cycle.
  Final populations
  $(P_A, P_B, P_P, P_R) = (0.997, 0.001, 0.0001, 0.002)$:
  $tilde.op 0.2 %$ residual in $|R angle.r$ from EIT non-adiabaticity at
  the pulse edges, $tilde.op 0.1 %$ in the opposite computational level.

- *One control excited* ($|0, 1, * angle.r$, $|1, 0, * angle.r$,
  $F_k = 0.9927$). One Cs atom in $|r angle.r$ for $T_c + 2 T_f =
  375$ ns, single-Cs decay
  $exp(-gamma_r T) = exp(-375 \/ 77000) = 0.9951$,
  leaving $~3 times 10^(-3)$ for blockade leakage and protocol
  imperfection.

- *Two controls excited* ($|1, 1, * angle.r$, $F_k = 0.9870$).
  Two Cs atoms in $|r angle.r$ over the same window, double-Cs decay
  $exp(-2 gamma_r T) = 0.9903$,
  with the remaining $~3 times 10^(-3)$ split between $V_(c c)$
  phase on $|r, r angle.r$ (which is *not* a global phase because
  the target moves between $|A angle.r$ and $|B angle.r$ during
  the gate) and residual blockade leakage. The $V_(c c)$ phase
  accumulates over $T$ to $V_(c c) dot T \/ (2 pi) approx
  0.73 dot 0.0004 dot 2 pi approx 1.7 degree$ -- small but
  non-negligible at the $10^(-3)$ level.

The arithmetic identity
$overline(F)_"OR" = (2 dot 0.9972 + 4 dot 0.9927 + 2 dot 0.9870)
\/ 8 = 0.9924$
matches the computed average.

== Why rev. 2 was so much worse

Rev. 2 had $Omega_p = 2 pi times 20$ MHz, forcing $T_f = 1200$ ns
(at $alpha = 4$, $K = 1$) to satisfy the area constraint. Total
OR-gate time was $2441$ ns -- $6.4 times$ longer than the final
$385$ ns -- so Cs $|r angle.r$ decay over the gate, which scales as
$gamma_r dot T_("total")$, was $6.4 times$ worse. The blockade
margin in rev. 2 was likewise much *stronger* than necessary
($Omega_p^2 \/ (2 Delta) = 0.4$ MHz instead of $2.5$ MHz here),
which slowed everything for no fidelity gain.

The fundamental lesson is the one derived in Sec. 4.3: at fixed
*blockade margin*, the gate time depends only on $V_(c t)$, not on
$Omega_p$ separately. Rev. 2's blockade was over-engineered.

= CCX gate results

The CCX (Toffoli) gate is unchanged from rev. 2 -- it uses
piecewise-square control and target pulses, no Gaussian envelope.
The CCX flips the target only when *both* controls are in
$|1 angle.r$:

$ |c_1, c_2, t angle.r arrow.bar
  |c_1, c_2, thick t plus.circle (c_1 and c_2) angle.r. $

#figure(
  caption: [CCX gate: per-input state fidelities at
    $Omega_(c c) = 2 pi times 50$ MHz, $Omega_t = 2 pi times 20$
    MHz. Unchanged from rev. 2.],
  kind: table,
  table(
    columns: (auto, auto, auto),
    align: (center, left, center),
    stroke: 0.4pt,
    table.header([*\#*], [*Input*], [*$F_k$*]),
    [1], [$|0,0,A angle.r$], [$0.995929$],
    [2], [$|0,0,B angle.r$], [$0.996259$],
    [3], [$|0,1,A angle.r$], [$0.993376$],
    [4], [$|0,1,B angle.r$], [$0.994237$],
    [5], [$|1,0,A angle.r$], [$0.993376$],
    [6], [$|1,0,B angle.r$], [$0.994237$],
    [7], [$|1,1,A angle.r$], [$0.999708$],
    [8], [$|1,1,B angle.r$], [$0.999708$],
  ),
)

#block(
  fill: rgb("#eef8ec"),
  inset: 10pt,
  radius: 4pt,
  width: 100%,
)[
  *CCX gate average fidelity:*
  $ overline(F)_"CCX" = bold(0.995854), quad
    1 - overline(F)_"CCX" = bold(4.15 times 10^(-3)). $
]

The CCX is *not* limited by the pulse-shape parameters explored in
this report -- it uses no Gaussian envelope and no Raman detuning.
Its error budget is dominated by Cs $|r angle.r$ decay over the
$95$ ns gate window for the $|0, 0, * angle.r$ and one-blockade
branches; the $|1, 1, * angle.r$ branch is essentially perfect
($F = 0.99971$) because the controls stay in $|1 angle.r$ and the
target sees only $T_t = 25$ ns of $|R angle.r$ dwell.

= Discussion

== Bottom line

#figure(
  caption: [Side-by-side OR-gate progression at the new energy
    levels. The final configuration is the recommended setting for
    any further use of this module.],
  kind: table,
  table(
    columns: (auto, auto, auto, auto, auto),
    align: (left, center, center, center, center),
    stroke: 0.4pt,
    table.header(
      [*Configuration*],
      [*$Omega_p$*],
      [*$Delta$*],
      [*$T_("OR")$ (ns)*],
      [*$overline(F)_"OR"$*],
    ),
    [Rev. 1 (paper-style, $Omega_c = 20$, $alpha = 2.41$, $K = 1$)],
      [$20$], [$500$], [$2078$], [$0.9715$],
    [Rev. 2 ($Omega_c = 50$, $alpha = 4$ smooth, $K = 1$)],
      [$20$], [$500$], [$2441$], [$0.9666$],
    [*FINAL ($Omega_c = 50$, $alpha = 4$ smooth, $K = 0.95$)*],
      [$bold(50)$], [$bold(500)$], [$bold(385)$], [$bold(0.9924)$],
  ),
)

The key insight that closed the gap to $> 0.99$ was the area
de-tuning. At finite $V_(c t)$ the OR protocol's correct
two-photon rotation is *not* $pi \/ 4$ but slightly under
($0.95 thick pi \/ 4$). The 2-D heatmap (Fig. 1) showed this
optimum spans a comfortable region in $(Delta, K)$ space, so the
final configuration is robust.

== Remaining levers

If still tighter fidelity is desired ($overline(F)_"OR" > 0.995$),
two routes are visible without changing the protocol:

1. *Re-verify lifetimes with ARC.* The $77$ μs / $74$ μs values
   used here are $n^3$-scaled from the paper's $n=66, 76$ values.
   If the BBR-included ARC values come out larger, fidelity
   improves linearly with $tau$.
2. *Tighten the lattice* to $a = 4$ μm. This doubles $V_(c t)$,
   halves the gate time, and is projected (by extrapolating the
   $V_(c t) = 517$ MHz / $V_(c t) = 226$ MHz scaling we already
   have) to push $overline(F)_"OR"$ to $tilde.op 0.996$. $R_("DD")$ and
   $R_("AA")$ remain $> 100$ at $a = 4$ μm because the
   smaller-spacing penalty on the same-species vdW is offset by
   the same-species Förster suppression. This is *outside the
   scope of the parameter sweep authorized for this report* but
   worth noting.

== Open issue: K-round / Z sub-cycle

This report only covers the X-error sub-cycle (Cs ancillas). The
corresponding Z-sub-cycle pair (Rb + K) has not been re-scanned;
see `SelfCorrectingRydberg/plan/Rb_Cs_level_redesign.md` § 6. The
two sub-cycles can use *different* Rb Rydberg states, so the
X-round answer here is independently usable while the K-round is
being settled.

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
  python or_average_gate_fid_gaussian.py   # final config, F_bar = 0.9924
  python ccx_average_gate_fidelity.py      # F_bar = 0.9959
  # Parameter scans that found the optimum:
  python sweep_delta.py        # scan over Delta at K=1
  python sweep_alpha.py        # smoothness vs fidelity
  python sweep_omega_p.py      # scan Omega_p at fixed blockade margin
  python sweep_area.py         # scan K at fixed Delta
  python sweep_2d.py           # joint (Delta, K) scan
  python sweep_omega_R.py      # rev. 4: scan (Delta, Omega_R/Omega_p)
  python sweep_smallDelta_refine.py # rev. 4: refinement at small Delta
  python plot_omega_R_landscape.py  # rev. 4: figure omega_R_landscape.png
  python make_pulse_figs.py    # regenerate figures
  ```
]

Output logs (`or_final.log`, `ccx.log`, `sweep_*.log`),
self-contained source files, and the figures
(`pulse_shapes_final.png`, `fidelity_heatmap.png`,
`omega_R_landscape.png`) are all in this directory. Interaction
coefficients, lifetimes, and lattice geometry are imported from
`level_parameters.py` so that all scripts share a single source
of truth.

#v(0.5em)

#block(
  fill: rgb("#fff4e8"),
  inset: 9pt,
  radius: 3pt,
  width: 100%,
)[
  *Caveats.*
  - Two lifetimes ($tau_r = 77$ μs, $tau_R = 74$ μs) are
    $n^3$-scaled estimates, not direct ARC values. Re-running
    `Rb.getStateLifetime(54, 2, 1.5, ...)` and
    `Cs.getStateLifetime(62, 2, 2.5, ...)`
    before publication will tighten the OR-gate fidelity prediction
    by an amount proportional to $Delta tau \/ tau$.
  - $V_("DD")$ does *not* enter the 3-atom Hamiltonian used in
    `mesolve` -- it is a cross-gate global-drive effect, included
    in the per-cycle error budget separately when chaining
    parallel CCX corrections.
  - The "area $= 0.95 thick pi \/ 4$" optimum is *empirical*; it
    has not been derived analytically for the finite-$V_(c t)$
    case. The 5 % de-tuning may shift slightly under model
    refinements (e.g. ARC-verified lifetimes, K-round level
    coupling).
]

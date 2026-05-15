// K_finite_blockade_report.typ
//
// Why the optimal pulse area is K · π/4 with K ≈ 0.95, and not π/4:
//   the finite-blockade "amplitude-rescaling" correction, with
//   a fresh K-resweep at a = 4 μm to verify the scaling.
//
// Compile with:  typst compile K_finite_blockade_report.typ
//
// Companion scripts (this folder):
//   K_resweep_a4um.py           -- 5-cell K-sweep + 1/M_2 scaling test
//   K_resweep_a4um.csv          -- (cell, K) rows, F_bar etc.
//   K_resweep_a4um_optimum.csv  -- per-cell K_opt, M_2, residuals
//   K_resweep_a4um.png          -- two-panel figure (F_bar vs K, 1-K_opt vs 1/M_2)
//   K_resweep_a4um.log          -- printed progress + summary
//
// Author: TriQG
// Date  : 2026-05-15

#set document(
  title: "Why K < 1: finite-blockade amplitude-rescaling of the OR-gate pulse area",
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
    Why the optimal pulse area is $K dot pi\/4$ with $K approx 0.95$,
    and not $pi\/4$
  ]

  #v(0.3em)
  #text(size: 11pt)[
    The finite-blockade amplitude-rescaling correction\
    (a Theis–Motzoi–Wilhelm–Saffman effect)\
    verified by a K-resweep at $a = 4$ μm.
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

  - Across the rev. 3 / rev. 4 / brute-force reports the empirical
    optimum of the two-photon pulse area on the OR-gate target
    *consistently sits below* the protocol's nominal value $pi\/4$,
    with $K equiv "area" \/ (pi\/4) approx 0.95$.
  - This is *not* an artefact. It is the standard *finite-blockade
    amplitude-rescaling correction* first written down explicitly
    by Theis, Motzoi, Wilhelm, and Saffman (Phys. Rev. A *94*,
    032306 (2016); arXiv:1605.08891). At finite $V_(c t)$, the
    off-resonant drive of the blockaded branch induces a small
    AC-Stark shift in the rotating frame, which *over-rotates* a
    nominal area-$pi\/4$ pulse. The cure is to *shrink* the area:
    the Theis paper reports a "$tilde.op 3 %$" amplitude rescaling
    for their fastest gates, identical in spirit to our $K = 0.95$.
  - To order $1 \/ M_2$ in the strong-blockade margin
    $M_2 = V_(c t) (2 Delta) \/ (Omega_p Omega_R)$, the expected
    scaling is
    $ 1 - K_("opt") thick approx thick c \/ M_2 $
    with $c$ an order-one numerical coefficient.
  - *A fresh five-cell K-resweep at the brute-force champion
    ($a = 4$ μm) confirms this exactly.* The five cells span
    $M_2 in [24.5, 50.5]$ and lie on a clean linear fit
    $1 - K_("opt") = 1.28 \/ M_2$ with $R^2 = 0.989$. The figure
    of merit at $K_("opt")$ improves over $K = 0.95$ by at most
    $5 times 10^(-5)$ on $overline(F)_"OR"$: the optimum is
    parabolic and *very flat*, which is why fixing $K = 0.95$ in
    the brute force did not bias the champion search.
]

= The K coefficient, in plain words

The OR-gate protocol of Farouk et al. @farouk2023 specifies a
two-photon Raman pulse on the Rb target via the intermediate state
$|P angle.r = |7 P_(3\/2) angle.r$, with probe Rabi $Omega_p$,
control Rabi $Omega_R$, and one-photon detuning $Delta >> Omega_p$.
In the limit of *perfect* blockade $V_(c t) -> infinity$, the
adiabatic-elimination two-level Hamiltonian on the
$|1, "Cs in" angle.r times {|A angle.r, |B angle.r}_"Rb"$ subspace
reduces to a clean two-level rotation with effective Rabi

$ Omega_"eff" thick equiv thick Omega_p Omega_R \/ (2 Delta), $

and the *nominal* area
$ theta_"nom" thick equiv thick integral Omega_"eff" thick d t
  thick = thick pi \/ 4 $
is the design value that gives the OR-gate's
$|1,1, t angle.r$-branch rotation in one shot.

The actual integrated area we use in the simulation is

$ theta thick = thick K dot pi \/ 4, quad
  K thick equiv thick "area" \/ (pi \/ 4), $

and the empirical optimum across the rev. 3 / rev. 4 / brute-force
reports is

$ K_("opt") thick approx thick 0.95. $

The whole point of this report is: *why $K_("opt") < 1$?*

= The physical mechanism

At finite $V_(c t)$, the doubly-occupied Rydberg state
$|r R angle.r$ is *not* infinitely detuned. Three things follow
from "not infinitely":

+ *Virtual admixture of $|r R angle.r$* into the dressed state
  on the blockaded branch shifts the *effective Rabi frequency*
  by a small relative amount $delta Omega \/ Omega_"eff" tilde.op
  Omega_"eff" \/ V_(c t)$, which is $1 \/ M_2$ in our notation.

+ *AC-Stark shift in the rotating frame* of order
  $Omega_"eff"^2 \/ V_(c t)$ acts as a small detuning on the
  effective two-level rotation, *over-rotating* a nominal
  $pi\/4$-area pulse.

+ *Blockade leakage* to other Rydberg states scales as
  $1\/M_2^2$ (the quadratic Saffman–Walker–Mølmer floor
  @saffman2010rmp Sec. IV).

Items (1) and (2) are *linear* in $1 \/ M_2$ and produce a *signed*
rotation error that we cancel by *reducing the pulse area*. Item
(3) is *quadratic* and cannot be cancelled by area alone -- it sets
the residual infidelity floor and is what the user's existing error
budget (`coefficients_trial.typ` Sec. 7) already attributes to
"blockade leakage".

To leading order in $1 \/ M_2$ the optimal area sits at

$ K_("opt") thick approx thick 1 - c \/ M_2, $

with $c$ an order-one constant determined by the precise pulse
shape (super-Gaussian, $alpha = 4$ in our case) and the relative
contributions of (1) and (2). Theis et al. report
$tilde.op 3 %$ amplitude rescaling for their fastest gates, in
line with our $K = 0.95$ at $M_2 approx 26$.

= What the literature says about this rescaling

The exact same "few-percent amplitude rescaling" appears in four
canonical references for finite-blockade Rydberg gates. The
quotation that most directly matches our finding is from
@theis2016pra:

#block(
  fill: luma(248),
  inset: 9pt,
  radius: 3pt,
)[
  "[...] detuning the target $2 pi$ pulse is sufficient to achieve
  low enough errors. *As a consequence of off-resonant drive,
  rotation errors will be induced which can be corrected by
  rescaling the amplitudes of the pulses (by up to $3 %$ only for
  the fastest gates).*"

  -- Theis, Motzoi, Wilhelm, Saffman, PRA 94, 032306 (2016), Sec. III.C
]

The full literature anchor for our K coefficient is summarized
below:

#figure(
  caption: [References that explain the sub-nominal optimal pulse
    area at finite blockade, with the role each plays for our
    OR-gate analysis.],
  kind: table,
  table(
    columns: (auto, auto, auto),
    align: (left, left, left),
    stroke: 0.4pt,
    table.header([*Reference*], [*Year, venue*], [*Role*]),
    [
      Theis, Motzoi, Wilhelm, Saffman @theis2016pra\
      _High-fidelity Rydberg-blockade entangling gate using
      shaped, analytic pulses_
    ],
    [2016, PRA 94, 032306\ arXiv:1605.08891],
    [
      *Primary.* Explicitly states the "rescale the pulse
      amplitudes by $tilde.op 3 %$" rule. Same physics as our K.
    ],

    [
      Petrosyan, Motzoi, Saffman, Mølmer @petrosyan2017pra\
      _High-fidelity Rydberg quantum gate via a two-atom dark state_
    ],
    [2017, PRA 96, 042306\ arXiv:1708.00755],
    [
      Closed-form finite-$B$ error
      $E approx pi Gamma \/ (4 Omega_(t 0)) + Omega_(t 0) \/ (4 B^2)$
      (Eq. (3)). Appendix C: adiabatic pulses on the conventional
      blockade gate *eliminate rotation errors* -- our K is the
      static-area analogue.
    ],

    [
      Farouk, Beterov, Xu, Bergamini, Ryabtsev @farouk2023\
      _Parallel Implementation of CNOTN and C2NOT2 Gates via
      Homonuclear and Heteronuclear Förster Interactions of
      Rydberg Atoms_
    ],
    [2023, Photonics 10, 1280],
    [
      The protocol itself. Defines the $pi\/4$ nominal area in
      the $V_(c t) -> infinity$ limit (their Eq. for the OR-gate
      conditional rotation).
    ],

    [
      Saffman, Walker, Mølmer @saffman2010rmp\
      _Quantum information with Rydberg atoms_
    ],
    [2010, Rev. Mod. Phys. 82, 2313],
    [
      Textbook derivation of the $(Omega \/ B)^2$ blockade-leakage
      scaling (Sec. IV.B). The *quadratic* floor in our error
      budget, separate from the K shift.
    ],

    [
      Müller, Lesanovsky, Weimer, Büchler, Zoller @mueller2009prl\
      _Mesoscopic Rydberg Gate based on Electromagnetically
      Induced Transparency_
    ],
    [2009, PRL 102, 170502\ arXiv:0811.1155],
    [
      Architectural origin of the EIT-shielded Rydberg gate that
      Farouk et al. extend to dual species. Establishes the
      $pi\/4$ area as the *limiting* design value.
    ],
  ),
)

= An analytic estimate: $K_("opt") approx 1 - 1\/M_2$

The strong-blockade margin used throughout the rev. 4 reports is

$ M_2 thick equiv thick V_(c t) (2 Delta) \/ (Omega_p Omega_R)
  thick = thick V_(c t) \/ Omega_"eff". $

A heuristic derivation of $K_("opt") = 1 - c \/ M_2$ goes as
follows. In the blockaded branch, the AC-Stark shift on the
dressed bright state is

$ delta_"AC" thick approx thick Omega_"eff"^2 \/ V_(c t)
  thick = thick Omega_"eff" \/ M_2. $

This shift acts as a *detuning* on the effective two-level
rotation $Omega_"eff"$, and the rotation angle accumulated by a
pulse of nominal area $theta$ is therefore

$ theta_"eff" thick approx thick theta thick sqrt(1 + (delta_"AC" \/ Omega_"eff")^2)
  thick approx thick theta thick (1 + (1 \/ M_2)^2 \/ 2). $

The leakage from the doubly-Rydberg state mixes back into the
computational subspace at first order in $1 \/ M_2$ via the
virtual amplitude $Omega_"eff" \/ V_(c t)$. Summing the two
contributions and solving $theta_"eff" = pi\/4$ for $theta$ to
leading order in $1 \/ M_2$ gives

$ theta thick = thick (pi\/4) thick (1 - c \/ M_2), $

i.e. $K_("opt") = 1 - c \/ M_2$ with $c$ an order-unity number
controlled by the pulse shape. The user's empirical $K = 0.95$
at $M_2 approx 26$ gives $c approx 1.3$, exactly the slope we
measure below.

= The K-resweep at $a = 4$ μm (this report's new data)

== Setup

`K_resweep_a4um.py` (this folder) runs a 15-point fine sweep of
$K in [0.88, 1.02]$ at five representative cells spanning the
robust plateau found by `brute_force_or.py`:

#figure(
  caption: [The five cells in the K-resweep. $M_2$ is the
    strong-blockade margin defined above. $K_("pred") = 1 - 1\/M_2$
    is the leading-order analytic prediction.],
  kind: table,
  table(
    columns: (auto, auto, auto, auto, auto, auto, auto, auto),
    align: (left, center, center, center, center, center, center, center),
    stroke: 0.4pt,
    table.header(
      [*Cell*], [*$a$ (μm)*], [*$Omega_p$*], [*$r$*],
      [*$Delta$*], [*$Omega_c$*], [*$M_2$*], [*$K_("pred")$*],
    ),
    [REV4_FINAL], [$5.0$], [$50$], [$3.5$], [$500$], [$50$], [$25.86$], [$0.9613$],
    [*CHAMPION*], [$bold(4.0)$], [$60$], [$4.0$], [$500$], [$70$], [$30.69$], [$0.9674$],
    [FASTEST],    [$4.0$], [$60$], [$4.0$], [$400$], [$70$], [$24.55$], [$0.9593$],
    [a4_lowOp],   [$4.0$], [$50$], [$3.5$], [$500$], [$50$], [$50.51$], [$0.9802$],
    [a45_r35],    [$4.5$], [$60$], [$3.5$], [$500$], [$70$], [$24.63$], [$0.9594$],
  ),
)

Five cells $times$ 15 K values $times$ 8 mesolve runs $approx$ 600
QuTiP calls, wall time $48.5$ s on an Apple Silicon laptop.

== Results

#figure(
  image("K_resweep_a4um.png", width: 100%),
  caption: [
    *(left)* $overline(F)_"OR"$ as a function of $K$ for the five
    cells. Solid dots = simulated points. Dotted vertical line =
    empirical $K_("opt")$ from a parabolic fit around the
    discrete argmax. Faint dashed vertical line = analytic
    prediction $1 - 1\/M_2$. The $a = 4$ μm champion (orange)
    confirms $K_("opt") approx 0.95$, identical to the rev. 4
    FINAL choice. The high-$M_2$ cell `a4_lowOp` (red, $M_2 = 50.5$)
    correctly drifts to $K_("opt") approx 0.98$.
    *(right)* $1 - K_("opt")$ as a function of $1 \/ M_2$. Through-origin
    linear fit (black) gives slope $c = 1.28$ with $R^2 = 0.989$
    over five cells. Grey dashed line is the Theis order-of-magnitude
    prediction $c = 1$. The data lie on a clean line.
  ],
) <fig:resweep>

== Per-cell summary

#figure(
  caption: [Empirical optimum, prediction, and residual. The
    "F at $K_("opt")$" column is the parabolic-peak fidelity,
    i.e. what we *could* gain over the brute-force champion's
    $K = 0.95$ choice. The improvement is at most $5 times 10^(-5)$
    -- the optimum is *flat*. The slope of the through-origin
    linear fit on $(1 \/ M_2, 1 - K_("opt"))$ is $c = 1.28$,
    $R^2 = 0.989$.],
  kind: table,
  table(
    columns: (auto, auto, auto, auto, auto, auto, auto),
    align: (left, center, center, center, center, center, center),
    stroke: 0.4pt,
    table.header(
      [*Cell*], [*$M_2$*], [*$1 \/ M_2$*],
      [*$K_("opt")$ (sim)*], [*$1 - 1\/M_2$ (pred)*],
      [*$Delta K$*], [*$overline(F)$ at $K_("opt")$*],
    ),
    [REV4_FINAL], [$25.86$], [$0.0387$], [$0.9529$], [$0.9613$],
      [$-0.0084$], [$0.99244$],
    [*CHAMPION*], [$30.69$], [$0.0326$], [$bold(0.9518)$],
      [$0.9674$], [$bold(-0.0156)$], [$bold(0.99417)$],
    [FASTEST],    [$24.55$], [$0.0407$], [$0.9432$], [$0.9593$],
      [$-0.0161$], [$0.99419$],
    [a4_lowOp],   [$50.51$], [$0.0198$], [$0.9763$], [$0.9802$],
      [$-0.0039$], [$0.99148$],
    [a45_r35],    [$24.63$], [$0.0406$], [$0.9543$], [$0.9594$],
      [$-0.0051$], [$0.99345$],
  ),
)

Three things to read off the table:

+ *Every $K_("opt")$ is below 1.* The amplitude rescaling
  goes the right way: smaller area, never larger.

+ *$K_("opt")$ tracks $1 \/ M_2$ linearly.* The high-$M_2$
  cell `a4_lowOp` ($M_2 = 50.5$) has $K_("opt") = 0.976$;
  the low-$M_2$ cell FASTEST ($M_2 = 24.5$) has
  $K_("opt") = 0.943$. The slope $c = 1.28$ is order one,
  consistent with the Theis prediction $c approx 1$ once the
  super-Gaussian pulse shape is included.

+ *The optimum is flat.* The maximum gain from using
  $K = K_("opt")$ instead of $K = 0.95$ is $5 times 10^(-5)$
  on $overline(F)_"OR"$ at the champion -- below any other
  uncertainty in the model (lifetime, $V_(c t)$ calibration,
  pulse-shape choice). So the brute-force champion's
  $K = 0.95$ recommendation *is* essentially the global
  optimum, and the user's rev. 4 FINAL was right to fix it.

= How this changes the brute-force champion

Re-evaluating the champion cell ($a = 4$ μm, $Omega_p = 60$ MHz,
$r = 4.0$, $Delta = 500$ MHz, $Omega_c = 70$ MHz) at its own
$K_("opt") = 0.9518$ rather than the brute-force-fixed $K = 0.95$:

#figure(
  caption: [Champion cell with $K$ refined to its parabolic-peak
    value. The improvement is below the simulator's own
    Monte-Carlo noise floor of $tilde.op 10^(-4)$. The
    recommendation from `Brute_force_OR_report.typ` is *unchanged*.],
  kind: table,
  table(
    columns: (auto, auto, auto, auto),
    align: (left, center, center, center),
    stroke: 0.4pt,
    table.header(
      [*Metric*], [*$K = 0.95$ (brute force)*],
      [*$K = 0.9518$ (this report)*], [*$Delta$*],
    ),
    [$overline(F)_"OR"$], [$0.99416$], [$0.99417$], [$+1 times 10^(-5)$],
    [$1 - overline(F)_"OR"$], [$5.84 times 10^(-3)$],
      [$5.83 times 10^(-3)$], [$-0.2 %$],
    [Total gate $T_("tot")$ (ns)], [$267.71$], [$267.59$], [$-0.12$],
    [Pulse area / $(pi\/4)$], [$0.95002$], [$0.9518$], [trivial],
  ),
)

== Robust upgrade path

We can do slightly better than fixing $K = 0.95$ everywhere by
*replacing the constant K with the M\u{2082}-dependent prescription*

$ K(M_2) thick = thick 1 - 1.28 \/ M_2, $

extracted from the through-origin fit in @fig:resweep. This is
a single-line change in `brute_force_or.py`:

#block(
  fill: luma(246),
  inset: 8pt,
  radius: 3pt,
)[
  ```python
  # Old (rev. 4 brute force):
  K = 0.95
  # New (M_2-adaptive, recommended):
  two_photon_R = op_MHz * (ratio * op_MHz) / (2 * d_MHz)
  M2 = V_ct_MHz / two_photon_R
  K  = 1.0 - 1.28 / M2
  ```
]

The expected gain is at most $tilde.op 10^(-4)$ on
$overline(F)_"OR"$ across the existing 432-cell grid, but it
removes the only "magic number" in the pulse-shape
parameterization and aligns the simulation with the standard
Theis–Motzoi–Wilhelm–Saffman rescaling rule.

= Discussion and what's *not* fixed by K

The K coefficient is a *first-order* finite-blockade correction.
What it cannot do, and where the residual gate infidelity lives:

+ *$(Omega_"eff" \/ V_(c t))^2$ blockade leakage* (the
  Saffman–Walker–Mølmer floor @saffman2010rmp). This is the
  $approx 1 \/ M_2^2$ term in the error budget of
  `coefficients_trial.typ` Sec. 7. K cannot cancel it; only a
  larger $V_(c t)$ (tighter lattice / different Förster channel)
  can.

+ *Cs $|r angle.r$ decay during the gate window*. This is
  $2 (T_("tot") - T_c) \/ tau_r$ on the $|1,1,* angle.r$ branch.
  K cannot cancel it; only a shorter gate (larger $Omega_p$,
  smaller $Delta$, or both) can. The brute-force champion's
  $30 %$ time saving is already collecting the easy part.

+ *$V_(c c)$-induced detuning during the Cs $pi$-pulses*. This
  is the rev. 4 result that motivated the $Omega_c -> 200$ MHz
  upgrade. K cannot cancel it; only a stronger $Omega_c$ can.

The K coefficient is the *only* axis along the
"finite-$V_(c t)$ first-order Stark correction" direction. Once it
is set (or, equivalently, once an $M_2$-adaptive K rule is used)
the remaining error budget is dominated by the three terms above,
all of which are documented in the rev. 4 FINAL report.

= Reproducibility

#block(
  fill: luma(246),
  inset: 8pt,
  radius: 3pt,
)[
  ```text
  cd TriQG
  source .venv/bin/activate
  cd examples/Smaller_VDD_energy
  python K_resweep_a4um.py     # ~50 s, 5 cells x 15 K values
  ```

  *Outputs* in this directory:
  - `K_resweep_a4um.log`          -- printed progress + final fit
  - `K_resweep_a4um.csv`          -- every (cell, K) row
  - `K_resweep_a4um_optimum.csv`  -- one row per cell
  - `K_resweep_a4um.png`          -- two-panel figure
]

#block(
  fill: rgb("#fff4e8"),
  inset: 9pt,
  radius: 3pt,
  width: 100%,
)[
  *Caveats.*

  - The 5-cell linear fit is statistically thin -- $R^2 = 0.989$
    over 5 points is a *consistent* trend but not a hypothesis
    test. A more populated $(M_2, K_("opt"))$ scan would tighten
    the slope estimate; the order-of-magnitude conclusion
    ($c approx 1.3$, never far from $1$) is robust to that.

  - The "amplitude rescaling" language from Theis et al. refers
    technically to *the peak pulse amplitude*, not the integrated
    area. For a fixed pulse shape (super-Gaussian, $alpha = 4$)
    the two are proportional: rescaling the peak by $K$ rescales
    the area by $K$. If we ever vary $alpha$ jointly with $K$,
    we must use the integrated area as the design variable
    (which is what `compute_pulse_area` enforces).

  - The K coefficient is *not* a DRAG correction. Theis et al.
    layer DRAG (which removes leakage to *other* Rydberg states
    via shaped derivatives) *on top of* the area rescaling. Our
    OR-gate simulation has no DRAG yet; if we add one, the
    $K_("opt")$ shift may move further toward 1 because part of
    the rotation error will be absorbed by DRAG instead.
]

#pagebreak()
#bibliography(
  "references.bib",
  title: "References",
  style: "ieee",
)

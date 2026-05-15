// coefficients_trial.typ
//   Companion to Smaller_VDD_average_fidelity_report.typ.
//
//   Documents a series of coefficient / parameter trials run after the
//   FINAL (rev. 3) report:
//     - Sweep Omega_R independently of the Farouk 3.5 ratio
//     - Test the small-Omega_p regime
//     - ARC-verify lifetimes and scan Cs ancilla alternatives
//     - Tighten lattice to a = 4 um
//     - Diagnose the V_cc-during-pi-pulse error and propose Rev4
//
//   Compile with:  typst compile coefficients_trial.typ
//
//   Scripts (all in examples/Smaller_VDD_energy/):
//     sweep_omega_R.py, sweep_smallDelta_refine.py, plot_omega_R_landscape.py
//     sweep_smallOp.py
//     or_a4um.py, or_a4um_arclife.py, or_a4um_revcandidate2.py
//   And ARC-side scan in SelfCorrectingRydberg/python_scripts/:
//     scan_cs_alternatives.py  (-> proofread/appendix-a/cs_alternatives_scan.csv)

#set document(
  title: "Coefficients trial: Omega_R, lattice, and V_cc-pi error budget",
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
    Coefficients trial: $Omega_R$, lattice, and the $V_(c c)$-$pi$
    error budget
  ]

  #v(0.3em)
  #text(size: 11pt)[
    Companion to `Smaller_VDD_average_fidelity_report.typ`.  Documents
    the parameter trials run after rev. 3 FINAL and the Rev4 candidate.
  ]

  #v(0.4em)
  #text(size: 10pt)[
    TriQG `examples/Smaller_VDD_energy` · compiled 2026-05-14
  ]
]

#v(0.6em)

#block(
  fill: rgb("#eef8ec"),
  inset: 10pt,
  radius: 4pt,
  width: 100%,
)[
  *Headline.*  After running six parameter trials beyond the FINAL config:

  - The Farouk ratio $Omega_R \/ Omega_p = 3.5$ at $Delta = 500$ MHz
    sits on the *broadest* $overline(F)_"OR" > 0.99$ plateau.  A
    secondary peak at $(Delta = 225$ MHz, ratio $= 2.5)$ matches the
    fidelity in $2.1 times$ shorter gate but is a knife-edge.
  - $Omega_p = 20$ MHz at $Delta = 200$ MHz cannot reach $0.99$ in any
    $(Omega_R, K)$ corner -- gate stretches as $1\/Omega_p^2$, decay
    dominates.
  - ARC corrects $tau_r$ from $77 mu s$ ($n^3$-scaled) to $bold(83.9 mu s)$.
    No other Cs state simultaneously matches the F\"orster channel
    and $R_("AA") >= 100$: the present Cs $62 D_(5\/2)$ is uniquely pinned.
  - Lattice tightening to $a = 4$ μm doubles $V_(c t)$ but quadruples
    $V_(c c)$.  Naive transplant (FINAL config $->$ a = 4) loses
    $0.6 %$ on $|1,1,* angle.r$.  Diagnosis: $V_(c c)$ detunes
    $|r,r angle.r$ during the Cs $pi$-pulses, error $tilde.op 4 (V_(c c)\/Omega_c)^2$.
  - *Rev4 candidate:* $a = 4$ μm, $Omega_p = 50$, $Omega_R = 175$,
    $Delta = 500$ MHz, $K = 0.98$, *$Omega_c = 200$ MHz* gives
    $overline(F)_"OR" = bold(0.99357)$ vs FINAL's $0.99282$
    (with ARC lifetimes throughout).  Gate $381$ ns vs FINAL $385$ ns.
    A marginal $0.075 %$ gain at essentially the same gate time --
    not the $0.996$ that lattice tightening alone seemed to promise.
]

= Why $Delta in [400, 600]$ MHz with $Omega_R = 3.5 thick Omega_p$

The OR-gate has two competing errors:

#figure(
  caption: [Error budget at the FINAL config ($Omega_p = 50$,
    $Omega_R = 175$, $Omega_c = 50$ MHz, $K = 0.95$, $alpha = 4$,
    $a = 5$ μm).  $M_2 = V_(c t) thin (2 Delta)\/(Omega_p Omega_R)$
    is the two-photon blockade margin.],
  kind: table,
  table(
    columns: (auto, auto, auto, auto),
    align: (left, center, center, center),
    stroke: 0.4pt,
    table.header(
      [*$Delta\/(2 pi)$ MHz*], [*$M_2$*], [*Leakage $tilde.op 1\/M_2^2$*],
      [*Decay $2 T_g\/tau_r$*],
    ),
    [300], [15.5], [0.42 %], [0.27 %],
    [400], [20.7], [0.23 %], [0.37 %],
    [*500 (FINAL)*], [25.9], [0.15 %], [*0.46 %*],
    [600], [31.0], [0.10 %], [0.55 %],
    [800], [41.4], [0.058 %], [0.73 %],
    [1000], [51.7], [0.037 %], [0.92 %],
  ),
)

The plateau $Delta in [400, 800]$ MHz balances the two.  Below $400$,
leakage takes off; above $800$, decay catches up.  With $tau_r$ updated
to $83.9 mu s$, the decay column is slightly tighter than what the
FINAL report quoted, but the plateau structure is unchanged.

= Sweeping $Omega_R$ independently of the Farouk ratio

`sweep_omega_R.py` scans $(Delta, Omega_R\/Omega_p)$ over
$Delta in [200, 1000]$ MHz and ratio $in [1.0, 5.0]$ at fixed
$Omega_p = 50$ MHz, $K = 0.95$, $alpha = 4$.

#figure(
  image("omega_R_landscape.png", width: 100%),
  caption: [
    (a) Coarse $(Delta, Omega_R \/ Omega_p)$ scan at $Omega_p\/(2 pi)
    = 50$ MHz, $K = 0.95$.  White contour: $overline(F)_"OR" = 0.99$.
    Red contour: $overline(F)_"OR" = 0.992$.  Red dot: FINAL config.
    Cyan star: small-$Delta$ alternative.  A dressed-state
    anti-resonance at ratio $= 2$ drops $overline(F)$ to
    $tilde.op 0.90$ across every $Delta$.
    (b) Refinement at $Delta = 225$ MHz showing the
    $(Omega_R\/Omega_p, K)$ landscape near the small-$Delta$ peak --
    knife-edge ridge: only $4\/36$ cells with $overline(F)_"OR" > 0.99$.
  ],
) <fig:omega_R>

Three findings:

+ The Farouk ridge (ratio $approx 3.0-3.5$, $Delta in [400, 700]$ MHz)
  is the *broadest* $overline(F)_"OR" > 0.99$ plateau.
+ A *second* peak exists at $Delta = 225$ MHz, ratio $= 2.5$,
  $K = 0.95$ giving $overline(F)_"OR" = 0.9920$ in a *$184$ ns*
  gate -- $2.1 times$ faster than FINAL.  But it is a knife-edge:
  only $4$ of $36$ cells in a $(Omega_R\/Omega_p, K)$ refinement
  grid exceed $0.99$, with $tilde.op plus.minus 0.1$ tolerance on
  ratio and $tilde.op plus.minus 0.02$ on $K$.
+ A pathological dip at ratio $= 2$ collapses $overline(F)$ to
  $tilde.op 0.90$ across every $Delta$.  Attributed to a $2:1$
  Floquet-like resonance in the EIT dressed-state structure.

= Small $Omega_p$ trial: the gate-time penalty

If $Omega_p$ is dropped from $50$ to $20$ MHz to make adiabatic
elimination cleaner, the area constraint forces

$ T_f prop Delta \/ Omega_p^2 thick arrow.r thick T_f(Omega_p = 20)
  = 6.25 thick T_f(Omega_p = 50) $

at the same $Delta$.  `sweep_smallOp.py` scans
$Delta in [100, 1000]$ MHz $times$ ratio $in [2.0, 6.0]$ at
$Omega_p = 20$ MHz.  *Zero of 42 cells reach $overline(F)_"OR" > 0.99$.*

#figure(
  caption: [Best cells in the $Omega_p = 20$ MHz scan.  Gate time is
    multiples of $T_f$ which scales as $1 \/ Omega_p^2$.  All cells
    are decay-limited.],
  kind: table,
  table(
    columns: (auto, auto, auto, auto, auto, auto),
    align: (center, center, center, center, center, center),
    stroke: 0.4pt,
    table.header(
      [*$Delta$ MHz*], [*ratio*], [*gate ns*], [*$M_2$*],
      [*decay floor*], [*$overline(F)$*],
    ),
    [150], [4.50], [704], [37.7], [1.83 %], [*0.9857*],
    [200], [3.50], [932], [64.7], [2.42 %], [0.9819],
    [300], [4.50], [1388], [75.4], [3.60 %], [0.9784],
    [500], [4.50], [2301], [125.7], [5.97 %], [0.9670],
  ),
)

The dimensionless parameter controlling the $|1,1,* angle.r$ floor
is $2 T_g \/ tau_r$, not $Omega_p\/(2 Delta)$.  At
$Omega_p = 20$ MHz, $Delta = 200$ MHz, the $700-900$ ns gate alone
costs $tilde.op 2 %$ on $|1,1,* angle.r$ from Cs decay.  *No setting
of $Omega_R$ or $K$ rescues this:* blockade can be made arbitrarily
deep, but the ancilla still decays during the long dwell.

Optimal $Omega_p$ in this family is *$25-50$ MHz*, where the
gate is fast enough that $T_g \/ tau_r$ does not bind.

= ARC lifetime verification and Cs ancilla scan

== Lifetime values

$n^3$-scaled estimates were off by $tilde.op 10 %$:

#figure(
  kind: table,
  table(
    columns: (auto, auto, auto),
    align: (left, center, center),
    stroke: 0.4pt,
    table.header([*State*], [*$n^3$ estimate*], [*ARC ($300$ K, BBR)*]),
    [Cs $62 D_(5\/2)$ ($|r angle.r$)], [$77$ μs], [*$bold(83.9 mu s)$*],
    [Rb $54 D_(3\/2)$ ($|R angle.r$)], [$74$ μs], [*$bold(83.4 mu s)$*],
    [Rb $7 P_(3\/2)$ ($|P angle.r$)], [$0.131$ μs], [$0.131$ μs (unchanged)],
  ),
)

Updating the simulation lifetimes lifts FINAL's $overline(F)_"OR"$
from $0.9924$ to $0.99282$.

== Scan for longer-lived Cs ancillas

`scan_cs_alternatives.py` enumerates Cs $n D_J$ for
$n in [50, 95]$, $J in {3\/2, 5\/2}$, with Rb $54 D_(3\/2)$ fixed.
For each candidate it computes (i) BBR-included lifetime, (ii)
F\"orster defect for the canonical channel
$ |"Rb" thick 54 D_(3\/2);"Cs" thick n D_J angle.r
  arrow.r |"Rb" thick 55 P_(3\/2);"Cs" thick (n - 2) F_(J')angle.r, $
(iii) $tilde C_3^"Rb-Cs"$, (iv) Cs--Cs $C_6$, (v) $R_("AA")$ at
$r_("AA") = 7.07$ μm.

Filters: $tau_r >= 80 mu s$, $|Delta_F| <= 200$ MHz, $R_("AA") >= 30$.

#figure(
  caption: [The F\"orster defect is monotonic in $n$.  Cs $62 D_(5\/2)$
    sits at a single zero-crossing for the Rb $54 D_(3\/2)$ channel.
    Lifetime grows as $n^3$, $|C_6^"CsCs"|$ as $tilde.op n^(11)$, so
    $R_("AA")$ falls faster than $tau$ rises.  *No alternative survives
    all three filters.*],
  kind: table,
  table(
    columns: (auto, auto, auto, auto, auto, auto),
    align: (left, center, center, center, center, center),
    stroke: 0.4pt,
    table.header(
      [*Cs state*], [*$tau_r$ μs*], [*$|Delta_F|$ MHz*],
      [*$V_(c t)$ MHz*], [*$R_("AA")$*], [*verdict*],
    ),
    [$60 D_(5\/2)$], [$76.9$], [$1431$], [$211$], [$422$], [F\"orster fail],
    [*$62 D_(5\/2)$ (current)*], [*$83.9$*], [*$5.2$*], [*$226$*], [*$309$*],
    [*pass*],
    [$63 D_(5\/2)$], [$87.5$], [$654$], [$233$], [$266$], [F\"orster fail],
    [$70 D_(5\/2)$], [$115.5$], [$4200$], [$290$], [$96$], [both fail],
    [$80 D_(5\/2)$], [$163.3$], [$7299$], [$381$], [$27$], [both fail],
    [$95 D_(5\/2)$], [$253.4$], [$9790$], [$541$], [$5.5$], [both fail],
  ),
)

*Conclusion.*  With Rb $54 D_(3\/2)$ fixed and the canonical D--D $arrow.r$
P--F F\"orster channel, Cs $62 D_(5\/2)$ is the unique near-resonant
ancilla.  Longer-lived candidates require either a different F\"orster
channel (different $Delta n$ or partner $L$), an F-state Cs ancilla,
or a different Rb data state -- each of which decouples from the
current redesign.

= Lattice tightening: a = 4 μm

Atomic constants $tilde C_3$, $C_6^"RbRb"$, $C_6^"CsCs"$ are
$a$-independent; only the geometric V's rescale.

#figure(
  caption: [Interaction strengths and selectivity ratios across $a$.
    Both ratios stay above $100$ down to $a = 4$ μm.],
  kind: table,
  table(
    columns: (auto, auto, auto, auto, auto, auto),
    align: (center, center, center, center, center, center),
    stroke: 0.4pt,
    table.header(
      [*$a$ (μm)*], [*$V_(c t)\/(2 pi)$*], [*$V_("DD")\/(2 pi)$*],
      [*$V_(c c)\/(2 pi)$*], [*$R_("DD")$*], [*$R_("AA")$*],
    ),
    [$3.50$], [$659.7$ MHz], [$4.934$ MHz], [$-6.188$ MHz], [$133.7$], [$106.6$],
    [$4.00$], [$bold(441.9)$ MHz], [$2.214$ MHz], [$-2.777$ MHz], [$199.6$], [$159.1$],
    [$4.50$], [$310.4$ MHz], [$1.092$ MHz], [$-1.370$ MHz], [$284.2$], [$226.6$],
    [*$5.00$ FINAL*], [*$226.3$ MHz*], [$0.581$ MHz], [$-0.728$ MHz],
    [*$389.8$*], [*$310.8$*],
  ),
)

Going from $a = 5 -> 4$ μm: $V_(c t)$ grows $1.95 times$, but $V_(c c)$
grows $3.81 times$ (steeper $1\/r^6$ scaling).  Both $R$ ratios halve.

== Naive transplant fails

The "obvious" recipe -- transplant FINAL config to $a = 4$, retune $Delta$
to keep $M_2$ constant -- *does not improve $overline(F)$*.  Three
trials at $a = 4$ μm with ARC lifetimes:

#figure(
  caption: [Trials at $a = 4$ μm with K refined.  Same $Omega_p$,
    $Omega_R$, $Omega_c = 50$ MHz throughout.],
  kind: table,
  table(
    columns: (auto, auto, auto, auto, auto, auto, auto),
    align: (left, center, center, center, center, center, center),
    stroke: 0.4pt,
    table.header(
      [*Case*], [*$Delta$ MHz*], [*$K$*], [*gate ns*],
      [*$M_2$*], [*$overline(F)$*], [*$F_(1\,1\,*)$*],
    ),
    [REF ($a = 5$ FINAL)], [$500$], [$0.95$], [$385$], [$25.9$],
    [*$0.99282$*], [$0.9878$],
    [B: same $Delta$], [$500$], [$0.98$], [$396$], [$50.5$], [$0.99188$], [$0.9819$],
    [C: $M_2$-matched], [$250$], [$0.97$], [$206$], [$25.3$], [$0.99122$], [$0.9902$],
  ),
)

Two distinct things go wrong:

- *Case B* (same $Delta$) loses $0.6 %$ on $|1,1,* angle.r$.  $M_2$
  is *twice* as big as FINAL, so it cannot be blockade leakage.
- *Case C* (smaller $Delta$ to keep gate short) gains on
  $|1,1,* angle.r$ but loses $0.8 %$ on $|0,0,* angle.r$.

== Diagnosis: $V_(c c)$ during the Cs $pi$-pulses

$V_(c c)$ as a static interaction on $|r, r angle.r$ is a *global
phase* for computational basis inputs.  It does *not* show up in
state fidelity for separable inputs (this is what the Hamiltonian
docstring asserts; we verified by hand for Case B).

But during the Cs $pi$-pulses at the gate boundaries, *both* controls
for a $|1, 1, t angle.r$ input rotate $|1 angle.r -> |r angle.r$
simultaneously.  The intermediate $|r, r angle.r$ sees an extra
energy $V_(c c)$ that the laser drive does not compensate, so the
$pi$-rotation on the $|r, r angle.r$ component is detuned.  Per-pulse
imperfect rotation $tilde.op (V_(c c) \/ Omega_c)^2$; with two $pi$-pulses
and two controls $tilde.op 4 (V_(c c) \/ Omega_c)^2$ total on $|1, 1, * angle.r$.

#figure(
  caption: [Predicted vs simulated $V_(c c)$-induced infidelity on
    $|1,1,*angle.r$ at $a = 4$ μm, $Delta = 500$ MHz, $K = 0.98$.
    Decay floor is $0.95 %$.  Residual = simulated infidelity minus
    decay floor.  The $(V_(c c)\/Omega_c)^2$ prediction tracks the
    simulated residual to a factor of $tilde.op 1.5$.],
  kind: table,
  table(
    columns: (auto, auto, auto, auto, auto, auto),
    align: (center, center, center, center, center, center),
    stroke: 0.4pt,
    table.header(
      [*$Omega_c$ (MHz)*], [*$T_c$ ns*], [*$4(V_(c c)\/Omega_c)^2$*],
      [*$1 - F_(1\,1\,*)$ simul*],
      [*residual*], [*$overline(F)$*],
    ),
    [$50$], [$10.0$], [$1.24 %$], [$1.81 %$], [$0.86 %$], [$0.99188$],
    [$75$], [$6.67$], [$0.55 %$], [$1.42 %$], [$0.47 %$], [$0.99287$],
    [$100$], [$5.00$], [$0.31 %$], [$1.28 %$], [$0.33 %$], [$0.99322$],
    [$150$], [$3.33$], [$0.14 %$], [$1.18 %$], [$0.23 %$], [$0.99348$],
    [*$200$*], [*$2.50$*], [*$0.08 %$*], [*$1.14 %$*], [*$0.19 %$*],
    [*$0.99357$*],
  ),
)

The $|0, 0, * angle.r$ branch is *completely untouched* by $Omega_c$
($F = 0.9992$ in all five rows): no Cs in $|r angle.r$, no $V_(c c)$
effect.  This confirms the diagnosis.

= Rev4 candidate

Combining ARC lifetimes, $a = 4$ μm, $Delta = 500$ MHz (preserving
$|0, 0, * angle.r$ EIT), $K = 0.98$ (sub-$pi\/4$ optimum at the new
$V_(c t)$), $Omega_c = 200$ MHz (suppressing $V_(c c)$-$pi$ error):

#figure(
  caption: [Rev4 candidate side-by-side with FINAL.  Same $Omega_p =
    50$, $Omega_R = 175$ MHz; ARC lifetimes throughout.  The "gain"
    column is the fractional improvement on each metric.],
  kind: table,
  table(
    columns: (auto, auto, auto, auto),
    align: (left, center, center, center),
    stroke: 0.4pt,
    table.header(
      [*Metric*], [*FINAL (rev. 3)*], [*Rev4 candidate*], [*$Delta$*],
    ),
    [$a$ (μm)], [$5.00$], [$4.00$], [$-1.0$],
    [$V_(c t)\/(2 pi)$ (MHz)], [$226$], [$442$], [$+1.95 times$],
    [$Omega_c\/(2 pi)$ (MHz)], [$50$], [$bold(200)$], [$+4 times$],
    [$Delta\/(2 pi)$ (MHz)], [$500$], [$500$], [unchanged],
    [$K$], [$0.95$], [$0.98$], [$+0.03$],
    [$T_c$ (ns)], [$10.0$], [$2.5$], [$0.25 times$],
    [Total gate (ns)], [$385$], [$381$], [$0.99 times$],
    [$M_2$], [$25.9$], [$50.5$], [$+1.95 times$],
    [$overline(F)_"OR"$], [$0.99282$], [*$bold(0.99357)$*], [$+0.075 %$],
    [$1 - overline(F)_"OR"$], [$7.18 times 10^(-3)$],
    [$6.43 times 10^(-3)$], [$0.89 times$],
  ),
)

Per-input fidelities at Rev4:
$F_(0\,0\,*) = 0.9992$,
$F_(0\,1\,*) = F_(1\,0\,*) = 0.9933$,
$F_(1\,1\,*) = 0.9886$.

== What Rev4 actually buys you

Two things, both modest:

- *$overline(F)_"OR"$ ceiling* rises from $0.99282$ to $0.99357$.  Real,
  but $0.075 %$ -- well below the $0.996$ that lattice tightening
  alone seemed to promise before the $V_(c c)$-$pi$ error was
  understood.
- *Gate time* does *not* speed up.  The original hope ("$3 times$
  faster at same $overline(F)$") died because $Delta$ must stay
  $tilde.op 500$ MHz to preserve $|0, 0, * angle.r$ EIT.  We could
  only collect the $V_(c c)$-$pi$ savings, not the decay savings.

What Rev4 *costs* is a $4 times$ stronger Cs control laser
($2 pi times 200$ MHz vs the FINAL's $50$).  At $Omega_c\/(2 pi) =
200$ MHz, off-resonant excitation of nearby Cs Rydberg states
(fine-structure partner $62 D_(3\/2)$, neighboring $n D$) should be
re-checked at the $10^(-4)$ level before committing.

= Error budget formula

The simulations make the OR-gate error budget legible:

$ 1 - overline(F)_"OR" thick approx thick
  underbrace(frac(2 T_g, tau_r), "decay floor on" |1\,1\,*angle.r)
  thin + thin
  underbrace(frac(1, M_2^2), "blockade leakage")
  thin + thin
  underbrace(4 thick (V_(c c)\/Omega_c)^2, "V"_(c c)*"-"*pi" on Cs ctrls")
  thin + thin
  underbrace((Omega_R^2 \/ (4 Delta^2))^"-ish", "EIT residue on" |0\,0\,*angle.r). $

Each term has its own knob:

#figure(
  kind: table,
  table(
    columns: (auto, auto, auto),
    align: (left, left, left),
    stroke: 0.4pt,
    table.header([*Term*], [*Reduce by*], [*FINAL value*]),
    [decay], [larger $V_(c t)$ (smaller $Delta$, smaller $T_f$), or
              longer $tau_r$], [$0.92 %$],
    [blockade], [larger $M_2$, i.e. larger $V_(c t)$ or larger $Delta$],
                [$0.15 %$],
    [$V_(c c)$-$pi$], [larger $Omega_c$ (shorter $T_c$)],
                     [$0.09 %$ (small at $a = 5$)],
    [EIT residue], [larger $Delta$ (so $Omega_R \/ Delta$ smaller)],
                   [$0.30 %$],
  ),
)

The knobs partially conflict: smaller $Delta$ shrinks decay AND
$M_2$ but enlarges EIT residue.  Larger $V_(c t)$ shrinks decay AND
blockade but enlarges $V_(c c)$ (faster) so $V_(c c)$-$pi$ grows.
The Rev4 candidate exploits the one knob (larger $Omega_c$) that
has no partner.

= Reproducibility

#block(
  fill: luma(246),
  inset: 8pt,
  radius: 3pt,
)[
  ```text
  cd TriQG/examples/Smaller_VDD_energy
  source ../../.venv/bin/activate
  # Omega_R sweep + small-Delta refinement
  python sweep_omega_R.py
  python sweep_smallDelta_refine.py
  python plot_omega_R_landscape.py     # writes omega_R_landscape.png
  # Small-Omega_p trial
  python sweep_smallOp.py
  # a = 4 um trials
  python or_a4um.py
  python or_a4um_arclife.py
  python or_a4um_revcandidate2.py
  # ARC ancilla scan
  cd ../../../SelfCorrectingRydberg
  ../TriQG/.venv/bin/python python_scripts/scan_cs_alternatives.py
  ```

  Output logs `sweep_*.log`, `or_a4um*.log`, and the CSV
  `proofread/appendix-a/cs_alternatives_scan.csv` document each trial.
]

#v(0.5em)

#block(
  fill: rgb("#fff4e8"),
  inset: 9pt,
  radius: 3pt,
  width: 100%,
)[
  *Caveats and open work.*

  - $Omega_c\/(2 pi) = 200$ MHz: re-check off-resonant excitation of
    Cs $62 D_(3\/2)$ (fine-structure partner) and nearest D, S
    states before committing.  If the splitting is $< 1$ GHz, a
    composite (BB1) $pi$-pulse at $Omega_c = 100$ MHz may be safer
    than raw $200$ MHz.
  - The $(V_(c c)\/Omega_c)^2$ prediction agrees with simulation
    within a factor of $tilde.op 1.5$ -- there is residual physics
    (precise pulse shape, transient $|r, 1 angle.r$, $|1, r angle.r$
    populations) not captured by the leading-order estimate.
  - F-state Cs ancillas and alternative F\"orster channels
    ($Delta n_("Cs") in {-1, 0, +1}$, partner $L in {2, 4}$) were
    *not* scanned -- they remain a route to longer $tau_r$
    without abandoning the Rb $54 D_(3\/2)$ data state.
  - Spin-echo schemes that decouple $V_(c c)$ on $|r, r angle.r$
    during the target Raman were *not* implemented; for
    computational basis inputs they reduce to a no-op, but they
    matter for process fidelity and superposition inputs.
]

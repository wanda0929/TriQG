// Brute_force_OR_K1_addendum.typ
//   Addendum to Brute_force_OR_report.typ.
//   Same 5-D sweep but with K = 1 (canonical pi/4 area) instead of K = 0.95.
//
// Compile with:  typst compile Brute_force_OR_K1_addendum.typ
//
// Scripts:
//   brute_force_or_K1.py             (5-D sweep, K = 1)
//   brute_force_K_comparison_plot.py (K = 0.95 vs K = 1 figure)
//
// Outputs / inputs:
//   brute_force_results_K1.csv, brute_force_robust_K1.csv,
//   brute_force_or_K1.log,
//   brute_force_K_comparison.png

#set document(
  title: "Addendum: brute-force OR sweep at K = 1",
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
    Addendum: brute-force OR sweep at $K = 1$
  ]

  #v(0.25em)
  #text(size: 10.5pt)[
    Does the canonical $pi \/ 4$ area $K = 1$ admit any
    $overline(F)_"OR" > 0.99$ cells?
  ]

  #v(0.3em)
  #text(size: 9pt)[
    TriQG `examples/Smaller_VDD_energy` · companion to
    `Brute_force_OR_report.typ` · 2026-05-15
  ]
]

#v(0.4em)

#block(
  fill: rgb("#eef8ec"),
  inset: 10pt,
  radius: 4pt,
  width: 100%,
)[
  *Answer:* Yes, but the high-fidelity plateau is *3.6 times smaller*
  than at $K = 0.95$ and the peak fidelity is $0.0013$ lower.

  - *$K = 1$ pass count*: $bold(61 \/ 432)$ cells ($14 %$) clear
    $overline(F)_"OR" > 0.99$, vs $221 \/ 432$ ($51 %$) at $K = 0.95$.
  - *$K = 1$ champion*: $a = 4.0$ μm, $Omega_p = 2 pi times 60$ MHz,
    $bold(Omega_R \/ Omega_p = 3.0)$, $Delta = 2 pi times 400$ MHz,
    $Omega_c = 2 pi times 70$ MHz $arrow.r overline(F)_"OR" = bold(0.99280)$
    in $bold(228) thick "ns"$, also fully robust ($5 \/ 5$ neighbors).
  - *No $K = 1$ cell passes at $a >= 5.0$ μm* -- a plateau that was
    half-full at $K = 0.95$ goes empty at $K = 1$. So the rev. 4
    setting ($a = 5.0$, $Omega_p = 50$, ratio $= 3.5$, $Delta = 500$)
    *does need* $K < 1$ to clear $0.99$.
  - *Lesson*: the optimal $K$ is not a property of the protocol; it
    is a *compensation knob* whose value depends on the other five
    knobs. $K = 0.95$ helps when $a$ is loose; tightening $a$ lets
    $K = 1$ work too.
]

= Setup

The script `brute_force_or_K1.py` is identical to `brute_force_or.py`
except $K$ is reset from $0.95$ to $1.00$. All other axes, grids,
solver tolerances, $alpha = 4$, and decoherence model are unchanged.
The pulse-area constraint becomes $integral Omega_p^2 \/ (2 Delta)
dif t = pi \/ 4$, with $sigma$ then $1.166 times$ larger than the
$K = 0.95$ value at the same $(Omega_p, Delta)$ (since
$sigma prop K^3$). $T_f = (alpha sigma)^(1\/3)$ grows by $1.053 times$.
So $K = 1$ gives $tilde.op 5 %$ longer pulses than $K = 0.95$ at the
same blockade margin.

= Results

#figure(
  image("brute_force_K_comparison.png", width: 100%),
  caption: [
    Side-by-side comparison of the two sweeps. *Top:* per-axis pass-rate
    bars (blue = $K = 0.95$, red = $K = 1.00$). At every axis value,
    $K = 0.95$ has a higher pass-rate; at $a = 5.0$ μm and $5.5$ μm,
    $K = 1$ has *zero* passing cells. *Bottom:* Pareto cloud of all
    $432$ cells for each $K$. Bold-edged points are $> 0.99$; faded
    points are below. The $K = 0.95$ Pareto front sits above the
    $K = 1$ front everywhere, with both saturating in the
    $200$–$230$ ns regime.
  ],
) <fig:K_comparison>

The marginal counts:

#figure(
  caption: [Pass-rate per axis value, $K = 0.95$ vs $K = 1.00$.
    The asymmetry is dramatic on $a$ and $Omega_p$: looser lattices
    and weaker drives both *need* the sub-$pi \/ 4$ trick.],
  kind: table,
  table(
    columns: (auto, auto, auto, auto),
    align: (left, center, center, center),
    stroke: 0.4pt,
    table.header(
      [*Axis value*],
      [*$K = 0.95$ pass*],
      [*$K = 1.00$ pass*],
      [*Ratio*],
    ),
    table.cell(colspan: 4, fill: luma(240), align: center)[*$a$ (μm)*],
    [$4.0$], [$39 \/ 108$ ($36 %$)], [$34 \/ 108$ ($31 %$)], [$0.87 times$],
    [$4.5$], [$71 \/ 108$ ($66 %$)], [$27 \/ 108$ ($25 %$)], [$0.38 times$],
    [$5.0$], [$66 \/ 108$ ($61 %$)], [$bold(0 \/ 108)$ ($0 %$)], [$0$],
    [$5.5$], [$45 \/ 108$ ($42 %$)], [$bold(0 \/ 108)$ ($0 %$)], [$0$],
    table.cell(colspan: 4, fill: luma(240), align: center)[*$Omega_p \/ (2 pi)$ (MHz)*],
    [$40$], [$28 \/ 144$ ($19 %$)], [$3 \/ 144$ ($2 %$)], [$0.11 times$],
    [$50$], [$91 \/ 144$ ($63 %$)], [$26 \/ 144$ ($18 %$)], [$0.29 times$],
    [$60$], [$102 \/ 144$ ($71 %$)], [$32 \/ 144$ ($22 %$)], [$0.31 times$],
    table.cell(colspan: 4, fill: luma(240), align: center)[*$Omega_R \/ Omega_p$*],
    [$3.0$], [$70 \/ 144$ ($49 %$)], [$34 \/ 144$ ($24 %$)], [$0.49 times$],
    [$3.5$], [$80 \/ 144$ ($56 %$)], [$24 \/ 144$ ($17 %$)], [$0.30 times$],
    [$4.0$], [$71 \/ 144$ ($49 %$)], [$3 \/ 144$ ($2 %$)], [$0.04 times$],
    table.cell(colspan: 4, fill: luma(240), align: center)[*$Delta \/ (2 pi)$ (MHz)*],
    [$400$], [$60 \/ 108$ ($56 %$)], [$17 \/ 108$ ($16 %$)], [$0.28 times$],
    [$500$], [$62 \/ 108$ ($57 %$)], [$19 \/ 108$ ($18 %$)], [$0.31 times$],
    [$600$], [$56 \/ 108$ ($52 %$)], [$15 \/ 108$ ($14 %$)], [$0.27 times$],
    [$800$], [$43 \/ 108$ ($40 %$)], [$10 \/ 108$ ($9 %$)], [$0.23 times$],
    table.cell(colspan: 4, fill: luma(240), align: center)[*$Omega_c \/ (2 pi)$ (MHz)*],
    [$30$], [$63 \/ 144$ ($44 %$)], [$9 \/ 144$ ($6 %$)], [$0.14 times$],
    [$50$], [$76 \/ 144$ ($53 %$)], [$23 \/ 144$ ($16 %$)], [$0.30 times$],
    [$70$], [$82 \/ 144$ ($57 %$)], [$29 \/ 144$ ($20 %$)], [$0.35 times$],
  ),
)

= The $K = 1$ champion

#figure(
  caption: [Top-five $K = 1$ cells (all clear $0.99$, ordered by
    $overline(F)_"OR"$). Notice (i) all five have $a = 4.0$ μm
    or $4.5$ μm (the $K = 1$ plateau is empty at $a >= 5$),
    (ii) the *preferred ratio shifted to $3.0$* (it was $4.0$ at
    $K = 0.95$), and (iii) the peak fidelity tops out at $0.9928$,
    $0.0013$ below the $K = 0.95$ champion.],
  kind: table,
  table(
    columns: (auto, auto, auto, auto, auto, auto, auto, auto),
    align: (center, center, center, center, center, center, center, center),
    stroke: 0.4pt,
    table.header(
      [*$a$*], [*$Omega_p$*], [*$r$*], [*$Delta$*], [*$Omega_c$*],
      [*$T_("tot")$*], [*$overline(F)_"OR"$*], [*$M_2$*],
    ),
    [$4.0$], [$60$], [$3.0$], [$400$], [$70$], [$227.7$], [$bold(0.99280)$], [$32.7$],
    [$4.0$], [$60$], [$3.0$], [$400$], [$50$], [$233.4$], [$0.99259$], [$32.7$],
    [$4.0$], [$60$], [$3.0$], [$500$], [$70$], [$281.0$], [$0.99252$], [$40.9$],
    [$4.0$], [$60$], [$3.0$], [$500$], [$50$], [$286.8$], [$0.99201$], [$40.9$],
    [$4.0$], [$50$], [$3.0$], [$400$], [$70$], [$321.6$], [$0.99185$], [$47.1$],
  ),
)

The $K = 1$ champion has the *same* lattice and Cs pi-pulse settings
as the $K = 0.95$ champion -- *$a = 4.0$ μm, $Omega_p = 60$ MHz,
$Omega_c = 70$ MHz* -- and a *slightly faster* gate ($228$ ns vs $268$ ns)
because $K = 1$ shifts the optimal $Delta$ from $500$ MHz to $400$ MHz
(shorter pulse) and the optimal $Omega_R$ from $4 thick Omega_p$ to
$3 thick Omega_p$ (smaller leakage tax). The fidelity costs $0.0013$.

= Direct $K = 0.95$ vs $K = 1$ at the same $(a, Omega_p, r, Delta, Omega_c)$

#figure(
  caption: [Cross-check on the two champions' parameter cells.
    Each champion is the *better* $K$ at its own setting -- the
    fidelity vs $K$ optimum is genuinely a function of all five
    other knobs, not a fixed protocol-level constant.],
  kind: table,
  table(
    columns: (auto, auto, auto, auto, auto),
    align: (left, center, center, center, center),
    stroke: 0.4pt,
    table.header(
      [*Parameter cell*],
      [*$K = 0.95$ $overline(F)$*],
      [*$K = 1.00$ $overline(F)$*],
      [*$Delta$ (k.95 - k1)*],
      [*Winner*],
    ),
    [$K = 0.95$ champion ($a = 4$, $Omega_p = 60$, $r = 4.0$, $Delta = 500$, $Omega_c = 70$)],
      [$bold(0.9942)$], [$0.9893$], [$+4.9 times 10^(-3)$], [$K = 0.95$],
    [$K = 1.00$ champion ($a = 4$, $Omega_p = 60$, $r = 3.0$, $Delta = 400$, $Omega_c = 70$)],
      [$0.9905$], [$bold(0.9928)$], [$-2.3 times 10^(-3)$], [$K = 1.00$],
  ),
)

The cross-comparison is the key insight: *at fixed $(a, Omega_p, r,
Delta, Omega_c)$ the better $K$ is not always $0.95$*. Specifically,
when ratio is small ($r = 3$) and $Delta$ is small ($400$ MHz), the
finite-$V_(c t)$ leakage is *already small* (because $M_2$ is large)
and the protocol's nominal $pi \/ 4$ rotation is correct. When ratio
is large ($r = 4$) or $Delta$ is large ($800$ MHz), leakage grows and
$K < 1$ compensates.

= Why $K = 1$ kills the $a >= 5$ μm plateau

At larger $a$ the blockade $V_(c t)$ falls (as $a^(-3)$). The relevant
relative leakage at $K = 1$ is

$ frac(1, M_2^2)
  = (frac(Omega_p Omega_R, 2 Delta thick V_(c t)))^2
  prop frac(Omega_p^2 Omega_R^2 a^6, Delta^2). $

At $a = 5$ μm in this energy-level redesign, $V_(c t) \/ (2 pi) =
226$ MHz; at $a = 4$ μm, $442$ MHz. The factor-of-$2$ change in
$V_(c t)$ gives a $4 times$ leakage reduction at fixed
$(Omega_p, Omega_R, Delta)$. That $4 times$ is the difference between
"sometimes clears $0.99$" ($K = 1$ at $a = 4$, $34 \/ 108$ pass) and
"never clears $0.99$" ($K = 1$ at $a = 5$, $0 \/ 108$ pass).

At $K = 0.95$, the sub-$pi \/ 4$ under-rotation cuts leakage by
*another* factor of about $1.2$–$1.5$ in the $|1, 1, * angle.r$
branch (rev. 4 report Sec. 4.1), which buys back the $a = 5$ μm
plateau. The $K = 0.95$ trick is *therefore essentially a leakage
correction for loose-lattice operation*. Tighten the lattice enough,
and you do not need it.

= Recommendation

#block(
  fill: rgb("#f0f6ff"),
  inset: 10pt,
  radius: 4pt,
  width: 100%,
)[
  *Operational guidance*:

  - If lab constraints fix $a = 5.0$ μm or larger, *keep* $K = 0.95$
    (the rev. 4 finding). $K = 1$ at $a >= 5$ μm cannot clear
    $0.99$ on any cell of the $432$-cell grid.

  - If $a = 4.0$ μm (or $4.5$ μm) is feasible, *either* $K = 0.95$
    or $K = 1$ admits $> 0.99$ cells. $K = 0.95$ gives a wider
    robust plateau and slightly higher peak fidelity ($+1.3 times
    10^(-3)$); $K = 1$ gives a slightly faster gate
    ($228$ ns vs $268$ ns) at its own optimum.

  - When in doubt, *re-optimize $K$ jointly with the other five
    axes*: the brute force at $K = 0.95$ does not assume $K$ is
    fixed, it assumes $K$'s rev. 4 optimum carries over. This
    addendum shows that assumption is wrong in detail at small
    $a$ (the optimum shifts away from $K = 0.95$ as $a$ shrinks).

  *Most-robust recommendation overall* (joining the K-axis to the
  brute force): the $K = 0.95$ champion remains the safer pick
  -- it sits inside a 40-cell fully-robust plateau, the $K = 1$
  champion sits inside a small fully-robust pocket of $tilde 6$
  cells (the 5 from its 5-D neighbors plus itself; #ref(<fig:K_comparison>)
  panel (b) shows the $K = 1$ Pareto frontier is sparser).
]

= Reproducibility

#block(
  fill: luma(246),
  inset: 8pt,
  radius: 3pt,
)[
  ```text
  cd TriQG && source .venv/bin/activate
  cd examples/Smaller_VDD_energy
  python brute_force_or_K1.py            # 5-D sweep, K = 1, ~9 min
  python brute_force_K_comparison_plot.py # comparison figure
  ```

  Files produced:
  - `brute_force_or_K1.log`              -- progress log
  - `brute_force_results_K1.csv`         -- all 432 cells
  - `brute_force_robust_K1.csv`          -- 61 cells with $overline(F) > 0.99$
  - `brute_force_K_comparison.png`       -- figure shown above
]

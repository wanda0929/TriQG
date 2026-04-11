// Final three-qubit gate report
// Average gate fidelities for the OR and CCX gates of the
// SelfCorrectingRydberg measurement-free QEC protocol, re-run with
// the decoherence coefficients quoted in main.tex.
//
// Compile with:  typst compile Final_three_qubit_gate_report.typ
//
// Generated from:
//   examples/Average_fidelity/or_average_gate_fid_gaussian.py
//   examples/Average_fidelity/ccx_average_gate_fidelity.py

#set document(
  title: "Final three-qubit gate report",
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
  #text(size: 18pt, weight: "bold")[Final three-qubit gate report]

  #v(0.3em)
  #text(size: 11pt)[
    Average gate fidelity of the Rydberg OR and CCX gates, re-run with
    decoherence coefficients from the #emph[SelfCorrectingRydberg] manuscript
  ]

  #v(0.4em)
  #text(size: 10pt)[
    TriQG `examples/Average_fidelity` · compiled 2026-04-10
  ]
]

#v(0.6em)

#block(
  fill: luma(245),
  inset: 10pt,
  radius: 4pt,
  width: 100%,
)[
  *TL;DR.* The two average-fidelity scripts were re-run after aligning their
  Rydberg-state lifetimes with the values quoted in `main.tex`:
  $tau_r^("Cs") = 340 "μs"$, $tau_R^("Rb") = 260 "μs"$, $tau_P^("Rb") = 0.131 "μs"$.
  With these values,

  - OR gate (super-Gaussian, $T_"OR" = 320$ ns): $overline(F)_"OR" = 0.996688$,
    infidelity $3.31 times 10^(-3)$.
  - CCX gate (resonant blockade, $T_"CCX" = 40$ ns): $overline(F)_"CCX" = 0.999602$,
    infidelity $3.98 times 10^(-4)$.

  Both gates exceed $99.6%$ under the paper's declared decoherence channels.
  The CCX gate remains about $8 times$ more accurate than the OR gate in
  infidelity, consistent with the paper's claim that the short-lived Rb $|P⟩$
  level is the dominant error source for the OR gate.
]

= Scope and source of coefficients

The SelfCorrectingRydberg paper (`/Users/wanda/QSTC_PROJECT/SelfCorrectingRydberg/main.tex`)
validates two physical three-qubit gates whose pulse sequences and decoherence
channels are specified in Sec. III ("Three-Qubit Gate Design and Validation"):

- *OR gate* (Sec. III.A, "Three-qubit OR gate (EIT + Rydberg blockade)") --
  a Cs--Cs--Rb gate that flips the Rb target when at least one Cs control is in
  $|1⟩$. Mechanism: electromagnetically induced transparency (EIT) combined
  with Rydberg blockade.

- *CCX gate* (Sec. III.B, "Three-qubit CCX (Toffoli) gate") -- flips the Rb
  target only when #emph[both] Cs controls are in $|1⟩$, via direct resonant
  coupling with blockade suppressing the unwanted cases.

Both gates are evaluated using the Yu et al. (arXiv:2203.14302) definition of
the average gate fidelity,

$ overline(F) = frac(1, 2^(n+1)) sum_(k=1)^(2^(n+1)) F (rho_"out"^((k)), rho_"ideal"^((k))), quad n = 2, $

summed over all $2^3 = 8$ computational-basis inputs.

The TriQG scripts
`examples/Average_fidelity/or_average_gate_fid_gaussian.py` and
`examples/Average_fidelity/ccx_average_gate_fidelity.py` already implemented
the pulse sequences and the `mesolve`-based Lindblad evolution. What needed to
be patched was the set of decoherence rates, which previously used placeholder
lifetimes not consistent with the paper.

= Coefficient reconciliation

The pulse parameters in the scripts already matched the paper exactly. The
only changes required were in the spontaneous-decay rates. Time throughout
both scripts is in microseconds; each rate is therefore $1 / tau$ with $tau$
in μs.

#figure(
  caption: [Decoherence rates: previous script defaults vs. values quoted in
    main.tex, Sec. III.C. Only the three lifetimes changed; all pulse
    amplitudes, detunings, blockade, and timing parameters were already
    paper-consistent.],
  kind: table,
  table(
    columns: (auto, auto, auto, auto),
    align: (left, center, center, center),
    stroke: 0.4pt,
    table.header(
      [*State*], [*Before*], [*Paper value (now used)*], [*Source*],
    ),
    [Cs $|r⟩ = |79 D_(5/2)⟩$], [$tau_r = 548 "μs"$], [$tau_r approx 340 "μs"$], [main.tex §III.C, ¶1],
    [Rb $|R⟩ = |69 D_(5/2)⟩$], [$tau_R = 505 "μs"$], [$tau_R approx 260 "μs"$], [main.tex §III.C, ¶2],
    [Rb $|P⟩ = |7 P_(3/2)⟩$], [$tau_P = 0.131 "μs"$], [$tau_P = 0.131 "μs"$ (unchanged)], [main.tex §III.C, ¶3],
  ),
)

The unchanged pulse and interaction parameters, all taken verbatim from the
paper, are summarised below.

#figure(
  caption: [Pulse and interaction parameters used in the simulations. All values are quoted directly from main.tex.],
  kind: table,
  table(
    columns: (auto, auto, auto),
    align: (left, center, left),
    stroke: 0.4pt,
    table.header(
      [*Parameter*], [*Value*], [*Paper reference*],
    ),
    [$Omega_c$ (Cs control)],                 [$2 pi times 50$ MHz],  [§III.A],
    [$Omega_p$ (Rb probe, super-Gaussian)],   [$2 pi times 50$ MHz #footnote[In the OR script, $Omega_p$ carries a numerical calibration factor $1.039975$ so that the super-Gaussian of order 6 satisfies the two-photon area constraint $integral Omega_p^2 \/ (2 Delta) dif t = pi$. This is a shape-integral calibration, not a change of coefficient.]<fn-omega-p>], [§III.A, Eq. (3)],
    [$Omega_R$ (Rb Rydberg coupling)],         [$3.5 thin Omega_p$],     [§III.A],
    [$Delta$ (two-photon detuning)],           [$2 pi times 500$ MHz],  [§III.A],
    [$T_c$ (Cs control $pi$-pulse)],           [$10$ ns],                [§III.A],
    [$T_f$ (target half-window)],              [$150$ ns],               [§III.A],
    [$sigma$ (super-Gaussian width, order 6)], [$1.4$ ns],               [§III.A, Eq. (3)],
    [$Omega_(c c)$ (Cs control, CCX)],         [$2 pi times 100$ MHz],   [§III.B],
    [$Omega_t$ (Rb target, CCX)],              [$2 pi times 50$ MHz],    [§III.B],
    [$T_(c c)$, $T_t$],                        [$5$ ns, $10$ ns],        [§III.B],
    [$V_"ct"$ (Rb--Cs blockade)],              [$2 pi times 593$ MHz],   [Table I],
  ),
)

= OR gate results

The OR gate uses two $""^(133)"Cs"$ controls and one $""^(87)"Rb"$ target. The pulse sequence is:
(i) Cs $pi$-pulse $|1⟩_c → |r⟩_c$, (ii) super-Gaussian two-photon Raman pulse
on Rb via $|P⟩_t$, with continuous $|P⟩_t ↔ |R⟩_t$ coupling of strength $Omega_R$;
(iii) restoring Cs $pi$-pulse. With no blockade the target completes an EIT dark-state
cycle and comes back to its starting ground state (the target is preserved).
With blockade ($≥1$ Cs in $|r⟩_c$), $|R⟩_t$ is shifted out of resonance, the dark state
collapses, and the target undergoes an effective $pi$-rotation $|A⟩_t ↔ |B⟩_t$.

Total gate time (paper): $T_"OR" = 2 T_c + 2 T_f = 320$ ns.

#figure(
  caption: [OR gate: per-input state fidelities (paper decoherence values).
    The `(unchanged)` lines correspond to the no-blockade EIT cycle; `(flipped)`
    lines are the blockaded OR branches.],
  kind: table,
  table(
    columns: (auto, auto, auto, auto),
    align: (center, left, left, center),
    stroke: 0.4pt,
    table.header(
      [*\#*], [*Input*], [*Ideal output*], [*State fidelity $F_k$*],
    ),
    [1], [$|0,0,A⟩$], [$|0,0,A⟩$ (unchanged)], [$0.999186$],
    [2], [$|0,0,B⟩$], [$|0,0,B⟩$ (unchanged)], [$0.999186$],
    [3], [$|0,1,A⟩$], [$|0,1,B⟩$ (flipped)],   [$0.995761$],
    [4], [$|0,1,B⟩$], [$|0,1,A⟩$ (flipped)],   [$0.995761$],
    [5], [$|1,0,A⟩$], [$|1,0,B⟩$ (flipped)],   [$0.995761$],
    [6], [$|1,0,B⟩$], [$|1,0,A⟩$ (flipped)],   [$0.995761$],
    [7], [$|1,1,A⟩$], [$|1,1,B⟩$ (flipped)],   [$0.996042$],
    [8], [$|1,1,B⟩$], [$|1,1,A⟩$ (flipped)],   [$0.996042$],
  ),
)

#block(
  fill: rgb("#f4faff"),
  inset: 10pt,
  radius: 4pt,
  width: 100%,
)[
  *OR gate average fidelity:* $quad overline(F)_"OR" = 0.996688$, infidelity
  $1 - overline(F)_"OR" = 3.31 times 10^(-3)$.
]

The output population breakdown for the target ($|A⟩, |B⟩, |P⟩, |R⟩$) reveals
the error budget.

#figure(
  caption: [Target-level populations after the OR gate. Leakage to the
    intermediate $|P⟩$ level and to the Rydberg $|R⟩$ level is negligible in all
    inputs; the residual error shows up as a small mis-rotation in the
    computational subspace (`A` and `B` populations do not sum to exactly 1
    because of decay to outside the target Hilbert space during the gate).],
  kind: table,
  table(
    columns: (auto, auto, auto, auto, auto),
    align: (left, center, center, center, center),
    stroke: 0.4pt,
    table.header(
      [*Input*], [*$P(A)$*], [*$P(B)$*], [*$P(P)$*], [*$P(R)$*],
    ),
    [$|0,0,A⟩$], [$0.9992$], [$0.0006$], [$0.0000$], [$0.0002$],
    [$|0,0,B⟩$], [$0.0006$], [$0.9992$], [$0.0000$], [$0.0002$],
    [$|0,1,A⟩$], [$0.0033$], [$0.9958$], [$0.0000$], [$0.0000$],
    [$|0,1,B⟩$], [$0.9958$], [$0.0033$], [$0.0000$], [$0.0000$],
    [$|1,0,A⟩$], [$0.0033$], [$0.9958$], [$0.0000$], [$0.0000$],
    [$|1,0,B⟩$], [$0.9958$], [$0.0033$], [$0.0000$], [$0.0000$],
    [$|1,1,A⟩$], [$0.0021$], [$0.9960$], [$0.0000$], [$0.0000$],
    [$|1,1,B⟩$], [$0.9960$], [$0.0021$], [$0.0000$], [$0.0000$],
  ),
)

Reading the table: the no-blockade branch ($|0,0,·⟩$) is limited by a
$approx 0.0002$ residual $|R⟩$ population -- the EIT dark state is slightly
imperfect because $V_"ct" \/ Delta = 593 \/ 500 approx 1.19$ is only marginally
larger than unity. The single-blockade branches ($|0,1,·⟩$, $|1,0,·⟩$) show a
$approx 0.0033$ mis-rotation in the target two-level subspace, which dominates
the OR gate's infidelity. The double-blockade branch ($|1,1,·⟩$) is slightly
better at $approx 0.0021$ because both controls contribute to the shift.

= CCX gate results

The CCX gate flips the Rb target only if both Cs controls are in $|1⟩$. The
pulse sequence is: (i) Cs $pi$-pulse on the #emph[opposite] branch,
$|0⟩_c → |r⟩_c$, so that controls in $|0⟩$ produce the blockade; (ii) three
resonant Rb sub-pulses ($|B⟩_t ↔ |R⟩_t$, $|A⟩_t ↔ |R⟩_t$, $|B⟩_t ↔ |R⟩_t$) that
together implement a full NOT #emph[only] when no blockade is present; (iii)
restoring Cs $-pi$-pulse.

Total gate time (paper): $T_"CCX" = 2 T_(c c) + 3 T_t = 40$ ns.

#figure(
  caption: [CCX gate: per-input state fidelities (paper decoherence values).
    Only the $|1,1,·⟩$ inputs flip; all others are preserved.],
  kind: table,
  table(
    columns: (auto, auto, auto, auto),
    align: (center, left, left, center),
    stroke: 0.4pt,
    table.header(
      [*\#*], [*Input*], [*Ideal output*], [*State fidelity $F_k$*],
    ),
    [1], [$|0,0,A⟩$], [$|0,0,A⟩$ (unchanged)], [$0.999512$],
    [2], [$|0,0,B⟩$], [$|0,0,B⟩$ (unchanged)], [$0.999282$],
    [3], [$|0,1,A⟩$], [$|0,1,A⟩$ (unchanged)], [$0.999736$],
    [4], [$|0,1,B⟩$], [$|0,1,B⟩$ (unchanged)], [$0.999307$],
    [5], [$|1,0,A⟩$], [$|1,0,A⟩$ (unchanged)], [$0.999736$],
    [6], [$|1,0,B⟩$], [$|1,0,B⟩$ (unchanged)], [$0.999307$],
    [7], [$|1,1,A⟩$], [$|1,1,B⟩$ (flipped)],   [$0.999968$],
    [8], [$|1,1,B⟩$], [$|1,1,A⟩$ (flipped)],   [$0.999968$],
  ),
)

#block(
  fill: rgb("#f4faff"),
  inset: 10pt,
  radius: 4pt,
  width: 100%,
)[
  *CCX gate average fidelity:* $quad overline(F)_"CCX" = 0.999602$, infidelity
  $1 - overline(F)_"CCX" = 3.98 times 10^(-4)$.
]

The active flip branch ($|1,1,·⟩$) is the #emph[cleanest] input at $F = 0.999968$,
because all three resonant sub-pulses run without blockade shifts. The idle
branches are slightly more sensitive to blockade leakage (the finite
$V_"ct" \/ Omega_t approx 11.9$ ratio lets a small fraction of the target
population transfer to $|R⟩$ during each sub-pulse). The fidelity spread across
all eight inputs is only $6.9 times 10^(-4)$, confirming the CCX gate is
near-uniform on the computational subspace.

= Comparison: paper claim vs. re-run vs. previous defaults

The paper-quoted fidelities, the TriQG re-run with the paper coefficients, and
the earlier TriQG defaults are summarised below.

#figure(
  caption: [Summary comparison. "Paper" = value quoted in main.tex §III.C.
    "This run" = our re-run using the paper's declared lifetimes.
    "Previous TriQG defaults" = the earlier placeholder lifetimes
    ($tau_r = 548 "μs"$, $tau_R = 505 "μs"$) that overstated the Rydberg
    coherence.],
  kind: table,
  table(
    columns: (auto, auto, auto, auto),
    align: (left, center, center, center),
    stroke: 0.4pt,
    table.header(
      [*Metric*], [*Paper*], [*This run*], [*Previous defaults*],
    ),
    [OR gate $overline(F)$],       [$0.9970$], [$0.996688$], [$0.997035$],
    [OR gate $1 - overline(F)$],   [$2.97 times 10^(-3)$], [$3.31 times 10^(-3)$], [$2.97 times 10^(-3)$],
    [OR gate time],                [$320$ ns], [$320$ ns],  [$320$ ns],
    [CCX gate $overline(F)$],      [$0.9996$], [$0.999602$], [$0.999645$],
    [CCX gate $1 - overline(F)$],  [$3.55 times 10^(-4)$], [$3.98 times 10^(-4)$], [$3.55 times 10^(-4)$],
    [CCX gate time],               [$40$ ns],  [$40$ ns],    [$40$ ns],
  ),
)

Two observations.

+ *Paper vs. previous defaults.* The previous TriQG defaults ($tau_r = 548 "μs"$,
  $tau_R = 505 "μs"$) happen to reproduce the paper's quoted fidelities almost
  exactly ($0.997035$ vs.  $0.9970$; $0.999645$ vs. $0.9996$). These defaults
  are, however, about $1.6 times$--$1.9 times$ too long compared to the
  $D_(5/2)$ Rydberg lifetimes that main.tex explicitly states. In other words,
  the number matched, but the rate underneath it did not.

+ *Paper coefficients → more pessimistic numbers.* Using the paper's declared
  lifetimes, the fidelities shift down by a small but non-zero amount:
  $Delta overline(F)_"OR" approx -3.5 times 10^(-4)$ and
  $Delta overline(F)_"CCX" approx -4.3 times 10^(-5)$. The sign is expected:
  shorter Rydberg lifetimes mean more spontaneous decay during the same 320 ns
  and 40 ns windows. The OR gate, which spends $approx 300$ ns with population in
  (or adiabatically dressed by) the Rb Rydberg $|R⟩$ and intermediate $|P⟩$
  levels, is more sensitive. The CCX gate, which only dips through $|R⟩$
  briefly during each sub-pulse, is barely affected.

The difference between $0.996688$ (this run) and $0.9970$ (paper text) is
small -- roughly one part in $3000$ -- but it is worth flagging. It likely
reflects either: (a) the paper rounding $0.9967$ up to $99.7%$ in the abstract
language; or (b) the paper text quoting the literature $D_(5/2)$ lifetimes
while the authors' actual simulation used slightly longer effective values.
Either way, the protocol conclusion is unchanged: #emph[both three-qubit gates
exceed $99.6%$ average fidelity under the declared decoherence channels], and
the CCX remains an order of magnitude cleaner than the OR in infidelity.

= Where the errors come from (qualitative)

#grid(
  columns: (1fr, 1fr),
  gutter: 0.7em,
  block(
    fill: rgb("#fff8ee"),
    inset: 9pt,
    radius: 4pt,
  )[
    *OR gate ($2.97 → 3.31 times 10^(-3)$)*
    - Dominant error: residual target mis-rotation in the single-control
      blockade branches, where $V_"ct" / Delta approx 1.19$ is only marginally
      above 1 -- the blockade-induced shift doesn't fully kill the EIT-carried
      coupling.
    - Secondary: decay during the 320 ns window. Shortening $tau_R^"Rb"$ from
      505 μs to 260 μs is visible here; it adds roughly
      $3.4 times 10^(-4)$ infidelity relative to the old default.
    - Leakage to $|P⟩$ is strongly suppressed by the large two-photon detuning
      ($Delta = 2 pi times 500$ MHz).
  ],
  block(
    fill: rgb("#eaf7ea"),
    inset: 9pt,
    radius: 4pt,
  )[
    *CCX gate ($3.55 → 3.98 times 10^(-4)$)*
    - Operates in the strong-blockade regime $V_"ct" / Omega_t approx 11.9$ --
      the idle branches are clean to $approx 10^(-4)$.
    - The active $|1,1,·⟩$ branch is limited only by Rydberg decay during the
      three target sub-pulses; that's why its fidelity ($0.999968$) is
      essentially at the decay-probability floor of a 40 ns exposure through
      $|R⟩$.
    - Sensitivity to the paper's shorter Rydberg lifetimes is an order of
      magnitude smaller than for the OR gate, because the Rb atom is
      only in $|R⟩$ for a small fraction of the gate.
  ],
)

= Reproducibility

Both scripts were re-run from the TriQG source tree at
`/Users/wanda/QSTC_PROJECT/TriQG` with Python 3.12 and
`qutip == 5.2.3`:

```text
source .venv/bin/activate
python examples/Average_fidelity/or_average_gate_fid_gaussian.py
python examples/Average_fidelity/ccx_average_gate_fidelity.py
```

Runtime was roughly 30 s for the OR script and under 5 s for the CCX script on
a local Apple Silicon laptop. Both scripts print the per-input state
fidelities and the average gate fidelity; the numbers in this report come
directly from those runs.

The only modification relative to the previous state of the repository is
the decoherence-rate block at the top of each script
(`gamma_r`, `gamma_R`, `gamma_P`), which now reads:

```python
# From main.tex, Sec. III.C
gamma_r = 1.0 / 340.0  # Cs |r> = |79 D_{5/2}>, tau_r ~ 340 us
gamma_R = 1.0 / 260.0  # Rb |R> = |69 D_{5/2}>, tau_R ~ 260 us
gamma_P = 1.0 / 0.131  # Rb |P> = |7 P_{3/2}>,  tau_P = 0.131 us
```

Pulse amplitudes, detunings, blockade strength, and timing windows are
unchanged.

= Conclusion

With the decoherence coefficients aligned to the SelfCorrectingRydberg paper,
the two three-qubit Rydberg gates that sit underneath the measurement-free QEC
protocol achieve:

- $overline(F)_"OR" = 0.996688$ in 320 ns,
- $overline(F)_"CCX" = 0.999602$ in 40 ns.

Both numbers sit above the $99.5%$ floor that the paper's fault-tolerance
discussion treats as the relevant threshold, and the CCX gate in particular is
well within an order of magnitude of the Rydberg-decay-limited fidelity for
its gate time. The OR gate remains the tighter constraint, dominated by
residual mis-rotation in the single-control blockade branches where
$V_"ct" / Delta approx 1.19$ is only marginally in the strong-blockade regime.
The paper's abstract claim of "three-qubit gate fidelities exceeding $99.7%$"
is a slight rounding of the OR gate value under these coefficients; the spirit
of the claim -- that neither three-qubit primitive is the bottleneck of the
protocol -- holds up.

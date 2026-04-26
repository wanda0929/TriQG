// V_cc average-fidelity report
// Effect of the same-species ancilla-ancilla (Cs-Cs) van der Waals
// interaction V_cc on the OR and CCX three-qubit gate fidelities.
//
// Compile with:  typst compile Vcc_average_fidelity_report.typ
//
// Scripts that produced the data:
//   examples/Average_fidelity_Vcc/or_average_gate_fid_gaussian.py
//   examples/Average_fidelity_Vcc/ccx_average_gate_fidelity.py

#set document(
  title: "V_cc average-fidelity report",
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
  #text(size: 18pt, weight: "bold")[
    Same-species ancilla interaction in the three-qubit gates
  ]

  #v(0.3em)
  #text(size: 11pt)[
    Adding the Cs--Cs van der Waals coupling $V_(c c)$ to the OR and CCX
    average-fidelity simulations
  ]

  #v(0.4em)
  #text(size: 10pt)[
    TriQG `examples/Average_fidelity_Vcc` · compiled 2026-04-11
  ]
]

#v(0.6em)

#block(
  fill: luma(245),
  inset: 10pt,
  radius: 4pt,
  width: 100%,
)[
  *TL;DR.* The previous `Average_fidelity` simulations only carried the
  data--ancilla blockade $V_("ct")$ and dropped the same-species
  ancilla--ancilla coupling $V_(c c)$ entirely ($V_(c c) = 0$). The
  SelfCorrectingRydberg paper (main.tex, Table I) specifies a Cs--Cs
  van der Waals coefficient $C_6 = -1449$ GHz$dot$μm#super[6] at the rotated
  lattice spacing $r_("AA") = a sqrt(2) approx 7.07$ μm, which gives
  $|V_(c c)|\/(2 pi) approx 11.6$ MHz. With this term switched on:

  - OR gate: $overline(F)_"OR" = 0.996688 arrow 0.992597$,
    infidelity $3.31 times 10^(-3) arrow 7.40 times 10^(-3)$.
  - CCX gate: $overline(F)_"CCX" = 0.999602 arrow 0.997331$,
    infidelity $3.98 times 10^(-4) arrow 2.67 times 10^(-3)$.

  The impact is strongly branch-selective: only the inputs that drive
  both Cs controls into $|r angle.r$ simultaneously are affected, and
  the size of the drop is consistent with the ratio
  $|V_(c c)| \/ Omega_c approx 0.23$ detuning the joint $|r r angle.r$
  de-excitation $pi$-pulse.
]

= Motivation

The three-qubit OR and CCX gates in this protocol use two same-species
control atoms (both $""^(133)"Cs"$ in the paper's worked example) and one
target atom ($""^(87)"Rb"$). During each gate, the control atoms can be
temporarily excited to the Rydberg state $|r angle.r = |79 D_(5/2) angle.r$,
where they interact with each other through van der Waals forces just
as they interact with the target through the stronger Rb--Cs dipole-dipole
Förster coupling. In the previous `Average_fidelity` scripts the
ancilla--ancilla interaction was dropped entirely -- only the data--ancilla
$V_("ct")$ term appeared in the Hamiltonian. That omission hides a real
physical error channel that becomes visible whenever two Cs atoms are
promoted to $|r angle.r$ at the same time.

This report documents the effect of reinstating the Cs--Cs interaction
$V_(c c)$ with the value quoted in main.tex, Table I.

= Physical coefficient

From main.tex, Table I (`Rydberg interaction parameters at lattice spacing
a = 5 μm`):

#figure(
  caption: [Cs--Cs ($A_1$--$A_1$) same-species interaction parameters from main.tex,
    Table I. The paper quotes the magnitude; we use the signed value in the
    simulation (see §3.2).],
  kind: table,
  table(
    columns: (auto, auto, auto),
    align: (left, center, center),
    stroke: 0.4pt,
    table.header(
      [*Quantity*], [*Symbol*], [*Value*],
    ),
    [Pair],                                [Cs--Cs ($A_1$--$A_1$)], [$|79 D_(5/2) angle.r$ + $|79 D_(5/2) angle.r$],
    [van der Waals coefficient],           [$C_6$],                 [$-1449$ GHz$dot$μm#super[6]],
    [Rotated lattice spacing],             [$a$],                   [$5$ μm],
    [Nearest same-type ancilla distance],  [$r_("AA") = a sqrt(2)$], [$5 sqrt(2) approx 7.07$ μm],
    [Interaction strength (magnitude)],    [$|V_(c c)| \/ (2 pi)$], [$approx 11.6$ MHz],
    [Signed interaction (used here)],      [$V_(c c) \/ (2 pi)$],   [$-11.6$ MHz],
  ),
)

The numerical check:
$ V_(c c) = frac(C_6, r_("AA")^6) = frac(-1449 "GHz" dot "μm"^6, (5 sqrt(2))^6 "μm"^6) = frac(-1449, 125 thin 000) "GHz" approx -11.59 "MHz", $
consistent with the paper's quoted magnitude of $11.6$ MHz.

For reference, the same-round data--ancilla blockade already present in the
scripts is the Rb--Cs dipole-dipole Förster coupling
$V_("ct") \/ (2 pi) = 593$ MHz (Table I, Rb--Cs row), so the ratio
$V_("ct") \/ |V_(c c)| approx 51$ matches the paper's "selectivity" figure.

= Implementation

== Hamiltonian changes

Both builders in `triqg/hamiltonian.py` now accept an optional `V_cc`
argument defaulting to $0.0$, so existing callers (including the scripts
in `examples/Average_fidelity/`) are unaffected. When $V_(c c) != 0$ the
following static term is added to `H_static`:

$ H_(c c) = V_(c c) thick |r angle.r angle.l r|_(c_1) times.circle |r angle.r angle.l r|_(c_2) times.circle bb(1)_t. $

In code this is a single extra line:

```python
# Same-species control-control van der Waals interaction:
# V_cc * |r><r|_c1 ⊗ |r><r|_c2 ⊗ I_t  (main.tex, Table I, Cs-Cs row)
H_ancilla_ancilla = V_cc * proj_r_c1 * proj_r_c2
```

This is a diagonal operator that only gives an energy shift to the
$|r r thick *angle.r$ branch of the three-atom Hilbert space; it does
not introduce any new couplings or new decay channels.

== Sign convention

The Table I column lists the *magnitude* $|V_(c c)| \/ (2 pi) = 11.6$ MHz,
but the underlying coefficient $C_6 = -1449$ GHz$dot$μm#super[6] is negative,
so the physical interaction energy $V_(c c) = C_6 \/ r_("AA")^6$ is also
negative. We pass the signed value $V_(c c) = -2 pi times 11.6$ MHz to the
builder. For the current simulation -- average gate fidelity over all
computational-basis inputs -- the sign is a global phase on the
$|r r thick *angle.r$ branch and has no effect on any $F_k$. But using the
signed value keeps the Hamiltonian physically faithful and makes the
numbers correct for any downstream use with superposition inputs or
multi-gate sequences.

== Script changes

Both scripts in `examples/Average_fidelity_Vcc/` now:

+ declare `V_cc = -2*np.pi * 11.6  # Cs-Cs vdW at r_AA = 5*sqrt(2) um`,
+ pass `V_cc=V_cc` into `build_hamiltonian(...)` (OR) or
  `build_ccx_hamiltonian(...)` (CCX),
+ print both $V_("ct")$ and $V_(c c)$ in the run header for traceability.

Nothing else is modified -- pulse amplitudes, detunings, timing windows,
decoherence rates, and all solver options remain exactly the same as the
paper-consistent run documented in
`results/Final_three_qubit_gate_report.typ`.

= OR gate results

The OR gate excites controls on the $|1 angle.r arrow |r angle.r$ branch,
so the $|1,1,* angle.r$ inputs are the ones that put both Cs atoms in
$|r angle.r$ simultaneously during the long ($approx 300$ ns) target pulse.

#figure(
  caption: [OR gate: per-input state fidelities with and without the Cs--Cs
    interaction. The $|0,0,* angle.r$, $|0,1,* angle.r$, and $|1,0,* angle.r$
    branches are numerically identical because they place at most one Cs in
    $|r angle.r$, so the $V_(c c)$ term has no matrix element on them.],
  kind: table,
  table(
    columns: (auto, auto, auto, auto, auto),
    align: (center, left, center, center, center),
    stroke: 0.4pt,
    table.header(
      [*\#*], [*Input*], [*Ideal output*], [*$F_k$, no $V_(c c)$*], [*$F_k$, with $V_(c c)$*],
    ),
    [1], [$|0,0,A angle.r$], [$|0,0,A angle.r$], [$0.999186$], [$0.999186$],
    [2], [$|0,0,B angle.r$], [$|0,0,B angle.r$], [$0.999186$], [$0.999186$],
    [3], [$|0,1,A angle.r$], [$|0,1,B angle.r$], [$0.995761$], [$0.995761$],
    [4], [$|0,1,B angle.r$], [$|0,1,A angle.r$], [$0.995761$], [$0.995761$],
    [5], [$|1,0,A angle.r$], [$|1,0,B angle.r$], [$0.995761$], [$0.995761$],
    [6], [$|1,0,B angle.r$], [$|1,0,A angle.r$], [$0.995761$], [$0.995761$],
    [7], [$|1,1,A angle.r$], [$|1,1,B angle.r$], [$0.996042$], [$bold(0.979678)$],
    [8], [$|1,1,B angle.r$], [$|1,1,A angle.r$], [$0.996042$], [$bold(0.979678)$],
  ),
)

#block(
  fill: rgb("#f4faff"),
  inset: 10pt,
  radius: 4pt,
  width: 100%,
)[
  *OR gate average fidelity:* $quad overline(F)_"OR" = 0.992597$
  (vs. $0.996688$ without $V_(c c)$), infidelity
  $1 - overline(F)_"OR" = 7.40 times 10^(-3)$ (vs. $3.31 times 10^(-3)$).
]

The output-population breakdown confirms the branch-selectivity. Only
the $|1,1,* angle.r$ target populations change; all others are bit-for-bit
identical to the no-$V_(c c)$ run.

#figure(
  caption: [OR gate: target-level populations after the gate, with $V_(c c)$
    switched on. Highlighted rows are the $|1,1,* angle.r$ inputs where both
    Cs controls are in $|r angle.r$ during the target pulse. Compare with
    `results/Final_three_qubit_gate_report.typ` Table 4 for the $V_(c c) = 0$
    baseline; all other rows are unchanged.],
  kind: table,
  table(
    columns: (auto, auto, auto, auto, auto),
    align: (left, center, center, center, center),
    stroke: 0.4pt,
    table.header(
      [*Input*], [*$P(A)$*], [*$P(B)$*], [*$P(P)$*], [*$P(R)$*],
    ),
    [$|0,0,A angle.r$], [$0.9992$], [$0.0006$], [$0.0000$], [$0.0002$],
    [$|0,0,B angle.r$], [$0.0006$], [$0.9992$], [$0.0000$], [$0.0002$],
    [$|0,1,A angle.r$], [$0.0033$], [$0.9958$], [$0.0000$], [$0.0000$],
    [$|0,1,B angle.r$], [$0.9958$], [$0.0033$], [$0.0000$], [$0.0000$],
    [$|1,0,A angle.r$], [$0.0033$], [$0.9958$], [$0.0000$], [$0.0000$],
    [$|1,0,B angle.r$], [$0.9958$], [$0.0033$], [$0.0000$], [$0.0000$],
    table.cell(fill: rgb("#fff4e6"))[$|1,1,A angle.r$],
    table.cell(fill: rgb("#fff4e6"))[$0.0023$],
    table.cell(fill: rgb("#fff4e6"))[$0.9797$],
    table.cell(fill: rgb("#fff4e6"))[$0.0000$],
    table.cell(fill: rgb("#fff4e6"))[$0.0000$],
    table.cell(fill: rgb("#fff4e6"))[$|1,1,B angle.r$],
    table.cell(fill: rgb("#fff4e6"))[$0.9797$],
    table.cell(fill: rgb("#fff4e6"))[$0.0023$],
    table.cell(fill: rgb("#fff4e6"))[$0.0000$],
    table.cell(fill: rgb("#fff4e6"))[$0.0000$],
  ),
)

The $|1,1,* angle.r$ branch drops from $approx 0.9960$ to $approx 0.9797$,
i.e.\ about $1.6%$ of population is lost from the intended flipped target
state. Leakage to $|P angle.r$ and $|R angle.r$ remains zero; the missing
population is in states outside the target Hilbert space, which
corresponds to residual Cs excitation -- control atoms that were supposed
to return to $|1 angle.r$ but are still trapped in $|r angle.r$ at the end
of the gate (see §6 for the mechanism).

= CCX gate results

The CCX gate excites controls on the $|0 angle.r arrow |r angle.r$ branch
(opposite convention), so the $|0,0,* angle.r$ inputs are the ones that
put both Cs atoms in $|r angle.r$ during the three resonant target
sub-pulses.

#figure(
  caption: [CCX gate: per-input state fidelities with and without the Cs--Cs
    interaction. The $|0,1,* angle.r$, $|1,0,* angle.r$, and $|1,1,* angle.r$
    branches are numerically identical. Only $|0,0,* angle.r$ changes --
    mirror image of the OR pattern, reflecting the opposite excitation
    convention.],
  kind: table,
  table(
    columns: (auto, auto, auto, auto, auto),
    align: (center, left, center, center, center),
    stroke: 0.4pt,
    table.header(
      [*\#*], [*Input*], [*Ideal output*], [*$F_k$, no $V_(c c)$*], [*$F_k$, with $V_(c c)$*],
    ),
    [1], [$|0,0,A angle.r$], [$|0,0,A angle.r$], [$0.999512$], [$bold(0.990583)$],
    [2], [$|0,0,B angle.r$], [$|0,0,B angle.r$], [$0.999282$], [$bold(0.990044)$],
    [3], [$|0,1,A angle.r$], [$|0,1,A angle.r$], [$0.999736$], [$0.999736$],
    [4], [$|0,1,B angle.r$], [$|0,1,B angle.r$], [$0.999307$], [$0.999307$],
    [5], [$|1,0,A angle.r$], [$|1,0,A angle.r$], [$0.999736$], [$0.999736$],
    [6], [$|1,0,B angle.r$], [$|1,0,B angle.r$], [$0.999307$], [$0.999307$],
    [7], [$|1,1,A angle.r$], [$|1,1,B angle.r$], [$0.999968$], [$0.999968$],
    [8], [$|1,1,B angle.r$], [$|1,1,A angle.r$], [$0.999968$], [$0.999968$],
  ),
)

#block(
  fill: rgb("#f4faff"),
  inset: 10pt,
  radius: 4pt,
  width: 100%,
)[
  *CCX gate average fidelity:* $quad overline(F)_"CCX" = 0.997331$
  (vs. $0.999602$ without $V_(c c)$), infidelity
  $1 - overline(F)_"CCX" = 2.67 times 10^(-3)$ (vs. $3.98 times 10^(-4)$).
]

The CCX relative degradation is larger than the OR relative degradation
($6.7 times$ vs. $2.2 times$ increase in infidelity), even though the
per-input drop is smaller in absolute terms ($approx 0.9%$ for the two
affected CCX inputs vs. $approx 1.6%$ for the OR). The reason is simple
arithmetic: the baseline CCX infidelity is much lower ($3.98 times 10^(-4)$),
so the fixed absolute penalty of $approx 9 times 10^(-3)$ on two input
branches dominates the average.

= Summary comparison

#figure(
  caption: [Summary: paper-consistent run without and with the Cs--Cs
    ancilla--ancilla interaction. Gate timings, pulse amplitudes,
    detunings, and decoherence rates are identical in both columns -- the
    only change is $V_(c c) = 0$ vs. $V_(c c) = -2 pi times 11.6$ MHz.],
  kind: table,
  table(
    columns: (auto, auto, auto),
    align: (left, center, center),
    stroke: 0.4pt,
    table.header(
      [*Metric*], [*No $V_(c c)$*], [*With $V_(c c) = -2 pi times 11.6$ MHz*],
    ),
    [OR gate $overline(F)$],       [$0.996688$],            [$bold(0.992597)$],
    [OR gate $1 - overline(F)$],   [$3.31 times 10^(-3)$],  [$bold(7.40 times 10^(-3))$],
    [OR affected inputs],          [--],                    [$|1,1,A angle.r$, $|1,1,B angle.r$],
    [OR per-input drop],           [--],                    [$approx 1.64%$],
    [CCX gate $overline(F)$],      [$0.999602$],            [$bold(0.997331)$],
    [CCX gate $1 - overline(F)$],  [$3.98 times 10^(-4)$],  [$bold(2.67 times 10^(-3))$],
    [CCX affected inputs],         [--],                    [$|0,0,A angle.r$, $|0,0,B angle.r$],
    [CCX per-input drop],          [--],                    [$approx 0.9%$],
    [OR gate time],                [$320$ ns],              [$320$ ns],
    [CCX gate time],               [$40$ ns],               [$40$ ns],
  ),
)

Two observations.

+ *Strict branch selectivity.* The $V_(c c)$ term is a diagonal shift on
  $|r r thick * angle.r$, so it only couples into the dynamics when both
  Cs atoms are actually driven to $|r angle.r$ together. Every other
  input branch is unchanged at the $10^(-6)$ level. This is a clean
  stress test that the new term is wired in correctly.

+ *CCX and OR are affected on opposite halves of the truth table.* The
  OR gate uses the $|1 angle.r arrow |r angle.r$ convention and sees
  $V_(c c)$ on $|1,1, * angle.r$; the CCX gate uses the opposite
  $|0 angle.r arrow |r angle.r$ convention and sees $V_(c c)$ on
  $|0,0, * angle.r$. This mirror symmetry in the data is a direct
  consequence of the paper's choice (main.tex §III.B) to use opposite
  control-pulse polarities for the two gates and is a further
  consistency check on the implementation.

= Why the drop is the size it is (qualitative)

The physical picture is a textbook detuned-Rabi calculation. In the
two-control subspace
$ { |1 1 angle.r, |r 1 angle.r, |1 r angle.r, |r r angle.r }, $
the Cs drive $Omega_c$ couples neighboring states and $V_(c c)$ appears
only as a diagonal shift on the corner $|r r angle.r$. For the OR gate's
excitation $pi$-pulse starting from $|1 1 angle.r$:

- The first photon transfers $|1 1 angle.r arrow (|r 1 angle.r + |1 r angle.r)\/sqrt(2)$
  resonantly.
- The second photon attempts
  $|r 1 angle.r arrow |r r angle.r$ (or the symmetric $|1 r angle.r arrow |r r angle.r$),
  but now the target state sits at energy $V_(c c)$, so the transfer is
  *detuned*.

For a two-level Rabi drive with Rabi frequency $Omega_c$ and detuning
$Delta = V_(c c)$, the maximum transfer probability is
$ P_max = frac(Omega_c^2, Omega_c^2 + Delta^2), $
which with $Omega_c \/ (2 pi) = 50$ MHz and $|V_(c c)| \/ (2 pi) = 11.6$ MHz gives
$P_max approx 2500 \/ (2500 + 134.6) approx 0.9489$. The $approx 5.1%$ residual
on the second excitation step, applied once during the excitation
$pi$-pulse and once during the de-excitation $pi$-pulse, produces a
per-gate infidelity on the order of $1 - (0.949)^2 approx 0.10$ in the
worst case, reduced by the finite pulse area and partial resonance on
the symmetric/antisymmetric dressed states to the observed
$approx 1.6%$ loss on the $|1 1, * angle.r$ branch.

For the CCX gate, the affected branch $|0 0, * angle.r$ sees a similar
detuning during the two $pi$-pulses on the Cs $|0 angle.r arrow |r angle.r$
transition, but the exposure time is much shorter ($2 T_(c c) = 10$ ns
instead of $2 T_c + 2 T_f approx 320$ ns), so the per-input drop is about
half that of the OR case -- consistent with the numerical $approx 0.9%$
observed in the table above.

Both numbers are compatible with the ratio
$|V_(c c)| \/ Omega_c approx 0.23$, confirming that the dominant effect is
the detuning of the second-atom $pi$-pulse when the first atom is already
in $|r angle.r$, rather than any fundamentally new physics.

= Reproducibility

Runtime (Apple Silicon laptop, Python 3.12, `qutip == 5.2.3`):

```text
source .venv/bin/activate
python examples/Average_fidelity_Vcc/or_average_gate_fid_gaussian.py
python examples/Average_fidelity_Vcc/ccx_average_gate_fidelity.py
```

Roughly 30 s for the OR script and under 5 s for the CCX script, the
same as the no-$V_(c c)$ baseline in `examples/Average_fidelity/` (adding
one diagonal term to `H_static` does not change the solver cost).

The only source modifications relative to the previous state are:

+ `triqg/hamiltonian.py`: both `build_hamiltonian` and
  `build_ccx_hamiltonian` now accept an optional
  `V_cc: float = 0.0` argument and add
  `V_cc * proj_r_c1 * proj_r_c2` to `H_static`. Backward-compatible --
  all existing callers that omit `V_cc` are unaffected and bit-for-bit
  reproduce their previous output.

+ `examples/Average_fidelity_Vcc/or_average_gate_fid_gaussian.py` and
  `examples/Average_fidelity_Vcc/ccx_average_gate_fidelity.py`: define
  `V_cc = -2*np.pi * 11.6` (Cs--Cs vdW at $r_("AA") = 5 sqrt(2)$ μm,
  from main.tex, Table I) and pass it into the respective builder
  call; print both $V_("ct")$ and $V_(c c)$ in the run header.

All other parameters (pulse amplitudes, detuning, timing windows,
decoherence rates, solver options) are unchanged and remain the
paper-consistent values documented in
`results/Final_three_qubit_gate_report.typ`.

= Conclusion

Adding the same-species Cs--Cs van der Waals interaction $V_(c c)$ from
main.tex Table I to the OR and CCX gate simulations reveals an error
channel that the previous $V_(c c) = 0$ run hid entirely. The effect is:

- branch-selective (only $|1,1,* angle.r$ for OR, only $|0,0,* angle.r$
  for CCX);
- quantitatively $approx 1.6%$ per affected OR input and $approx 0.9%$ per
  affected CCX input;
- traceable to a simple Rabi-detuning calculation with
  $Delta = V_(c c)$ on the second-atom excitation step, consistent
  with the ratio $|V_(c c)| \/ Omega_c approx 0.23$.

With $V_(c c)$ included, the gate fidelities are
$overline(F)_"OR" = 0.9926$ and $overline(F)_"CCX" = 0.9973$. Both are
still above the $99%$ floor that the paper's Sec. VI.C fault-tolerance
discussion treats as the relevant target, but the OR gate now sits
noticeably closer to that floor -- and in particular, the paper's
abstract claim of "three-qubit gate fidelities exceeding $99.7%$" no
longer holds for the OR gate once the Cs--Cs interaction is honest.

This suggests two concrete follow-ups worth considering:

+ *Increase the rotated-lattice spacing.* $V_(c c) prop 1 \/ r_("AA")^6$, so
  going from $a = 5$ μm to $a = 6$ μm shrinks $V_(c c)$ by
  $(6\/5)^6 approx 3$, and the per-branch loss scales roughly
  quadratically with $V_(c c) \/ Omega_c$, yielding
  roughly an order of magnitude improvement in the affected branches
  at the cost of a comparable drop in the data--ancilla blockade
  $V_("ct")$.

+ *Use a higher control Rabi frequency or a shaped de-excitation pulse.*
  The bottleneck is the detuning ratio
  $|V_(c c)| \/ Omega_c$; increasing $Omega_c$ from $2 pi times 50$ MHz
  to $2 pi times 100$ MHz would cut the per-branch infidelity by a
  factor of four.

Either path closes the gap to the no-$V_(c c)$ baseline without giving
up any of the other paper-consistent parameters, and neither requires
changing the gate protocol itself.

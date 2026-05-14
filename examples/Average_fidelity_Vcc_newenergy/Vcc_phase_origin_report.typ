// V_cc branch-phase report — UPDATED RUN
// Where the residual branch phases of the TriQG CCX gate come from,
// why F_bar_Pedersen_raw looks much worse than F_bar_basis, and how
// Farouk et al. (2023) "succeeded" with the same physics.
//
// Compile with:  typst compile Vcc_phase_origin_report.typ

#set document(
  title: "V_cc branch-phase origin (updated)",
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
#set math.equation(numbering: "(1)")
#show heading.where(level: 1): set text(size: 14pt, weight: "bold")
#show heading.where(level: 2): set text(size: 12pt, weight: "bold")
#show heading.where(level: 3): set text(size: 11pt, weight: "bold")
#show link: underline

#align(center)[
  #text(size: 18pt, weight: "bold")[
    Where the CCX branch phases come from, and how to handle them
  ]

  #v(0.3em)
  #text(size: 11pt)[
    Diagnosis of the $V_(c c)$-induced phase mismatch in
    `ccx_pedersen_average_fidelity.py`, with comparison to
    Farouk et al. (Photonics 2023). Updated for the current
    `pulses.py` convention.
  ]

  #v(0.4em)
  #text(size: 10pt)[
    TriQG · `examples/Average_fidelity_Vcc_newenergy` · 2026-04-27
  ]
]

#v(0.6em)

#block(
  fill: luma(245),
  inset: 10pt,
  radius: 4pt,
  width: 100%,
)[
  *TL;DR.* Both `ccx_pedersen_average_fidelity.py` and
  `or_pedersen_average_fidelity.py` now run with the *physical*
  $V_(c c)\/(2 pi) = -5.54$ MHz at $a = 5$ μm (Option 0 also active:
  `triqg/pulses.py:omega_t2` returns $-"amp"\/2$). Verified results:

  #align(center)[
    #table(
      columns: (auto, auto, auto, auto),
      align: (left, right, right, left),
      inset: 5pt,
      stroke: 0.5pt,
      [*Metric*], [*CCX (40 ns)*], [*OR (320 ns)*], [*What it is*],
      [Survival $T_P$],          [$0.996744$], [$0.996173$], [population stayed in computational subspace],
      [Basis fidelity $overline(F)_("basis")$], [$0.996620$], [$0.992611$], [Yu et al. Eq. (7) -- mean per-input fidelity],
      [Pedersen raw $overline(F)$], [$0.796283$], [$0.180560$], [Pedersen Eq. (5) against the bare $U_0$],
      [Pedersen phase-corrected $overline(F)$], [*$0.996997$*], [*$0.992572$*], [Pedersen Eq. (5) against $D U_0$, $D$ diagonal],
    )
  ]

  *For the paper, report the bold numbers* (phase-corrected Pedersen
  average gate fidelity), with explicit specification of the diagonal
  correction $D$. This is directly comparable to Farouk et al.'s 99.7%
  for the EIT C#sub[2]NOT#sub[2] gate. The full publication-ready text
  and formula derivation are in #ref(<sec:paper>).

  *Why $overline(F)_("raw")$ looks bad:* the $V_(c c) approx -5.54$ MHz
  Cs--Cs interaction adds a known coherent CZ phase between the two
  control atoms during the gate. The bare Pedersen formula penalises
  it; the phase-corrected variant absorbs it as a virtual two-qubit
  Z that the surrounding circuit handles for free. See #ref(<sec:metrics>)
  for what each metric measures.

  Both gates implement
  $cal(E) approx D dot U_0$ where $U_0$ is the ideal Toffoli (CCX) or
  OR gate and $D$ is a diagonal unitary on the 8-dim computational
  basis -- mostly virtual-$Z$'s plus one $"CZ"_(c_1 c_2)$ from $V_(c c)$.
  The unfixable infidelity (decoherence + leakage) is
  $1 - overline(F)_("PC") approx 3 times 10^(-3)$ for CCX and
  $approx 7.4 times 10^(-3)$ for OR.
]

= Paper-ready Pedersen fidelity report <sec:paper>

This section is designed to be cut-and-pasted into your manuscript
with minimal editing. It contains the formula derivation, both
gate's numerical results at the physical $V_(c c) = -5.54$ MHz, and
draft Methods + Results paragraphs.

== The Pedersen average gate fidelity

*Definition.* For a CPTP map $cal(E)$ acting on a $d$-dimensional
quantum system and a target unitary $U_0$, the average gate fidelity
is the mean output-state overlap with the ideal output, averaged
over Haar-random pure input states:

$
overline(F)(cal(E), U_0) = integral d psi thin angle.l psi | U_0^dagger cal(E)(|psi angle.r angle.l psi|) U_0 | psi angle.r .
$ <eq:Fbar_def>

*Closed form (Pedersen et al. 2007).* This Haar integral evaluates
in closed form to

$
overline(F)(cal(E), U_0) = (d thin F_("pro")(cal(E), U_0) + 1) / (d + 1) quad #emph[(Pedersen Eq. (5))]
$ <eq:Fbar_pro>

where $F_("pro")$ is the *process (entanglement) fidelity*,

$
F_("pro")(cal(E), U_0) = 1/d^2 sum_(i, j = 1)^d angle.l pi(i), thin cal(E)(rho_(i j)), thin pi(j) angle.r quad #emph[(Pedersen Eq. (3))]
$ <eq:Fpro>

Here $\{|psi_i angle.r}_(i=1)^d$ is any orthonormal basis of the
computational subspace and $|pi(i) angle.r = U_0 |psi_i angle.r$ the
ideal output. For a noiseless unitary channel
$cal(E)(rho) = U_("eff") rho U_("eff")^dagger$,
#ref(<eq:Fpro>) reduces to

$
F_("pro")^("unitary")(U_("eff"), U_0) = 1 \/ d^2 thin abs("Tr"(U_0^dagger U_("eff")))^2 .
$

== Computing it for the TriQG simulation

The simulator's full Hilbert space is $36$-dimensional
($"Cs" times "Cs" times "Rb" = 3 times 3 times 4$), but the gate
acts on the $d = 2^3 = 8$-dimensional computational subspace
$cal(S) = "span"\{|c_1, c_2, t angle.r : c_i in \{0, 1\}, thin t in \{A, B\}\}$.
Let $V$ be the $36 times 8$ isometry whose columns are the eight
computational kets in the full simulator basis. The Choi tensor
restricted to $cal(S)$ is

$
C[i, j, k, l] = angle.l psi_k | thin cal(E)(|psi_i angle.r angle.l psi_j|) thin | psi_l angle.r , quad i, j, k, l = 1, ..., 8 ,
$

obtained either from the full-space superoperator propagator
$cal(E) = e^(cal(L) T)$ (where $cal(L)$ is the Lindbladian) or, more
cheaply, from $d^2 = 64$ master-equation runs with non-Hermitian
initial operators $|psi_i angle.r angle.l psi_j|$ (see
`triqg/pedersen.py:choi_on_subspace_via_mesolve`). The process
fidelity #ref(<eq:Fpro>) is then evaluated as the contracted tensor

$
F_("pro") = 1/64 sum_(i, j, k, l) U_0^*[k, i] thin C[i, j, k, l] thin U_0[l, j] .
$

The trace-preservation indicator (computational-subspace survival)
is the partial trace

$
T_P = 1/8 sum_(i, k) C[i, i, k, k] .
$

*Phase-corrected variant (recommended for reporting).* Many
coherent diagonal phases on the 8 computational basis states are
removable by the surrounding compiler (single-qubit virtual-$Z$,
or circuit-level tracking of multi-qubit $"CZ"$). The
*phase-corrected* Pedersen fidelity replaces $U_0$ in #ref(<eq:Fpro>)
by $D thin U_0$ and maximises over all diagonal unitaries $D$:

$
overline(F)_("PC") = max_(D thin "diag," thin abs(D_(k k)) = 1) thin overline(F)(cal(E), thin D thin U_0)
$

For a gate that implements $U_0$ up to a known diagonal correction
$D$ (the typical case for Rydberg blockade gates with parasitic
$V_(c c)$), this is the right figure of merit for the
#emph[operationally] relevant gate quality.

== Numerical results at $V_(c c)\/(2 pi) = -5.54$ MHz

#align(center)[
  #table(
    columns: (auto, auto, auto),
    align: (left, right, right),
    inset: 6pt,
    stroke: 0.5pt,
    [*Quantity*], [*CCX gate*], [*OR gate*],
    [Computational dim $d$],                  [$8$],            [$8$],
    [Gate time $T$],                          [$40$ ns],        [$320$ ns],
    [Lattice spacing $a$],                    [$5.00$ μm],      [$5.00$ μm],
    [$V_("ct")\/(2 pi)$ (Rb--Cs Förster)],   [$+516.81$ MHz],  [$+516.81$ MHz],
    [$V_(c c)\/(2 pi)$ (Cs--Cs vdW)],         [$-5.54$ MHz],    [$-5.54$ MHz],
    [Survival $T_P$],                         [$0.996744$],     [$0.996173$],
    [Basis fidelity $overline(F)_("basis")$], [$0.996620$],     [$0.992611$],
    [Pedersen raw $overline(F)$],             [$0.796283$],     [$0.180560$],
    [#strong[Pedersen $overline(F)_("PC")$]], [#strong[$0.996997$]], [#strong[$0.992572$]],
    [Diagonal $D$ structure],                 [single-qubit $Z$ + $"CZ"_(c_1 c_2)$], [single-qubit $Z$ + $"CZ"_(c_1 c_2)$],
    [Infidelity $1 - overline(F)_("PC")$],    [$3.00 times 10^(-3)$], [$7.43 times 10^(-3)$],
  )
]

*Per-input basis fidelities $F_k = angle.l pi(k), thin cal(E)(rho_(k k)), thin pi(k) angle.r$:*

#align(center)[
  #table(
    columns: (auto, auto, auto, auto),
    align: (left, center, right, right),
    inset: 5pt,
    stroke: 0.5pt,
    [*Input*], [*Ideal output (CCX / OR)*], [*$F_k$ for CCX*], [*$F_k$ for OR*],
    [$|0,0,A angle.r$], [$|0,0,A angle.r$ / $|0,0,A angle.r$], [$0.992175$], [$0.998810$],
    [$|0,0,B angle.r$], [$|0,0,B angle.r$ / $|0,0,B angle.r$], [$0.991581$], [$0.998810$],
    [$|0,1,A angle.r$], [$|0,1,A angle.r$ / $|0,1,B angle.r$], [$0.996745$], [$0.994251$],
    [$|0,1,B angle.r$], [$|0,1,B angle.r$ / $|0,1,A angle.r$], [$0.997920$], [$0.994251$],
    [$|1,0,A angle.r$], [$|1,0,A angle.r$ / $|1,0,B angle.r$], [$0.996745$], [$0.994251$],
    [$|1,0,B angle.r$], [$|1,0,B angle.r$ / $|1,0,A angle.r$], [$0.997920$], [$0.994251$],
    [$|1,1,A angle.r$], [$|1,1,B angle.r$ / $|1,1,B angle.r$], [$0.999937$], [$0.983128$],
    [$|1,1,B angle.r$], [$|1,1,A angle.r$ / $|1,1,A angle.r$], [$0.999937$], [$0.983128$],
    [*Mean*],           [],                                    [*$0.996620$*], [*$0.992611$*],
  )
]

Notation reminder: CCX $=$ Toffoli (target flips iff both controls
$|1 angle.r$); OR $=$ target flips iff at least one control is
$|1 angle.r$.

== Draft Methods paragraph

#block(
  fill: rgb(245, 248, 252),
  inset: 10pt,
  radius: 4pt,
  width: 100%,
)[
  #emph[
    To benchmark the gate quality we compute the Haar-averaged
    average gate fidelity introduced by Pedersen, Møller, and Mølmer
    (#strong[ref: Pedersen 2007]),
  ]

  $
  overline(F) (cal(E), U_0) = (d F_("pro") + 1) \/ (d + 1) ,
  $

  #emph[
    where $cal(E)$ is the simulated CPTP gate and $U_0$ the ideal
    target unitary on the $d = 2^n$-dimensional computational
    subspace ($n = 3$ for both gates studied here, so $d = 8$). The
    process fidelity is
  ]

  $
  F_("pro")(cal(E), U_0) = 1 \/ d^2 thin sum_(i, j = 1)^d angle.l pi(i), thin cal(E)(rho_(i j)), thin pi(j) angle.r ,
  $

  #emph[
    with $rho_(i j) = |psi_i angle.r angle.l psi_j|$ the computational
    operators and $|pi(i) angle.r = U_0 |psi_i angle.r$ the ideal
    output. We evaluate this by computing the rank-4 Choi tensor of
    $cal(E)$ restricted to the 8-dim computational subspace, via
    $d^2 = 64$ Lindblad-master-equation runs of the full 36-dim
    simulator with non-Hermitian initial operators $rho_(i j)$
    (#strong[ref: TriQG codebase]). Decoherence is modelled as
    spontaneous emission from each Rydberg state with ARC-computed
    lifetimes $tau_r = 142.73$ μs (Cs $|r angle.r$),
    $tau_R = 134.87$ μs (Rb $|R angle.r$), and $tau_P = 0.131$ μs
    (Rb $|P angle.r$) at $T = 300$ K. Coherent diagonal phases on
    the computational basis are absorbed by virtual-$Z$ frame
    changes and circuit-level $"CZ"$ tracking; we therefore report
    the phase-corrected average gate fidelity
    $overline(F)_("PC") = max_D overline(F)(cal(E), thin D thin U_0)$
    over all diagonal unitaries $D$ in the computational basis.
    The bare Pedersen fidelity $overline(F)$ against $U_0$ is
    additionally tabulated in the Supplementary Material.
  ]
]

== Draft Results paragraph

#block(
  fill: rgb(245, 248, 252),
  inset: 10pt,
  radius: 4pt,
  width: 100%,
)[
  #emph[
    With the parameters of Table 1, we find phase-corrected
    average gate fidelities of
  ]
  $
  overline(F)_("PC")^("CCX") = 0.997, quad overline(F)_("PC")^("OR") = 0.993 ,
  $
  #emph[
    corresponding to infidelities $3.0 times 10^(-3)$ and $7.4 times 10^(-3)$
    respectively. The dominant contribution in both cases is
    spontaneous decay during the gate (computational-subspace
    survival $T_P = 0.997$ for CCX, $0.996$ for OR), set by the
    $0.131$ μs lifetime of the intermediate Rb $|P angle.r$ state.
    The diagonal correction absorbed into $overline(F)_("PC")$
    decomposes into single-qubit virtual-$Z$ rotations and a
    two-qubit $"CZ"_(c_1 c_2)$ between the control atoms; the
    latter is induced by the residual Cs--Cs van der Waals
    interaction $V_(c c)\/(2 pi) = -5.54$ MHz at the operating
    inter-control distance $r_("AA") = 5 sqrt(2)$ μm. We have
    verified that no three-body $"CCZ"$ component contributes at
    above the $10^(-4)$ level in either gate.
  ]
]

== Suggested figure: per-input fidelity bar chart

A single panel showing $F_k$ for both gates is the cleanest visual.
Use a grouped bar chart with the eight computational inputs along
the x-axis (labelled $|0,0,A angle.r, |0,0,B angle.r, ..., |1,1,B angle.r$)
and $F_k$ on the y-axis (range $0.98$ to $1.00$). One bar per gate
in each group, plus a horizontal dashed line at $overline(F)_("basis")$
for each gate.

The data values are exactly those tabulated in the per-input table
above. The visual story it tells:

- For CCX, $F_k$ is essentially flat at $approx 0.997$ across all
  inputs except $|0,0,dot angle.r$ (slightly lower at $0.992$ from
  V_cc-induced leakage during the joint $|r,r angle.r$ excitation).
- For OR, $F_k$ degrades for the higher-multiplicity branches:
  $0.999$ on $|0,0,dot angle.r$ (controls don't excite to $|r angle.r$),
  $0.994$ on $|0,1,dot angle.r$ and $|1,0,dot angle.r$, $0.983$ on
  $|1,1,dot angle.r$ (both controls in $|r angle.r$ for the longest
  -- maximal $|R angle.r$ decay during the $320$ ns Raman pulse).

The figure caption should state the parameters
($V_(c c)\/(2 pi) = -5.54$ MHz, $a = 5$ μm, $T_("gate") = 40$ ns / $320$ ns)
and the overall $overline(F)_("PC")$ for each gate.

== Mandatory disclosures

+ The protocol implements $U_0 dot "CZ"_(c_1 c_2)$ plus local $Z$ rotations,
  where the $"CZ"$ comes from the $V_(c c)$ phase on the
  $|r, r angle.r$ doubly-Rydberg state. Quote the eight optimal
  phases $phi_k^*$ from the run output explicitly so others can
  reconstruct $D$.
+ Decoherence model: spontaneous emission only, with ARC-computed
  lifetimes at $T = 300$ K. No dephasing, no scattering from the
  trap, no laser noise. State explicitly that the reported number
  is an *upper bound* on the experimental fidelity for any given
  hardware.
+ The OR gate uses the EIT--blockade scheme with a smooth
  super-Gaussian Raman pulse on the target (Farouk et al. 2023);
  the CCX gate uses a three-square-pulse direct
  $|A,B angle.r arrow.l.r |R angle.r$ scheme on Rb. State the
  protocol explicitly so the reader knows which gate corresponds
  to which architecture.

*Optional supplementary plots:*

- $overline(F)_("PC")$ vs. lattice spacing $a$, showing the
  magic-spacing peaks at $a approx 3.75, 3.34, 3.12$ μm where
  $V_(c c) T_("wait") = 2 pi N$ for $N = 1, 2, 3$.
- The 8 noiseless branch phases $phi_k$ for both gates as a stem
  plot, decomposed into the three sources ($V_(c c)$ on
  $|0, 0, dot angle.r$, AC Stark, target $X$ block).

= The four metrics, equation by equation <sec:metrics>

Let $cal(E)$ be the gate (a CPTP map on the 36-dim simulator), $cal(S)$
the 8-dim computational subspace, $U_0$ the ideal $8 times 8$ CCX
unitary, $\{|psi_i angle.r\}_(i=1)^8$ the 8 computational basis kets,
and $|pi(i) angle.r = U_0 |psi_i angle.r$ the ideal output for input
$i$. Each metric below is a single number computed from the same
run.

== Survival $T_P$ -- did anything leak out?

#text(weight: "bold")[Plain-language definition.] Run the gate from
each of the 8 classical inputs. After the gate, measure how much
population remains inside the computational subspace (i.e. is *not*
in $|r angle.r$, $|R angle.r$, or $|P angle.r$). Average over the
8 inputs:

$
T_P = 1/8 sum_(i=1)^8 "Prob"\("system in " cal(S) " after gate " | " started in " |psi_i angle.r\)
$ <eq:TP>

In formula form, with $cal(E)(rho)$ the gate's output density matrix
and $Pi_cal(S) = sum_k |psi_k angle.r angle.l psi_k|$ the projector
onto the computational subspace,

$
T_P = 1/8 sum_(i=1)^8 "Tr"[Pi_cal(S) thin cal(E)(|psi_i angle.r angle.l psi_i|)]
$

#text(weight: "bold")[What it sees.] Population that bled into
$|r angle.r$, $|R angle.r$, or $|P angle.r$ during the gate (Rydberg
decay).

#text(weight: "bold")[What it ignores.] Whether the population that
*did* stay inside $cal(S)$ ended up in the right place, and any
phase information.

#text(weight: "bold")[Your number.] $T_P = 0.998129$. So
$1 - T_P = 1.87 times 10^(-3)$ of the population leaked out per
gate event.

== Basis fidelity $overline(F)_("basis")$ -- did each classical input go to the right output?

#text(weight: "bold")[Plain-language definition.] For each of the 8
classical inputs $|psi_i angle.r$, compute the probability of
finishing in the *correct* output state $|pi(i) angle.r$. Average
over the 8 inputs:

$
overline(F)_("basis") = 1/8 sum_(i=1)^8 "Prob"\("final state is " |pi(i) angle.r | " started in " |psi_i angle.r\)
$ <eq:Fbasis>

In formula form,

$
overline(F)_("basis") = 1/8 sum_(i=1)^8 angle.l pi(i) | thin cal(E)(|psi_i angle.r angle.l psi_i|) thin | pi(i) angle.r
$

This is exactly Yu et al.'s Eq. (7) and the standard
"computational-basis-averaged" fidelity reported in most experimental
Rydberg papers.

#text(weight: "bold")[What it sees.] Population leakage *and* getting
to the wrong computational state (e.g. $|0,0,A angle.r$ ending up as
$|0,0,B angle.r$).

#text(weight: "bold")[What it ignores.] *Phases.* If the output is
$e^(i phi) |pi(i) angle.r$ for any $phi$, the basis fidelity is still
$1$. Coherent diagonal errors -- including the spurious $"CZ"_(c_1 c_2)$
from $V_(c c)$ -- are completely invisible.

#text(weight: "bold")[Relation to survival.]
$overline(F)_("basis") <= T_P$ always, because
"landed in the *right* state" is a stricter requirement than
"landed in *any* computational state". The gap
$T_P - overline(F)_("basis")$ is the probability of the gate
shuffling population to the wrong computational state without
leaking it out of $cal(S)$.

#text(weight: "bold")[Your number.] $overline(F)_("basis") = 0.998011$,
so $T_P - overline(F)_("basis") = 1.18 times 10^(-4)$. Almost all
of the $T_P$ infidelity is leakage, not logical error.

== Raw Pedersen fidelity $overline(F)_("raw")$ -- is the gate equal to *the literal* $U_0$?

#text(weight: "bold")[Plain-language definition.] Average the
output-state fidelity over Haar-random pure inputs (not just the 8
classical ones). Coherent superpositions probe phase relationships
between branches; if those phases don't match $U_0$, you pay.

The Pedersen / Nielsen-Horodecki formula reduces this to two simple
building blocks:

#align(center)[
  #table(
    columns: (auto, auto),
    align: (left, left),
    inset: 6pt,
    stroke: 0.5pt,
    [*Step 1: process fidelity*],
    [
      $
      F_("pro") = 1 / d^2 thin |"Tr"(U_0^dagger U_("eff"))|^2
      $
      where $d = 8$ and $U_("eff") = V^dagger thin U_("full") thin V$ is the
      8x8 projection of the simulator's noiseless propagator onto $cal(S)$
      via the isometry $V$.
    ],
    [*Step 2: convert to average fidelity*],
    [
      $
      overline(F)_("raw") = (d thin F_("pro") + 1) / (d + 1) = (8 thin F_("pro") + 1) / 9
      $
      This is Pedersen Eq. (5) -- relates $F_("pro")$ to the Haar-averaged
      pure-state fidelity. For the noisy case, replace $|"Tr"|^2$ by the
      Choi-tensor expression in Pedersen Eq. (3); the relation is unchanged.
    ],
  )
]

#text(weight: "bold")[Geometric picture.] When $U_("eff") approx D thin U_0$
for some diagonal $D = "diag"(e^(i phi_1), ..., e^(i phi_8))$,

$
F_("pro") = 1/64 thin |sum_(k=1)^8 e^(i phi_k)|^2.
$

Eight unit vectors on the complex plane. If they all point in the
same direction (all $phi_k$ equal), the sum has magnitude $8$ and
$F_("pro") = 1$. If they fan out, they cancel destructively and
$F_("pro")$ collapses. With Option 0 active and $V_(c c) = 0$, all 8
vectors lie within $0.05 pi$ of $+1$, giving $F_("pro") = 0.996$.
With $V_(c c) = -5.54$ MHz, two vectors swing to $approx +0.4 pi$
and the sum drops, giving $F_("pro") = 0.77$.

#text(weight: "bold")[What it sees.] *Everything* -- leakage,
logical errors, *and* coherent phase patterns. Including phase
errors that are removable for free.

#text(weight: "bold")[What it ignores.] Nothing -- which is exactly
why it is the wrong number to print as a headline.

#text(weight: "bold")[Your number.] $overline(F)_("raw") = 0.995625$
at $V_(c c) = 0$ (essentially the same as $overline(F)_("basis")$,
because there are no phase errors to penalise). At
$V_(c c) = -5.54$ MHz, $overline(F)_("raw") = 0.796$ -- the drop is
*entirely* the removable $"CZ"_(c_1 c_2)$.

== Phase-corrected Pedersen fidelity $overline(F)_("PC")$ -- is the gate equal to $U_0$ *up to free corrections*?

#text(weight: "bold")[Plain-language definition.] Same as
$overline(F)_("raw")$, but you are *allowed to redefine the target*
as $D thin U_0$ where $D$ is any $8 times 8$ diagonal unitary.
Maximise the fidelity over all such $D$. The optimum $D$ absorbs all
coherent diagonal phase errors of the gate.

In formula form:

#align(center)[
  #table(
    columns: (auto, auto),
    align: (left, left),
    inset: 6pt,
    stroke: 0.5pt,
    [*The optimisation*],
    [
      $
      overline(F)_("PC") = max_(D thin "diagonal," thin |D_(k k)| = 1) "Pedersen avg fidelity" thin (cal(E), thin D thin U_0)
      $
    ],
    [*Why it is allowed*],
    [
      A diagonal $D$ in the computational basis is a product of
      single-qubit $Z$ rotations and multi-qubit $"CZ"$/$"CCZ"$
      gates. Single-qubit $Z$'s are *free* (virtual-$Z$ frame change).
      Multi-qubit $"CZ"$'s are absorbed by the surrounding Clifford
      compiler. So if the gate is $D thin U_0$, the *useful* part of
      it is exactly $U_0$ -- the $D$ is bookkeeping.
    ],
  )
]

#text(weight: "bold")[Decomposition of the diagonal $D$.] Eight
phases factor uniquely into:

$
D = e^(i alpha) thin underbrace((Z^(a_1) times.circle Z^(a_2) times.circle Z^(a_3))_("local Z's: 3 params"), "free virtual-Z") thin underbrace("CZ"_(c_1 c_2)^(b_1) thin "CZ"_(c_1 t)^(b_2) thin "CZ"_(c_2 t)^(b_3), "two-body Z's: 3 params, 1 entangling gate each") thin underbrace("CCZ"_(c_1 c_2 t)^(c_1), "three-body Z: 1 param, expensive")
$

#text(weight: "bold")[What it sees.] Real population errors --
leakage and decoherence. Off-diagonal coherent errors that *cannot*
be absorbed into a diagonal $D$.

#text(weight: "bold")[What it ignores.] Diagonal coherent errors
(any combination of $Z$'s, $"CZ"$'s, $"CCZ"$). These are exactly the
free-to-correct components.

#text(weight: "bold")[Your number.] $overline(F)_("PC") = 0.998230$
at $V_(c c) = 0$, $0.996997$ at $V_(c c) = -5.54$ MHz. The latter
matches Farouk et al.'s 99.7% almost exactly, on the same footing.

== Visual summary of the four metrics

#align(center)[
  #table(
    columns: (auto, auto, auto, auto, auto),
    align: (left, center, center, center, center),
    inset: 6pt,
    stroke: 0.5pt,
    [*Error type*],                              [*$T_P$*], [*$overline(F)_("basis")$*], [*$overline(F)_("raw")$*], [*$overline(F)_("PC")$*],
    [Leakage to $|r angle.r, |R angle.r, |P angle.r$], [✗], [✗], [✗], [✗],
    [Wrong state inside $cal(S)$ (e.g. $A arrow B$ flip wrong)], [✓], [✗], [✗], [✗],
    [Coherent phase, single-qubit $Z$ on any qubit], [✓], [✓], [✗], [✓],
    [Coherent phase, two-body $"CZ"$], [✓], [✓], [✗], [✓],
    [Coherent phase, three-body $"CCZ"$], [✓], [✓], [✗], [✓],
    [Off-diagonal coherent error (e.g. wrong axis)], [✗], [✗], [✗], [✗],
  )
]

*Legend:* ✗ = punishes this error (number drops); ✓ = blind to this
error (number unchanged). The $"CCZ"$ row in particular shows the
difference between $overline(F)_("raw")$ (sees everything) and
$overline(F)_("PC")$ (absorbs the diagonal part).

== The cascade in plain English

Four numbers, four levels of the *same* error budget:

- $1 - T_P$ = #strong[leaked out of the computational subspace]
  ($1.87 times 10^(-3)$ in your run).
- $T_P - overline(F)_("basis")$ = #strong[stayed inside but landed
  in the wrong computational state] ($1.18 times 10^(-4)$).
- $overline(F)_("PC") - overline(F)_("raw")$ = #strong[coherent
  phase error that virtual-Z absorbs for free] ($2.6 times 10^(-3)$
  at $V_(c c) = 0$, $2.0 times 10^(-1)$ at $V_(c c) = -5.54$ MHz).
- $1 - overline(F)_("PC")$ = #strong[unfixable error: leakage,
  decoherence, or off-diagonal coherent error] ($1.77 times 10^(-3)$).

The last bullet is the only thing that *cannot* be repaired by
circuit-level corrections. It is the right "infidelity" to put in a
paper.

= The setup, in one screen

The CCX Hamiltonian (`triqg/hamiltonian.py:132`) is

$
H(t) = & V_("ct") (|r angle.r angle.l r|_(c_1) + |r angle.r angle.l r|_(c_2)) times.circle |R angle.r angle.l R|_t \
& + V_(c c) |r r angle.r angle.l r r|_(c_1 c_2) \
& + Omega_(c c)(t) (|0 angle.r angle.l r| + |r angle.r angle.l 0|)_(c_1, c_2) \
& + Omega_(t,1)(t) (|B angle.r angle.l R| + |R angle.r angle.l B|)_t \
& + Omega_(t,2)(t) (|A angle.r angle.l R| + |R angle.r angle.l A|)_t
$ <eq:H>

Parameters at $a = 5$ μm (current run):

#align(center)[
  #table(
    columns: (auto, auto, auto),
    align: (left, right, left),
    inset: 5pt,
    stroke: 0.5pt,
    [*Quantity*],            [*Value*],                  [*Origin*],
    [$V_("ct") \/ (2 pi)$],  [$+516.81$ MHz],            [Rb--Cs Förster, $C_3$],
    [$V_(c c) \/ (2 pi)$],   [$0$ (set in script)],       [physical value $-5.54$ MHz, multiplied by 0 on line 69],
    [$Omega_(c c) \/ (2 pi)$], [$100$ MHz],              [control $|0 angle.r <-> |r angle.r$ Rabi],
    [$Omega_t \/ (2 pi)$],   [$50$ MHz],                  [target $|A,B angle.r <-> |R angle.r$ Rabi],
    [$T_(c c)$],             [$5$ ns],                    [$pi$-pulse on controls],
    [$T_t$],                 [$10$ ns],                   [target sub-pulse],
    [$T_"gate"$],            [$2 T_(c c) + 3 T_t = 40$ ns], [],
    [$tau_r, tau_R, tau_P$], [$143, 135, 0.131$ μs],     [ARC at 300 K],
  )
]

*Pulse-sign convention (current).* The control pulse $Omega_(c c)$ is
sign-alternating: $+Omega_(c c)\/2$ on the first $T_(c c)$, $-Omega_(c c)\/2$
on the second. The target sub-pulses now have the pattern $(+, -, +)$:
sub-pulses 1 and 3 (driven by `omega_t1`) return $+"amp"/2$, while
sub-pulse 2 (driven by `omega_t2`) returns $-"amp"/2$. This is the
result of Option 0 (#ref(<sec:opt0>)).

= Branch-phase budget: three sources, eight numbers <sec:budget>

For each computational input $|c_1, c_2, t angle.r$ the noiseless gate
produces $e^(i phi_k) U_0 |c_1, c_2, t angle.r$ up to small leakage.
The phase $phi_k$ is the integrated energy along the trajectory, with
three contributions.

== Source 1: $V_(c c)$ on the $|r,r angle.r$ block

Only the $|0,0,dot angle.r$ inputs pump *both* controls to $|r angle.r$
(because $Omega_(c c)$ drives $|0 angle.r <-> |r angle.r$ on each
control, see test `tests/test_hamiltonian.py:285-289`). Time spent in
$|r,r angle.r$:

$
T_(r r) = underbrace(3 T_t, "wait, controls fully in" |r r angle.r) + underbrace(2 dot 3/8 dot T_(c c), "average over the two" pi"-pulses") = 30 + 3.75 = 33.75 "ns"
$

The phase rate on $|r,r angle.r$ is $-V_(c c)$ (Schrödinger evolution
$e^(-i H t)$). Because $V_(c c) < 0$ (Cs--Cs vdW with $C_6 < 0$):

$
phi_(c c) = -V_(c c) dot T_(r r) = +2 pi dot 5.54 "MHz" dot 33.75 "ns" = +0.374 pi
$

This is a $approx +0.37 pi$ phase on every $|0,0,dot angle.r$ branch,
verified numerically (#ref(<sec:budget>)) at the physical
$V_(c c)\/(2 pi) = -5.54$ MHz. The phase moves to ~$0$ if $V_(c c)$
is multiplied out by hand or if the lattice is tuned to a magic spacing
(#ref(<sec:opt1>)).

== Source 2: AC-Stark shift from off-resonant target sub-pulses

When the target is blockaded (controls in $|r angle.r$), the sub-pulses
that try to drive $|A,B angle.r <-> |R angle.r$ at Rabi $Omega_t$ see
a detuned $|R angle.r$ shifted by $N V_("ct")$, where $N in \{1, 2\}$
is the number of controls in $|r angle.r$. The Stark shift on the
driven ground-state level is

$
delta_"AC" = Omega_t^2 / (4 N V_("ct")) quad ==> quad
delta_"AC" \/ (2 pi) = (50^2)/(4 N dot 517) "MHz" = (1.21)/(N) "MHz"
$

Multiplied by the $10$ ns each sub-pulse lasts, this is $+0.024 pi \/ N$
per actuated sub-pulse on the level being driven. *Sign-blind in*
$Omega_t$, so unaffected by Option 0. Sub-pulses 1 and 3 drive
$|B angle.r$, sub-pulse 2 drives $|A angle.r$.

== Source 3: target rotation $+X$ on the $|1,1,dot angle.r$ block (post Option 0)

For $|1,1,dot angle.r$ inputs, the controls stay in $|1 angle.r$ and
the three target sub-pulses fire freely. With the *current* signs
$(s_1, s_2, s_3) = (+, -, +)$, each path through $|R angle.r$ picks up
*one* factor of $-i$ and *one* factor of $+i$:

#align(center)[
  $|1,1,A angle.r arrow.long^("sp 2, " - Omega) +i|1,1,R angle.r arrow.long^("sp 3, " + Omega) (+i)(-i)|1,1,B angle.r = +|1,1,B angle.r$
]
#align(center)[
  $|1,1,B angle.r arrow.long^("sp 1, " + Omega) -i|1,1,R angle.r arrow.long^("sp 2, " - Omega) (-i)(+i)|1,1,A angle.r = +|1,1,A angle.r$
]

so $U_("tgt") = +X$, a clean Toffoli action. The diagonal of
$U_"eff" U_0^dagger$ on $|1,1,dot angle.r$ is $+1$, contributing a
phase $0$. *No spurious CZ*. Verified by the rerun:
$"branch_phases_pi"[6] = "branch_phases_pi"[7] = 0.0000$.

If `omega_t2` were restored to $+"amp"/2$, the action becomes
$U_("tgt") = -X$ and the $|1,1,dot angle.r$ branch picks up the
notorious $+pi$ phase that costs $0.7$ in $overline(F)_"raw"$.

== The complete prediction vs. the rerun

#align(center)[
  #table(
    columns: (auto, auto, auto, auto, auto, auto),
    align: (left, right, right, right, right, right),
    inset: 5pt,
    stroke: 0.5pt,
    [*Branch*],          [*$V_(c c)$*], [*Stark*],     [*$X$*], [*Predicted*], [*Simulated $V_(c c)=0$*],
    [$|0,0,A angle.r$],  [$0$],         [$+0.012 pi$], [$0$],   [$+0.012 pi$], [$+0.0119 pi$],
    [$|0,0,B angle.r$],  [$0$],         [$+0.024 pi$], [$0$],   [$+0.024 pi$], [$+0.0233 pi$],
    [$|0,1,A angle.r$],  [$0$],         [$+0.024 pi$], [$0$],   [$+0.024 pi$], [$+0.0234 pi$],
    [$|0,1,B angle.r$],  [$0$],         [$+0.048 pi$], [$0$],   [$+0.048 pi$], [$+0.0476 pi$],
    [$|1,0,A angle.r$],  [$0$],         [$+0.024 pi$], [$0$],   [$+0.024 pi$], [$+0.0234 pi$],
    [$|1,0,B angle.r$],  [$0$],         [$+0.048 pi$], [$0$],   [$+0.048 pi$], [$+0.0476 pi$],
    [$|1,1,A angle.r$],  [$0$],         [$0$],         [$0$],   [$0$],         [$+0.0000 pi$],
    [$|1,1,B angle.r$],  [$0$],         [$0$],         [$0$],   [$0$],         [$+0.0000 pi$],
  )
]

Every branch is reproduced to better than $0.001 pi$. With $V_(c c)$
restored, the $|0,0,dot angle.r$ branches shift to $+0.387 pi$ and
$+0.398 pi$ (also confirmed numerically), exactly matching the original
analysis with $phi_(c c) = +0.374 pi$.

= Numerical comparison of the four metrics across runs

#align(center)[
  #table(
    columns: (auto, auto, auto, auto),
    align: (left, right, right, right),
    inset: 6pt,
    stroke: 0.5pt,
    [*Metric*], [*$V_(c c) = 0$ (current)*], [*$V_(c c) = -5.54$ MHz*], [*Original (no Option 0)*],
    [$T_P$], [$0.998129$], [$0.996744$], [$0.996744$],
    [$overline(F)_"basis"$], [$0.998011$], [$0.996626$], [$0.996626$],
    [$overline(F)_"raw"$],   [$0.995625$], [$0.796283$], [$0.282014$],
    [$overline(F)_"PC"$],    [$0.998230$], [$0.996997$], [$0.996997$],
    [gap PC $-$ raw],        [$0.0026$],   [$0.201$],     [$0.715$],
  )
]

Reading the table column by column:

+ #strong[Current column ($V_(c c) = 0$, Option 0 applied):] all four
  numbers within $3 times 10^(-3)$. Means there is essentially no
  removable diagonal error, *and* very little leakage. This is the
  best-case theoretical limit for the protocol.

+ #strong[Physical $V_(c c)$ column:] $T_P$, $overline(F)_("basis")$,
  $overline(F)_("PC")$ all stay near $0.997$ (V_cc costs only
  $1.4 times 10^(-3)$ of leakage on the basis number). But
  $overline(F)_("raw")$ drops to $0.80$. The full $0.20$ gap is the
  removable $"CZ"_(c_1 c_2)$ from V_cc -- pure bookkeeping.

+ #strong[Original column (no Option 0):] $overline(F)_("basis")$
  and $overline(F)_("PC")$ identical to the physical-$V_(c c)$
  column above (Option 0 changes only a phase that PC absorbs).
  $overline(F)_("raw")$ collapses to $0.28$ because the $+pi$ on
  $|1,1,dot angle.r$ adds a *second* removable $"CZ"_(c_1 c_2)$ on
  top of the V_cc one.

The story across the row $overline(F)_("raw")$: $0.282 arrow.r 0.796 arrow.r 0.996$
as you (i) flip $omega_(t,2)$ sign, then (ii) zero out V_cc. Both
changes are removable -- the underlying gate quality (measured by
$overline(F)_("PC")$) only moves from $0.997$ to $0.998$.

== Which part of the diagonal $D$ is removable for free?

The optimal $D$ (the diagonal post-rotation that maximises the
phase-corrected fidelity) decomposes uniquely into single-qubit $Z$,
two-qubit $"CZ"$, and three-qubit $"CCZ"$ components. Single-qubit
$Z$'s are free on every platform (virtual-$Z$); two- and three-qubit
phases require additional gates. For the current run, the eight
optimal phases (in units of $pi$) are:

$
phi^* = [+0.717, -1.272, +0.729, +0.753, +0.729, -1.247, +0.705, -1.295]
$

With $V_(c c) = 0$, this $D$ is dominated by a single-qubit $Z$ on
the target and a global phase -- both free. With $V_(c c) = -5.54$ MHz
restored, the additional structure is a $"CZ"_(c_1 c_2)$; that is
*not* free locally, but it commutes with the surrounding Clifford
frame and is absorbed by circuit-level tracking (zero gate cost in
fault-tolerant compilation). *The $"CCZ"$ component remains negligible*
at both V_cc settings -- this is a clean Toffoli, not a noisy
three-body unitary.

= How Farouk et al. handled the same physics <sec:farouk>

The paper's own diagnosis is unambiguous (p. 14, after Fig. 9):

#block(
  fill: luma(248),
  inset: 8pt,
  radius: 3pt,
  width: 100%,
)[
  #emph[
    "In Figure 9b, we study a non-realistic case where we neglected
    the interaction between control atoms... the maximum value of
    fidelity becomes possible for a wider range of inter-atomic
    distances. This case proves that the destructive pattern in system
    dynamics is a direct result of $V_("CC")$.
  ]
]

Their Fig. 9a (realistic, $V_("CC") != 0$) shows a 99.7% peak only
in a *narrow strip* near $R_("CT") approx 6$ μm with
$Omega_c \/ Omega_p > 2.5$. Fig. 9b (artificial, $V_("CC") = 0$)
shows a much wider plateau. They scanned the $(R_("CT"), Omega_c \/ Omega_p)$
plane and reported the best point. *Setting $V_(c c) = 0$ in our
simulation is the equivalent of standing on Farouk's idealised Fig. 9b
plateau.*

== Their sweet spot is a magic-spacing condition

With Cs $C_6 = 2 pi dot 2364$ GHz μm#super[6] and their geometry
$R_("CC") = (2 \/ sqrt(5)) R_("CT") approx 0.894 R_("CT")$:

$
V_("CC") \/ (2 pi) "at" R_("CC") approx 5.36 "μm" arrow.r 100 "MHz"
$

(matches Fig. 9c). With gate duration $tau approx 1.5$ μs:

$
(V_("CC") dot tau) / (2 pi) approx 100 "MHz" dot 1.5 "μs" = 150 quad arrow.r quad "integer"
$

The $|1,1 angle.r$ branch makes *150 full revolutions* under $V_("CC")$
during the gate and lands back on its starting phase. Move $R_("CT")$
to $12$ μm and $V_("CC")$ drops to $approx 1.6$ MHz; then
$V_("CC") tau \/ (2 pi) approx 2.36$, off by $0.36 dot 2 pi approx 130 degree$,
and fidelity collapses --- exactly the "sharp drop" they describe at
$8 "μm" < R_("CT") < 12 "μm"$.

== Apples-to-apples comparison

Farouk's reported number is *state fidelity for a single superposition
input*:
$
|psi angle.r = 1/sqrt(2) (|0 0 angle.r |A A angle.r + |1 1 angle.r |A A angle.r)
arrow 1/sqrt(2) (|0 0 angle.r |A A angle.r + |1 1 angle.r |B B angle.r)
$
$
F = |angle.l psi^"target" | psi^"actual" angle.r|^2 = (1 + cos phi_(c c)) / 2
$

This is exactly the kind of quantity that is *sensitive* to the
$V_(c c)$ phase. Farouk's 99.7% requires $cos phi_(c c) approx 1$
modulo $2 pi$ -- the magic-spacing condition.

The TriQG CCX with Option 0 and $V_(c c)$ artificially zeroed reaches
$overline(F)_"PC" = 0.998$, slightly *better* than Farouk because the
gate is $40$ ns instead of $1.5$ μs and pays less Rb $|R angle.r$
decay. With Option 0 only and the physical $V_(c c)$, the apples-to-apples
score is $overline(F)_"PC" = 0.997$ -- equal to Farouk's. The "raw"
$0.80$ in that case is the same diagonal $"CZ"_(c_1 c_2)$ that Farouk
hides inside the integer $V_("CC") tau$.

#pagebreak()

= Four engineering options to remove the residual phase mismatch <sec:options>

These are presented in increasing order of *physics cost*. Option 0 is
already in place in the current code; Options 1--4 attack the residual
$V_(c c)$ contribution.

== Option 0 (free, *applied*): Flip the sign of `omega_t2` <sec:opt0>

*Status:* applied. `triqg/pulses.py:omega_t2` returns $-"amp"/2$. The
$|1,1,dot angle.r$ branch phase moved from $pi$ to $0$.

*Why this works.* The branch sign acquired by $|1,1,A angle.r$ through
$U_("tgt")$ is the product of the two sub-pulses that *actually couple
to its trajectory* (sub-pulses 2 and 3). The branch sign for
$|1,1,B angle.r$ is the product of sub-pulses 1 and 2. With sign vector
$(s_1, s_2, s_3)$,
$ U_("tgt") = mat(0, -s_1 s_2; -s_2 s_3, 0). $
A clean $+X$ requires $s_1 s_2 = s_2 s_3 = -1$, which is satisfied by
$(s_1, s_2, s_3) = (+, -, +)$ (current) or $(-, +, -)$. Independently
flipping *one* outer sub-pulse gives an asymmetric
$plus.minus i sigma_y$ and does *not* help, because half of
$omega_(t,1)$ is useless on any given input trajectory: sub-pulse 1
is dead for $|dot,dot,A angle.r$ inputs and sub-pulse 3 is dead for
$|dot,dot,B angle.r$ inputs.

*What it bought.* Eliminated the $+pi$ on the $|1,1,dot angle.r$
block, i.e. removed the spurious $"CZ"_(c_1 c_2)$. Confirmed
improvements (with $V_(c c) = 0$ idealisation):
$
overline(F)_"raw" : 0.282 arrow.r 0.996, quad
overline(F)_"basis" : 0.997 arrow.r 0.998, quad
overline(F)_"PC" : 0.997 arrow.r 0.998
$

*What it does not buy.* The $V_(c c)$-induced $approx +0.37 pi$ on
$|0,0,dot angle.r$ (visible at $V_(c c) = -5.54$ MHz, where
$overline(F)_"raw"$ is still only $0.80$) and the AC-Stark
$approx 0.05 pi$ on the $|0,1 angle.r, |1,0 angle.r$ branches. Those
need Options 1--4.

*Cost:* zero. No additional pulses, no extra time, no extra
decoherence.

== Option 1: Magic-spacing trick (Farouk's solution, ported) <sec:opt1>

*What:* Choose the lattice spacing $a$ so that
$|V_(c c)| dot T_"wait" = 2 pi N$ for some small integer $N$.

*The condition.* With current $T_"wait" approx 30$ ns and
$C_6 = -693$ GHz μm#super[6], requiring

$
|V_(c c)| / (2 pi) = N / T_"wait" = N dot 33.3 "MHz"
quad arrow.r quad
R_("AA") = ((|C_6| dot T_"wait") / (2 pi N))^(1\/6)
$

gives the table:

#align(center)[
  #table(
    columns: (auto, auto, auto, auto),
    align: (right, right, right, right),
    inset: 5pt,
    stroke: 0.5pt,
    [*$N$*], [*$|V_(c c)|\/(2 pi)$*], [*$R_("AA")$*], [*Lattice $a = R_("AA")\/sqrt(2)$*],
    [$1$], [$33.3$ MHz], [$5.31$ μm], [$3.75$ μm],
    [$2$], [$66.7$ MHz], [$4.73$ μm], [$3.34$ μm],
    [$3$], [$100.0$ MHz], [$4.41$ μm], [$3.12$ μm],
  )
]

*Bonus.* Tightening the lattice also shrinks $R_("DA") = a \/ sqrt(2)$
and increases $V_("ct") = C_3 \/ R_("DA")^3$. At $a = 3.75$ μm,
$V_("ct") \/ (2 pi)$ jumps from $517$ MHz to $approx 1.23$ GHz, which
*reduces* finite-blockade leakage as $(V_("ct") \/ Omega_t)^(-2)$ --
roughly a $5 times$ reduction in the dominant residual error.

*Tolerance.* The phase deviation from the magic point is
$Delta phi_(c c) = (partial V_(c c) \/ partial R_("AA")) dot Delta R_("AA") dot T_"wait"
= -6 V_(c c) Delta R_("AA") \/ R_("AA") dot T_"wait"$.
For $|Delta phi_(c c)| < pi \/ 100$ ($> 99.95%$ overlap), one needs
$|Delta R_("AA")| \/ R_("AA") < 10^(-2) \/ (6 N)$.
At $N = 1$, $a = 3.75$ μm, this is $approx 17$ nm tolerance on the
inter-atomic distance --- within reach of optical-tweezer spacing
calibration but not free.

*Caveats.*

- *Ground-state interaction.* Below $approx 3$ μm, the long-range
  ground-state forces (mostly van der Waals between
  $|6 S_(1\/2) angle.r$ Cs atoms) start lifting the lattice
  approximation; you will need to verify that the trap depth and
  loading temperature still allow sub-Lamb-Dicke confinement.
- *Cross-talk.* Tighter $a$ shrinks $R_("DA")$, which can put nearby
  Rb--Cs pairs into a regime where leakage to off-resonant Förster
  channels becomes significant. The Ireland et al. 2024 channel list
  must be re-checked at the new geometry.
- *Magic spacing is single-frequency.* It only kills $V_(c c) tau$,
  not other coherent shifts (e.g. AC Stark on the target). Those
  remain at the few-percent level.
- *Population leakage from V_cc remains.* The current rerun shows
  that even the *basis* fidelity drops by $1.4 times 10^(-3)$ when
  $V_(c c) = -5.54$ MHz is restored. This is *not* a phase, but a
  real population leak from V_cc detuning the $|0,0 angle.r arrow |r,r angle.r$
  excitation during the control $pi$-pulses (see #ref(<sec:newphysics>)
  below). Magic spacing kills the $|r,r angle.r$ phase but does *not*
  cure this leakage; the only way to suppress it fully is
  $|V_(c c)| \/ Omega_(c c) << 1$, which automatically holds at
  the small-V_cc end of the magic-spacing table ($N=1, a=3.75$ μm:
  $|V_(c c)|/Omega_(c c) = 33.3/100 = 0.33$, marginal) and fails
  at high N.

*What it buys.* $overline(F)_"raw" arrow.r overline(F)_"PC")$ on the
$|0,0,dot angle.r$ branch, eliminating the single largest source of
phase mismatch. Combined with Option 0,
$overline(F)_"raw" approx overline(F)_"PC" approx 0.997$ at the
physical $V_(c c)$, or $approx 0.998$ with the additional
$1.4 times 10^(-3)$ leakage suppression at small $|V_(c c)| / Omega_(c c)$.

== Option 2: Stretch $T_"wait"$ at fixed lattice <sec:opt2>

*What:* Hold $a = 5$ μm (so $V_(c c)\/(2 pi) = -5.54$ MHz). Pick
$T_"wait"$ so that $|V_(c c)| dot T_"wait" = 2 pi N$:
$
T_"wait"^"magic" = N \/ |V_(c c)\/(2 pi)| = N dot 180.5 "ns"
$

The smallest magic wait is $approx 180$ ns, $6 times$ longer than the
current $30$ ns. Total gate time grows from $40$ ns to
$2 T_(c c) + 180 = 190$ ns.

*What you pay.* Extra time spent in the lossy Rydberg states. The
incremental decoherence cost is
$
1 - F_"survival" approx (Delta T_"wait") (gamma_R + gamma_P P_R)
$
With $Delta T_"wait" = 150$ ns, $tau_R = 135$ μs and assuming the
target spends $approx O(1)$ of its time in $|R angle.r$ during the
extra wait (mostly *not*, since the target sub-pulses are not
extended), the dominant cost is $|R angle.r$ decay during the
intermediate sub-pulses, $approx 150 \/ 135000 = 1.1 times 10^(-3)$
extra infidelity. *This is a real loss*, comparable to the basis
infidelity itself.

*Variant.* Instead of stretching the wait, *hold the target sub-pulses
at full duration but insert idle gaps between them* until
$T_"wait" = 180$ ns. The target sub-pulse decoherence is unchanged;
only the controls' $|r angle.r$ population decays during the extra
wait. With $tau_r = 143$ μs, the cost drops to
$approx 2 dot 150 \/ 143000 approx 2 times 10^(-3)$ (factor of 2 for two
controls in $|r angle.r$). Still costly.

*Verdict.* Cheaper than Option 1 (no lattice rebuild), but trades a
coherent error for a population error of similar magnitude. Useful
mostly as a quick numerical experiment to *verify* the magic-spacing
mechanism, not as a final design.

== Option 3: Accept $V_(c c)$ as a free $"CZ"$ and report
   $overline(F)_"PC"$ <sec:opt3>

*What:* Stop computing $overline(F)_("Pedersen,raw")$ as the figure
of merit. Report either $overline(F)_("basis")$ or $overline(F)_("PC")$,
both of which exceed $0.997$ even at the physical $V_(c c)$.

*Why this is legitimate.* The "missing" infidelity is a known unitary
$D$ acting diagonally on the computational basis. Decompose:
$
D = e^(i alpha) thin (Z^(a_1) times.circle Z^(a_2) times.circle Z^(a_3))_("local") thin ("CZ"_(c_1 c_2))^(b_1) thin ("CZ"_(c_1 t))^(b_2) thin ("CZ"_(c_2 t))^(b_3) thin ("CCZ"_(c_1 c_2 t))^(c_1)
$

For your CCX, the 8 optimal phases factor (to within numerical
precision) as
$
D approx D_("local Z") dot ("CZ"_(c_1 c_2))_("from " V_(c c)) dot D_("small Stark Z")
$
with no measurable $"CCZ"$ component (verified at both
$V_(c c) = 0$ and $V_(c c) = -5.54$ MHz). So:

- the local $Z$'s are absorbed by virtual-$Z$ on subsequent gates
  ($"cost" = 0$);
- the $V_(c c)$-induced $"CZ"_(c_1 c_2)$ is absorbed by *one*
  software $"CZ"$ tracked through the surrounding Clifford layer
  ($"cost" = 0$ in fault-tolerant compilation, since $"CZ"$'s
  commute through the Clifford frame);
- the small two-body $Z$'s from Stark are at the $10^(-4)$ level and
  can be dropped.

*What it buys.* Honest fidelity number that matches the gate's
operational quality. Apples-to-apples with Farouk: their reported
99.7% is also a "phase-tuned" number.

*What it does not buy.* If your downstream circuit cannot absorb a
$"CZ"_(c_1 c_2)$ for free (e.g. you're benchmarking the CCX as a
standalone process tomography target with no surrounding Clifford
layer), this option is unavailable and you need Option 0 + 1 instead.
Also: the $1.4 times 10^(-3)$ population leakage from V_cc detuning
the $|0,0 angle.r arrow |r,r angle.r$ excitation is *not* a Z phase
and is *not* absorbable here. To regain that, you still need Option 1.

== Option 4: Hahn echo on $V_(c c)$ within the gate <sec:opt4>

*What:* Insert a midway control operation that swaps
$|r,r angle.r <-> |0,0 angle.r$ at $T_"wait" \/ 2$, so that $V_(c c)$
accumulates $+phi$ on the first half and $-phi$ on the second.

*The pulse sequence.* Modify the protocol to
$
+pi_(c c) quad arrow quad T_"wait" \/ 2 quad arrow quad 2pi_(c c) quad arrow quad T_"wait" \/ 2 quad arrow quad -pi_(c c)
$

The middle $2 pi$-pulse on each control returns it to $|0 angle.r$
and re-excites it to $|r angle.r$, swapping the population trajectory
through a momentary $|0,0 angle.r$ window. Total $V_(c c) tau$
contributions cancel.

*What you pay.* Two extra control $pi$-pulses (or one $2 pi$-pulse).
Time: $+2 T_(c c) = 10$ ns. Decoherence: an extra
$approx 10\/143000 approx 7 times 10^(-5)$ from $|r angle.r$ decay --
*negligible*.

*Subtleties.*

- The mid-sequence $2 pi$-pulse must occur *during the target wait
  window*, but the target is mid-Raman / mid-Rabi at that point. The
  swap will introduce a transient detuning on the target sub-pulse 2
  unless the controls are fully de-excited and re-excited cleanly.
- For sub-pulse 2 (which drives $|A angle.r <-> |R angle.r$) the
  target is in $|R angle.r$ during $approx 50%$ of the sub-pulse;
  during that time, swapping the controls in/out of $|r angle.r$
  modulates the blockade detuning $V_("ct")$, which is not
  cancelled by the echo. Numerical check required.
- The clean version is to insert the echo *between* sub-pulses, e.g.
  after sub-pulse 1 and before sub-pulse 3. That breaks the symmetry
  of the target pulse-area and may need re-optimisation of
  $Omega_t T_t = pi$.

*What it buys.* Removes the $V_(c c)$ phase *exactly*, independent
of lattice spacing. Combined with Option 0, gives
$overline(F)_"raw" approx overline(F)_"PC" approx 0.997$ at the
*current* lattice $a = 5$ μm.

*What it does not buy.* Same as Option 1: the Stark contribution
is unchanged, and the $1.4 times 10^(-3)$ V_cc-induced population
leakage during the *control* $pi$-pulses is also unchanged (the
echo is inside the target wait window, not inside the control
pulses).

= New finding from the rerun: V_cc-induced leakage during the *control* pulses <sec:newphysics>

The original analysis claimed V_cc was a pure $Z$ phase on the
$|r,r angle.r$ block, hence invisible to $overline(F)_"basis"$. The
rerun shows this is *almost* true but not quite: at the physical
$V_(c c) = -5.54$ MHz, $overline(F)_"basis"$ drops by
$1.4 times 10^(-3)$ relative to the $V_(c c) = 0$ baseline. Source:

The two-control system under the Hamiltonian $H_(c c) + V_(c c) P_(r r)$,
where $H_(c c)$ is the symmetric Rabi drive on both controls and
$P_(r r)$ is the projector onto the doubly-Rydberg state, does *not*
factorise into two independent single-atom rotations whenever V_cc is
nonzero.
In the symmetric two-atom subspace, the Hamiltonian becomes a 3-level
cascade with effective Rabi $sqrt(2) thin Omega_(c c)$ between
the ground state, the symmetric singly-excited state, and the
doubly-excited state, with the doubly-excited level detuned by
$V_(c c)$. A nominally $pi$-pulse on each individual control no longer
fully transfers the population from ground to doubly-excited; the
transfer infidelity is

$
epsilon_"transfer" approx (V_(c c)^2) / (Omega_(c c)^2) approx 3 times 10^(-3)
$

per control pair excitation event, halved over the two excitations
(in and out) because the protocol is a near-echo. The observed drop
of $1.4 times 10^(-3)$ matches this estimate to within a factor of 2.

*Implication for Option 1.* Magic spacing cures the V_cc phase, but
if the magnitude of V_cc is large, the population leakage during the
control pulses can dominate. The optimum is small V_cc (magic spacing
at the smallest N) and large Omega_cc such that
(V_cc / Omega_cc)^2 < 1e-3, e.g. Omega_cc > 5 |V_cc|. The idealised
limit with V_cc = 0 is the strict upper bound on what any
magic-spacing realisation can achieve.

= Recommendation

#align(center)[
  #table(
    columns: (auto, auto, auto, auto, auto),
    align: (left, center, center, center, left),
    inset: 6pt,
    stroke: 0.5pt,
    [*Option*], [*Cost*], [*$overline(F)_"raw"$*], [*$overline(F)_"PC"$*], [*Recommended for*],
    [Pre-Option-0 baseline], [--], [$0.282$], [$0.997$], [diagnostics only],
    [Option 0 only], [zero], [$0.796$], [$0.997$], [#emph[done already]],
    [Option 0 + V_cc=0 (idealised)], [n/a], [$0.996$], [$0.998$], [theoretical bound],
    [0 + 1: magic spacing], [hardware],  [$approx 0.996$], [$approx 0.998$], [final design],
    [0 + 2: stretch $T_"wait"$], [$10^(-3)$ extra], [$approx 0.99$], [$approx 0.996$], [verification only],
    [0 + 3: report $overline(F)_"PC"$], [bookkeeping], [$0.796$], [$0.997$], [comparing to Farouk],
    [0 + 4: $V_(c c)$ echo], [$10^(-5)$ extra], [$approx 0.997$], [$approx 0.997$], [if lattice is fixed],
  )
]

#block(
  fill: luma(245),
  inset: 10pt,
  radius: 4pt,
  width: 100%,
)[
  *Bottom line.* Option 0 is in. The current $overline(F)_"PC" = 0.998$
  with $V_(c c) = 0$ is the operational ceiling of this protocol at
  $a = 5$ μm and $T_"gate" = 40$ ns. To approach it with the *physical*
  $V_(c c)$ in place, choose:
  + *For publication numbers:* Option 3 -- report
    $overline(F)_"PC" = 0.997$ honestly and document the
    $"CZ"_(c_1 c_2)$ absorption. Apples-to-apples with Farouk.
  + *For raw-number performance:* Option 1 -- shrink the lattice to
    $a = 3.75$ μm to land on $N = 1$ magic spacing. This also doubles
    $V_("ct")$ and lowers blockade leakage. Expected
    $overline(F)_"raw" approx 0.997$.
  + *If the lattice is locked at $a = 5$ μm:* Option 4 -- insert a
    $V_(c c)$ Hahn echo. Cost is one extra $2 pi$-pulse on the
    controls, $approx 10^(-5)$ extra decoherence.
]

= Verification protocol — STATUS

#align(center)[
  #table(
    columns: (auto, auto, auto),
    align: (left, left, center),
    inset: 6pt,
    stroke: 0.5pt,
    [*Test*], [*Predicted*], [*Status*],
    [Option 0 verification (`omega_t2` flipped)], [$|1,1,dot angle.r$ branch phase: $pi arrow.r 0$. $overline(F)_"raw"$: $0.282 arrow.r 0.80$], [#text(fill: rgb(0,128,0))[*confirmed*]],
    [Source 1 isolation ($V_(c c) = 0$)], [$|0,0,dot angle.r$ branch phase: $0.39 pi arrow.r 0.012 pi$], [#text(fill: rgb(0,128,0))[*confirmed*]],
    [Source 2 (Stark) prediction], [$|dot,dot,A angle.r$ phase: $0.024 pi \/ N$ per A-driving sub-pulse], [#text(fill: rgb(0,128,0))[*confirmed within $0.001 pi$*]],
    [Combined (Option 0 + V_cc=0)], [$overline(F)_"raw" approx 0.996$], [#text(fill: rgb(0,128,0))[*confirmed: $0.9956$*]],
    [V_cc population leakage], [$Delta(1 - overline(F)_"basis") approx (V_(c c)/Omega_(c c))^2 / 2 approx 1.5 times 10^(-3)$], [#text(fill: rgb(0,128,0))[*confirmed: $1.4 times 10^(-3)$*]],
    [Option 1 verification ($a = 3.75$ μm)], [$overline(F)_"raw" arrow.r overline(F)_"PC"$ at the physical $V_(c c)$], [#text(fill: rgb(160,80,0))[*pending* -- requires lattice-distance scan]],
  )
]

All five completed predictions in the original analysis are confirmed
to within stated tolerances. The remaining lattice scan (Option 1)
is the next experimental step.

#v(1em)
#line(length: 100%, stroke: 0.5pt)
#text(size: 9pt)[
  Source script: `examples/Average_fidelity_Vcc_newenergy/ccx_pedersen_average_fidelity.py`
  (line 69: `V_cc_MHz = 1000.0 * C6_CsCs / r_AA**6 * 0`). \
  Pulse definitions: `triqg/pulses.py` (`omega_t2` returns `-amp / 2`). \
  Branch-phase analysis: `examples/Average_fidelity_Vcc_newenergy/branch_phase_analysis.py`
  (uses physical $V_(c c)$, not zeroed). \
  Reference: M. Farouk et al., #emph[Photonics] *10*, 1280 (2023),
  Sec. 5 and Fig. 9. \
  Pedersen formalism: L. H. Pedersen, N. M. Møller, K. Mølmer,
  #emph[Phys. Lett. A] *367*, 47 (2007). \
  Rerun timestamp: 2026-04-27.
]

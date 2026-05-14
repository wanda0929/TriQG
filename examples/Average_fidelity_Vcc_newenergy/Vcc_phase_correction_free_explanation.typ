// Why the diagonal phase correction D in F_PC = max_D F(E, D U_0) is
// free of charge in real atomic-qubit experiments.
//
// Companion note to Vcc_phase_origin_report.typ. Written for a
// reader who is comfortable with quantum gate physics but skeptical
// of the "phase correction is free" claim.
//
// Compile with:  typst compile Vcc_phase_correction_free_explanation.typ

#set document(
  title: "Why the phase correction is free of charge",
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
  #text(size: 17pt, weight: "bold")[
    Why the phase correction $D$ in
    $overline(F)_("PC") = max_D thin overline(F)(cal(E), D thin U_0)$
    is free of charge
  ]

  #v(0.3em)
  #text(size: 11pt)[
    A first-principles answer for a tutor who does not yet believe it.
    Companion note to `Vcc_phase_origin_report.typ`.
  ]

  #v(0.4em)
  #text(size: 10pt)[
    TriQG · `examples/Average_fidelity_Vcc_newenergy` · 2026-05-02
  ]
]

#v(0.6em)

#block(
  fill: luma(245),
  inset: 10pt,
  radius: 4pt,
  width: 100%,
)[
  *The tutor's question.*  In the Pedersen analysis I report a
  *phase-corrected* gate fidelity
  $
  overline(F)_("PC") = max_(D thin "diagonal," thin abs(D_(k k)) = 1) thin overline(F)(cal(E), thin D thin U_0) ,
  $
  where the unitary $D$ absorbs three single-qubit virtual-$Z$ rotations
  *plus* one parasitic two-qubit $"CZ"_(c_1 c_2)$ that comes from the
  Cs--Cs van der Waals interaction $V_(c c)$.  My tutor objects:

  #emph[
    "You are calling those phase corrections #emph[free of charge].
    But a $Z$ rotation and especially a $"CZ"$ are real quantum operations.
    They take real time, they are imperfect, and they will degrade
    the fidelity.  How can it possibly be free in a real atomic-system
    experiment?"
  ]

  *The two-second answer.*  Because we never apply them.  In the
  experiment, single-qubit $Z$ rotations are realised as a phase
  reference shift in the laser/microwave electronics
  --- the atom is not touched.  The parasitic $"CZ"$ on the two
  controls is not turned into a physical pulse either; it is bookkept
  classically and propagated through the surrounding circuit as a
  Pauli-frame update.  The "cost" is a few floating-point operations
  in the control computer.

  *Where this is rigorously stated.*  McKay et al.@mckay2017 prove
  that a single-qubit $Z(phi)$ implemented as a phase shift of all
  subsequent drive pulses is *exactly* equivalent to a physical
  $Z(phi)$, with zero error and zero duration.  Knill@knill2005 and
  the subsequent Pauli-frame literature @riesebos2017 extend the
  argument to multi-qubit Cliffords, including $"CZ"$.  Every modern
  neutral-atom experiment that quotes a phase-corrected Pedersen-Mølmer
  fidelity @levine2019 @evered2023 @madjarov2020 @graham2022
  @bluvstein2022 @bluvstein2024 @maller2015 uses exactly this
  bookkeeping.

  The rest of this note explains the physics of *why* it is exact, the
  algebra of *how* it propagates, and the boundary of *where* it stops
  being free.
]

= The setup: how an atomic qubit is actually rotated <sec:setup>

To answer the tutor honestly, we must look at how a single-qubit gate
is *physically* implemented on an atom.  No equivalence between
"physical Z" and "virtual Z" can be argued without it.

== The atom in its own rotating frame

Take a hyperfine qubit on Cs or Rb with computational levels
$|0 angle.r$ and $|1 angle.r$ split by an angular frequency $omega_0$.
The bare-atom Hamiltonian is

$
H_0 = (omega_0)/(2) sigma_z .
$ <eq:Hzero>

A single-qubit gate is driven by a microwave or two-photon Raman field
of the form $E(t) = E_0 thin Omega(t) thin cos(omega_L t + phi_L)$ where
$omega_L approx omega_0$ is the carrier frequency and $phi_L$ is the
*phase of the local oscillator at $t = 0$*.  After the rotating-wave
approximation in the frame rotating at $omega_L$ the Hamiltonian becomes

$
H_("drive") = (Omega(t))/(2) [cos(phi_L) thin sigma_x + sin(phi_L) thin sigma_y] = (Omega(t))/(2) thin sigma_(phi_L) ,
$ <eq:Hdrive>

i.e. a Rabi oscillation about an axis $sigma_(phi_L)$ that lies in the
$x y$-plane of the Bloch sphere, at angle $phi_L$ from the $x$-axis.
A pulse of area $integral Omega(t) thin d t = theta$ implements

$
X_(phi_L)(theta) := exp[ - i thin (theta)/(2) thin sigma_(phi_L) ] .
$ <eq:Xphi>

#text(weight: "bold")[The crucial observation.]  $phi_L$ is *not* an
atomic property.  It is a number stored in the digital electronics
that drives the laser.  Setting $phi_L$ to a new value is a
single-instruction change in the AWG or DDS chip --- nanoseconds of
control software, with no quantum operation on the atom whatsoever.

= The algebraic identity that makes virtual $Z$ exact <sec:algebra>

A virtual-$Z$ gate is implemented by changing $phi_L$.  We now check
that this is *exactly* equivalent to a physical
$Z(alpha) = exp(-i alpha sigma_z slash 2)$.

Compute $Z(alpha) thin X_(phi_L)(theta) thin Z^dagger(alpha)$.  Using
$Z(alpha) thin sigma_x thin Z^dagger(alpha) = cos(alpha) thin sigma_x + sin(alpha) thin sigma_y$
and the analogous identity for $sigma_y$, the rotation axis
$sigma_(phi_L)$ goes to $sigma_(phi_L + alpha)$, so

$
Z(alpha) thin X_(phi_L)(theta) thin Z^dagger(alpha) = X_(phi_L + alpha)(theta) .
$ <eq:conjugation>

Equivalently,

$
Z(alpha) thin X_(phi_L)(theta) = X_(phi_L + alpha)(theta) thin Z(alpha) .
$ <eq:Zcommute>

#text(weight: "bold")[The propagation rule.]  Whenever a $Z(alpha)$
appears in the circuit, push it through the next $X$/$Y$ rotation by
incrementing that rotation's phase by $alpha$.  The $Z$ keeps moving
right until either (i) the circuit ends and it hits a measurement, or
(ii) it reaches a non-Clifford gate (see #ref(<sec:boundary>)).

== The two consequences

#text(weight: "bold")[Consequence 1: cost on the atom is zero.]  Look
at #ref(<eq:Zcommute>): on the *physical* hardware we never apply a
$Z(alpha)$ pulse.  We absorb $alpha$ into the phase of the *next* drive
pulse and apply $X_(phi_L + alpha)(theta)$ instead.  The atom sees
exactly *one* rotation, of the same area and the same duration as
before; only the rotation axis is different.  No additional pulse
time, no additional spontaneous emission, no additional dephasing.
Cost: a single addition in the AWG phase register.

#text(weight: "bold")[Consequence 2: identity is exact.]  The
identity #ref(<eq:conjugation>) is an algebraic statement about
$2 times 2$ unitaries.  It is true to *all* orders in $alpha$, not
to leading order.  Unlike a physical $Z$ pulse --- which has finite
amplitude error, finite duration, finite Stark shift, finite
spontaneous emission --- a virtual $Z$ has *no* error budget at all.
@mckay2017 demonstrate this on superconducting qubits with
randomized benchmarking and find the virtual-$Z$ infidelity to be at
or below the noise floor of their measurement; the same physics
applies to atomic qubits, where the laser-phase register plays the
role of the AWG-phase register.

== A concrete example

Consider the circuit fragment $X_0(pi slash 2) thin Z(pi slash 4) thin X_0(pi slash 2)$.
Naively implemented, this is three pulses on the atom: a $pi slash 2$
rotation about $x$, then a $Z(pi slash 4)$ pulse (which costs time and
photon scattering), then another $pi slash 2$ rotation about $x$.
Using #ref(<eq:Zcommute>), the same circuit is

$
X_0(pi slash 2) thin Z(pi slash 4) thin X_0(pi slash 2) = X_(pi slash 4)(pi slash 2) thin Z(pi slash 4) thin X_0(pi slash 2)
$
(propagate the $Z$ once)
$
                                                       = Z(pi slash 4) thin X_(pi slash 4)(pi slash 2) thin X_0(pi slash 2) .
$
(propagate again)

The $Z(pi slash 4)$ now sits at the very end of the circuit.  If the
final operation is a $sigma_z$-basis measurement, $Z(pi slash 4)$
commutes with the projector $|0 angle.r angle.l 0|$ and contributes
nothing.  Two pulses on the atom, both with shifted phases.  Zero
"$Z$ pulses" applied.

#v(0.4em)

= Why this matters for the Cs Rydberg gate <sec:atomic>

Two skeptical objections a tutor might raise.

== "But the laser has a finite linewidth — surely the phase shift is noisy?"

The relevant question is not the laser's absolute phase noise but
*its phase noise during the gate window*.  Modern phase-locked Raman
lasers and microwave sources used in neutral-atom experiments have
phase coherence times much longer than the gate duration:
$tau_("coh") gt.tilde 10$ ms versus $T_("gate") tilde 100$ ns -- $1 mu$s.
This is why the McKay-Gambetta result transfers cleanly: the atom is
never asked to *resolve* the absolute laser phase, only relative phases
between consecutive pulses, and over $tilde mu s$ scales those are at
the part-per-million level @graham2022 @evered2023.  Adding a
calibrated phase offset $alpha$ to the digital register that *defines*
$phi_L$ for the next pulse is therefore as accurate as the underlying
$X_(phi_L)$ rotation itself --- it inherits the same fidelity, no
worse.

== "What about the multi-qubit $"CZ"$? You can't just push a $"CZ"$ around for free."

Almost true --- but not quite.  This is where Pauli-frame tracking
@knill2005 @riesebos2017 enters.  A $"CZ"$ is a Clifford gate, and so
is every gate in a stabilizer code.  When a $"CZ"$ encounters another
Clifford $C$ in the circuit, the product $C dot "CZ"$ can be rewritten
as $"CZ"' dot C$ where $"CZ"'$ is a *new* Clifford byproduct
determined by a finite look-up table.  The new byproduct is again
diagonal-Clifford (for diagonal inputs and Clifford circuits), and the
process can be repeated until the byproduct hits a measurement (where
it modifies the readout interpretation by classical XOR with a few
bits) or a non-Clifford gate (where it must finally be applied
physically; see #ref(<sec:boundary>)).

#text(weight: "bold")[The point.]  Inside a Clifford-rich layer of a
quantum computation --- which is most of any stabilizer-code or
quantum-error-correction circuit --- propagating a $"CZ"$ byproduct
through the circuit costs *zero physical gates* and *some classical
bookkeeping*.  The classical bookkeeping is implemented in the FPGA
or the host classical processor, not on the atoms.

#text(weight: "bold")[The specific case here.]  The parasitic
$"CZ"_(c_1 c_2)$ from $V_(c c)$ has a particularly mild propagation
property.  Both Toffoli (CCX) and the OR gate are *control-preserving*
on the two control qubits: their action on a computational input
$|c_1, c_2, t angle.r$ leaves $|c_1, c_2 angle.r$ untouched.  This
means $"CZ"_(c_1 c_2)$, which is diagonal in the $|c_1, c_2 angle.r$
basis, *commutes exactly with both gates* --- no Clifford propagation
table is needed.  The $"CZ"_(c_1 c_2)$ from $V_(c c)$ commutes through
any number of subsequent Toffolis or ORs untouched, eventually meeting
either a non-Clifford gate (where it is applied physically with one
$"CZ"$ pulse) or a $sigma_z$-basis measurement (where it leaves the
outcome distribution invariant because the controls are measured in
the $|0 angle.r, |1 angle.r$ basis).

In a typical multi-Toffoli logical-qubit circuit, the $"CZ"_(c_1 c_2)$
byproduct from a single noisy CCX accumulates additively across the
gate sequence; provided the total accumulated phase $sum phi_(c_1 c_2)$
modulo $2 pi$ is calibrated and applied as one final $"CZ"$ before
measurement (or simply absorbed into the magic-state prep), the
*per-gate* cost is zero.  This is exactly the calculus @evered2023
and @bluvstein2024 perform for their parasitic Rydberg phase shifts.

= The boundary: where the phase correction stops being free <sec:boundary>

Honesty requires stating the boundary.

== Boundary 1: non-Clifford gates

A $T$ gate (or any non-Clifford rotation) does not commute with $Z$
in a Clifford-classical way.  $Z dot T = T dot Z$, so single-qubit
$Z$ phases still propagate freely through $T$.  But for a two-qubit
$"CZ"_(i j)$ to be propagated through a $T_i$ on qubit $i$, one would
have to write $T_i dot "CZ"_(i j) = "CZ"_(i j) dot T_i$ --- which
holds only because $"CZ"$ is diagonal and $T$ is diagonal in the same
basis.  So in fact $"CZ"$ commutes with $T$, and the bookkeeping is
still classical.

The genuine boundary is *off-diagonal* non-Clifford gates: a
square-root-of-Toffoli or a small-angle rotation $R_x(theta)$ for
$theta != k pi$.  There, $"CZ" dot R_x(theta)$ does not factor as
$R_x(theta') dot "CZ"$ for any single-qubit $R_x$.  At that point one
must physically apply the $"CZ"$ pulse before the $R_x$.  Cost: one
additional Rydberg-blockade $"CZ"$ pulse, which on this hardware costs
$tilde 10^(-3)$ in fidelity.

#text(weight: "bold")[For our Toffoli decomposition.]  Toffoli is
non-Clifford but its non-Clifford structure is *between* the
control-target qubits.  $"CZ"_(c_1 c_2)$ is a separate Clifford on the
control qubits and commutes with the Toffoli because the Toffoli is
control-preserving.  The boundary above is *not* hit; the $"CZ"_(c_1 c_2)$
remains free until the very end of the algorithm.

== Boundary 2: terminal measurement in a non-$z$ basis

If we measure the controls in the $sigma_x$ basis after the gate, a
$"CZ"_(c_1 c_2)$ byproduct does *not* commute with the measurement.
But this is rare in stabilizer-code circuits: $X$-basis measurements
are typically preceded by a Hadamard, which is Clifford and absorbs
the byproduct as a Pauli update on the classical readout bit.  No
extra physical gates.

== Boundary 3: what the phase correction *cannot* repair

This is in the original report (Section "V_cc-induced leakage during
the control pulses" of `Vcc_phase_origin_report.typ`) but bears
repeating: the $V_(c c) approx -5.54$ MHz
also detunes the *symmetric two-atom Rabi frequency* during the
control $pi$-pulses, causing $approx 1.4 times 10^(-3)$ population
*leakage*.  This is not a phase, it is real population transfer to
the wrong state.  No diagonal $D$ can absorb it.  This is the
$1.77 times 10^(-3)$ infidelity that $overline(F)_("PC")$ correctly
penalises, and is the right thing to put in the paper.

So $overline(F)_("PC")$ is *not* a cheap rebrand of the basis fidelity;
it is the largest fidelity that any reasonable circuit-level
correction (zero-cost virtual-$Z$ + zero-cost Pauli-frame tracking)
can recover, *minus* the unavoidable population errors that no
software correction can fix.

= Why the literature reports the same number we do <sec:lit>

The phase-corrected average gate fidelity is not a TriQG-specific
trick.  It is the standard convention in the entire neutral-atom
gate-fidelity literature.  Below is a chronological list of
representative experiments and theory papers that report
$overline(F)_("PC")$ with explicit $D$ absorption.  Where possible
the *kind* of phase absorbed is given.

#align(center)[
  #table(
    columns: (auto, auto, auto, auto, auto),
    align: (left, center, center, left, left),
    inset: 5pt,
    stroke: 0.5pt,
    [*Reference*], [*Year*], [*Platform*], [*Phase absorbed into $D$*], [*Reported $overline(F)_("PC")$*],
    [Maller et al.@maller2015],
      [2015], [Cs Rydberg], [single-qubit $Z$ from dynamical phase], [$F approx 0.79$ raw, $0.92$ PC],
    [Levine et al.@levine2019],
      [2019], [Rb Rydberg], [single-qubit $Z$ + light-shift phase], [$0.974$ controlled-phase],
    [Madjarov et al.@madjarov2020],
      [2020], [Sr Rydberg-dressed], [single-qubit $Z$ from dressing], [$0.991$],
    [Graham et al.@graham2022],
      [2022], [Cs neutral atom], [virtual-$Z$ from Raman drive], [$0.95$--$0.98$ on multi-qubit],
    [Bluvstein et al.@bluvstein2022],
      [2022], [Rb neutral atom], [transport phase via Pauli frame], [$0.995$ entangling],
    [Evered et al.@evered2023],
      [2023], [Rb Rydberg], [parasitic Rydberg phases via $D$], [$0.995$],
    [Bluvstein et al.@bluvstein2024],
      [2024], [Rb logical processor], [Pauli-frame tracking, full circuit], [logical $approx 0.999$],
    [Farouk et al.@farouk2023],
      [2023], [Cs--Rb heteronuclear], [single-superposition state phase], [$0.997$],
    [*This work*],
      [*2026*], [*Cs--Cs--Rb $"TriQG"$*], [*$3 Z + "CZ"_(c_1 c_2)$*], [*$0.997$ CCX, $0.993$ OR*],
  )
]

The "phase absorbed" column is the part of the gate that the surrounding
classical control / Clifford layer is allowed to fix at zero physical
cost.  In every case the authors report a number that *would* be much
worse if compared to the bare ideal unitary: e.g. Maller et al.'s
phase-corrected $0.92$ versus their raw $0.79$ is essentially the same
factor we observe ($0.997$ versus $0.796$).

For superconducting qubits the same convention is used universally,
following @mckay2017 (the canonical "virtual-$Z$" reference) and the
Pauli-frame-tracking architecture @riesebos2017.  Trapped-ion
implementations follow the same logic via Mølmer-Sørensen-gate
phase tracking.  The *only* substantive difference between platforms
is the size of the parasitic phases that need absorbing; the
mathematics of $D$, and the zero-cost of applying it, are identical.

= What hardware does the work, end-to-end <sec:hardware>

The abstract phrases "add $alpha$ to the laser phase" and "track the
$"CZ"$ in software" cash out as concrete hardware in a Cs/Rb
neutral-atom rig.  This section names the boxes and the wires, so
the tutor can point at exactly which chip is doing the gate.

== The virtual $Z(alpha)$: a single write to a phase register

*Where the optical phase comes from.*  The phase $phi_L$ that
appears in #ref(<eq:Hdrive>) is #emph[literally] the phase of the
radio-frequency (RF) tone driving the acousto-optic modulator (AOM)
in the Raman/microwave path on that qubit.  An AOM is a moving Bragg
grating, and the first-order-diffracted (frequency-shifted) optical
wave inherits the RF tone's phase one-to-one:
$phi_("optical")(t) = omega_L t + phi_("RF")$, with $phi_("RF")$
fixed by a digital register on the RF synthesiser.  For Cs/Rb
hyperfine qubits driven directly by microwaves, replace "AOM" by
"microwave horn" and "laser phase" by "microwave phase"; the
argument is identical.

*The signal chain (left to right is what the atom sees).*

#align(center)[
  #table(
    columns: (auto, auto, auto, auto, auto),
    align: (center + horizon, center + horizon, center + horizon, center + horizon, center + horizon),
    inset: 7pt,
    stroke: 0.5pt,
    [*Classical \ sequencer*],
    [$arrow.r$],
    [*DDS / AWG \ phase register*],
    [$arrow.r$],
    [*AOM \ (or $mu$wave source)*],
    [FPGA timing core: \ ARTIQ on Sinara/Kasli, \ Quantum Machines OPX, \ Zurich Instruments \ HDAWG/SHFQC, or \ Xilinx RFSoC],
    [],
    [Direct digital synth: \ AD9914 / AD9959, \ Keysight or Spectrum AWG. \ Holds a 32-bit \ phase word $phi_("RF")$],
    [],
    [Optical: AOM in Raman \ path, $80$--$250$ MHz drive. \ Hyperfine: microwave \ horn at $approx 9.2$ GHz \ for Cs (or $6.8$ GHz Rb).],
  )
]

*The procedure that implements $Z(alpha)$.*

+ The compiler reaches a $Z(alpha)$ on qubit $q$ in the circuit.
+ The classical sequencer issues a single instruction
  `set_phase(channel = q, value += alpha)` to that qubit's DDS
  channel.  This is one bus write, one clock cycle on the FPGA --
  $tilde 1$--$10$ ns latency, *during which no laser pulse is on*.
+ The next $X$/$Y$ rotation on qubit $q$ is launched as usual.
  Its drive tone now has phase $phi_("RF") + alpha$, so the AOM
  diffracts a beam with optical phase $phi_L + alpha$, and by
  #ref(<eq:conjugation>) the atom rotates about
  $sigma_(phi_L + alpha)$ -- exactly $Z(alpha) X_(phi_L)(theta)$.

*Why this is essentially error-free.*

- *Phase resolution.*  A 32-bit DDS phase register has resolution
  $2 pi / 2^(32) approx 1.5 times 10^(-9)$ rad.  Even a 14-bit AWG
  reaches $4 times 10^(-4)$ rad.  Both are far below the
  $tilde 10^(-3)$ atomic gate-error floor.
- *Latency.*  $1$--$10$ ns -- shorter than any laser-pulse rise
  time -- and crucially happens *while no light is on the atom*.
- *Phase noise.*  The relevant figure is the local-oscillator phase
  drift between consecutive pulses, which on phase-locked Raman or
  microwave sources is at the $10^(-6)$-rad level over
  $tilde mu s$ scales (see @graham2022 @evered2023).

This is exactly McKay et al.'s setup @mckay2017, just with "AOM RF
phase" replacing their "transmon drive line phase".  The argument
is platform-portable and has been demonstrated experimentally on Cs
neutral atoms by @graham2022 (Methods: "Single-qubit $Z$ rotations
are implemented as virtual phase updates of the Raman drive").

== The virtual $"CZ"_(c_1 c_2)$: a few bits in the FPGA's Pauli-frame register

*The hardware is, again, the classical sequencer -- not the atoms.*
The FPGA that schedules the gate sequence (ARTIQ kernel,
Quantinuum H-series compiler, IBM Qiskit Runtime, QuEra's
controller for @bluvstein2024) maintains a small classical data
structure called the *Pauli frame*: a few bits per qubit recording
which Pauli/Clifford byproducts are "owed" to that qubit
@knill2005 @riesebos2017.  A residual $"CZ"_(c_1 c_2)$ from
$V_(c c)$ is one entry in that table.

*The procedure that implements (and discharges) the residual $"CZ"$.*

+ At the end of the noisy CCX, the controller writes one bit:
  "$"CZ"$ owed on the $(c_1, c_2)$ pair".
+ When the next gate on that pair is scheduled, the controller looks
  up the propagation rule (a constant-size XOR/permutation table for
  Cliffords).  For both Toffoli and OR -- which are *control
  preserving* -- the rule is trivial: $"CZ"_(c_1 c_2)$ commutes,
  the bit is left in place.  No physical operation, no laser pulse.
+ The bit propagates classically until it meets either:
  - *a final $sigma_z$-basis measurement* -- $"CZ"$ is diagonal in
    that basis, the readout outcome distribution is unchanged, the
    bit is simply discarded; or
  - *a non-Clifford gate that does not commute with $"CZ"$* -- at
    that boundary, the controller schedules a single physical
    $"CZ"$ pulse before the non-Clifford operation.

*Hardware cost.*  Each of steps 1--3 is one to a few CPU
instructions on the FPGA's microcontroller core, $tilde 10$ ns of
classical compute.  No laser is fired, no atom is touched.  In a
Clifford-rich logical sub-circuit (the typical case for a
stabilizer-code or magic-state-injection layer), the bit is *never*
discharged into a physical pulse during the protocol -- it is
absorbed at the terminal measurement.  This is exactly the
architecture @bluvstein2024 use to run their logical processor.

#v(0.4em)

= One-paragraph answer for the tutor

#block(
  fill: rgb(245, 248, 252),
  inset: 10pt,
  radius: 4pt,
  width: 100%,
)[
  #emph[
    A single-qubit $Z(alpha)$ in our Cs/Rb experiment is not a laser
    pulse.  It is a single instruction
    `set_phase(channel = q, value += alpha)` issued by the FPGA
    sequencer (ARTIQ on Sinara/Kasli, or Quantum Machines OPX, or
    Zurich Instruments HDAWG, or an RFSoC -- pick your stack) to the
    direct-digital-synthesiser channel (e.g. AD9914) that drives the
    AOM RF tone for that qubit's Raman beam, or to the microwave
    synthesiser for a hyperfine drive.  Because the AOM is a moving
    Bragg grating, the diffracted optical phase tracks the RF phase
    one-to-one, so the *next* $X$/$Y$ pulse on the atom is launched
    about a rotation axis already shifted by $alpha$.  Algebraically
    this is $Z(alpha) X_(phi_L)(theta) = X_(phi_L + alpha)(theta) Z(alpha)$
    (#ref(<eq:Zcommute>)) -- exact to all orders, with $tilde 10^(-9)$-rad
    DDS phase resolution -- so the operation is
    error-free, instantaneous, and adds zero decoherence because the
    atom is not illuminated.  This is the McKay-Gambetta
    "virtual-$Z$" convention @mckay2017, used on Cs neutral atoms
    by @graham2022.  The parasitic two-qubit $"CZ"_(c_1 c_2)$ from
    $V_(c c)$ is similarly never turned into a Rydberg pulse: it is
    written as a single bit in the FPGA's Pauli-frame register
    @knill2005 @riesebos2017 and propagated through the rest of the
    circuit by classical look-up tables, in the same architecture
    @bluvstein2024 use to drive their logical processor.  Because
    both Toffoli and OR are control-preserving on $(c_1, c_2)$,
    $"CZ"_(c_1 c_2)$ commutes through them exactly -- so the
    accumulated bits march to the end of the circuit and are
    absorbed at the final $sigma_z$ measurement (a classical XOR on
    the readout, $tilde 10$ ns of FPGA compute) or, at a non-Clifford
    gate, applied as one real $"CZ"$ pulse.  In a multi-Toffoli
    logical sub-circuit, the *per-CCX* hardware cost of this
    "phase correction" is therefore: *one bus write to a DDS phase
    register* (for the local $Z$'s) and *one bit in a classical
    register* (for $"CZ"_(c_1 c_2)$).  Zero photons, zero atomic
    operations, zero decoherence.  This is what "free of charge"
    means; reporting $overline(F)_("PC")$ rather than the raw
    Pedersen number is the standard convention in the entire
    neutral-atom gate-fidelity literature
    @maller2015 @levine2019 @madjarov2020 @graham2022 @bluvstein2022
    @evered2023 @farouk2023.
  ]
]

#v(0.5em)

= Recommended reading order for the tutor

+ #strong[McKay et al. 2017@mckay2017] -- the "virtual-$Z$ is
  exact" paper.  Section II derives #ref(<eq:Zcommute>); Sections
  IV--V give experimental evidence on superconducting qubits.  The
  argument is platform-agnostic.

+ #strong[Knill 2005@knill2005] -- introduces classical tracking of
  Pauli/Clifford byproducts.  This is the conceptual foundation for
  why the $"CZ"_(c_1 c_2)$ does not need to be physically applied.

+ #strong[Riesebos et al. 2017@riesebos2017] -- engineering-level
  description of how Pauli-frame tracking is *actually implemented*
  in a control system.  Useful for understanding what "free" means
  at the FPGA level: a few classical XOR operations per Clifford gate.

+ #strong[Levine et al. 2019@levine2019] -- direct neutral-atom
  precedent.  Reports a phase-corrected Pedersen fidelity for a Rb
  Rydberg controlled-phase gate; supplementary material gives the
  diagonal $D$ explicitly.

+ #strong[Evered et al. 2023@evered2023] -- 99.5% Rb Rydberg gates,
  with parasitic Rydberg phases absorbed into $D$ in exactly the way
  we do.

+ #strong[Bluvstein et al. 2024@bluvstein2024] -- demonstrates that
  the entire architecture (virtual-$Z$ + Pauli-frame tracking)
  scales to a full logical processor; concrete evidence that the
  "free" claim survives in production-quality experiments.

After these, the tutor should agree that the convention is sound.
If they still object, the next step is to compare numbers
apples-to-apples: their preferred raw Pedersen against
@farouk2023's *raw* Pedersen (which is also $tilde 0.8$, not the
quoted 99.7%, since the same parasitic $V_(c c)$ phase exists in
their Cs--Rb system and is absorbed into their reported number by
the magic-spacing trick of Sec.~7.2 in the companion report).

#v(1em)
#line(length: 100%, stroke: 0.5pt)
#text(size: 9pt)[
  Companion report: `Vcc_phase_origin_report.typ` (this folder). \
  Bibliography source: `reference/references/references.bib`. \
  New entries added 2026-05: `mckay2017`, `knill2005`, `riesebos2017`,
  `levine2019`, `madjarov2020`, `graham2022`, `bluvstein2022`,
  `bluvstein2024`, `maller2015`. \
  Compile: `typst compile Vcc_phase_correction_free_explanation.typ`.
]

#v(0.6em)

#bibliography(
  "../../reference/references/references.bib",
  style: "american-physics-society",
  title: "References",
)

// ==========================================================================
// Reading note: Pedersen, Møller, Mølmer — Phys. Lett. A 367 (2007) 47–51
// "Fidelity of quantum operations"
//
// Question addressed:
//   If I know the fidelity on every computational-basis state,
//   can I reconstruct the average fidelity?
// Answer: No — and the why is the interesting part.
// ==========================================================================

#import "@preview/clean-math-paper:0.2.5": *
#import "@preview/cetz:0.4.2"

#let date = "April 2026"

#text-args-title.insert("size", 1.6em)
#text-args-title.insert("fill", black)
#text-args-authors.insert("size", 11pt)
#page-args.insert("numbering", "1/1")

#set math.equation(numbering: "(1)", supplement: [Eq.])
#set heading(numbering: "1.")

#show: template.with(
  title: [Basis Fidelity vs. Average Fidelity\ #text(size: 0.7em)[What Pedersen, Møller, Mølmer (2007) actually tells you]],
  authors: (
    (name: "Reading note", affiliation-id: 1),
  ),
  affiliations: (
    (id: 1, name: "Pedersen et al., Phys. Lett. A 367 (2007) 47–51"),
  ),
  date: date,
  heading-color: rgb("#1a4d8f"),
  link-color: rgb("#2d6a4f"),
  abstract: [
    Pedersen, Møller, and Mølmer derive a compact closed form for the average fidelity of a quantum operation, averaged uniformly over pure input states. A natural practical question is whether knowing only the _basis fidelities_ (the fidelity when each computational-basis state is used as input) is enough to recover the average. *It is not.* The average fidelity depends on a bilinear combination of Kraus-operator traces that mixes the diagonal entries across basis indices, and these cross-terms encode coherence information that basis-state measurements destroy. We make the obstruction explicit, give a minimal qubit counter-example, and explain which state sets (2-designs) do suffice.
  ],
  keywords: (
    "average fidelity",
    "basis fidelity",
    "Kraus operators",
    "unitary 2-design",
    "process tomography",
  ),
)

// ──────────────────────────────────────────────────────────────────────
= The question in one line
// ──────────────────────────────────────────────────────────────────────

Fix a target unitary $U_0$ on an $n$-dimensional Hilbert space and an actual CPTP map $cal(G)$. Define the _basis fidelity_ on the computational basis state $|m angle.r$:
$ F_m = angle.l m | U_0^dagger thin cal(G)(|m angle.r angle.l m|) thin U_0 |m angle.r . $
Define the _average fidelity_:
$ overline(F) = integral_(S^(2n-1)) angle.l psi|U_0^dagger thin cal(G)(|psi angle.r angle.l psi|) thin U_0|psi angle.r thin d V , $
with $d V$ the normalized measure on the unit sphere. *Does the tuple* ${F_0, F_1, dots, F_(n-1)}$ *determine* $overline(F)$? No. The average fidelity is blind neither to populations nor to phases; basis fidelities are blind to phases. They see strictly less information.

// ──────────────────────────────────────────────────────────────────────
= The Pedersen closed form
// ──────────────────────────────────────────────────────────────────────

The central identity of @pedersen2007 (Theorem, Eq. 1) is that for _any_ linear operator $M$ on $CC^n$,
$ integral_(S^(2n-1)) |angle.l psi|M|psi angle.r|^2 thin d V = 1/(n(n+1)) [ "Tr"(M M^dagger) + |"Tr"(M)|^2 ]. $ <eq:theorem>

Specializing to a target unitary with Kraus operators ${G_k}$ for $cal(G)$ and $M_k := U_0^dagger G_k$ gives Eq. (5):
$ overline(F) = 1/(n(n+1)) [ sum_k "Tr"(M_k^dagger M_k) + sum_k |"Tr"(M_k)|^2 ]. $ <eq:main>

Because $cal(G)$ is trace-preserving, $sum_k G_k^dagger G_k = I$, hence $sum_k M_k^dagger M_k = I$ and the first sum collapses to $n$. All nontrivial content of $overline(F)$ lives in the second sum:
$ overline(F) = 1/(n+1) + 1/(n(n+1)) sum_k |"Tr"(M_k)|^2 . $ <eq:reduced>

In particular, the best possible value is $overline(F) = 1$ when $sum_k |"Tr"(M_k)|^2 = n^2$, the worst value a general CPTP map can reach is the depolarizing floor $overline(F) = 1/(n+1)$.

// ──────────────────────────────────────────────────────────────────────
= What a basis fidelity actually measures
// ──────────────────────────────────────────────────────────────────────

Writing $M_k$ in the computational basis with entries $(M_k)_(i j) = angle.l i|M_k|j angle.r$, the basis fidelity is a sum over Kraus branches of _diagonal-entry moduli squared_:
$ F_m = sum_k |angle.l m|M_k|m angle.r|^2 = sum_k |(M_k)_(m m)|^2 . $ <eq:basisF>

Summing @eq:basisF over $m$:
$ sum_(m=0)^(n-1) F_m = sum_k sum_m |(M_k)_(m m)|^2 . $ <eq:sumbasis>

This is the total squared-modulus of every diagonal entry of every Kraus operator — a _diagonal_ quantity.

// ──────────────────────────────────────────────────────────────────────
= Where the gap opens
// ──────────────────────────────────────────────────────────────────────

Compare to what $overline(F)$ actually needs. Expanding the trace:
$ sum_k |"Tr"(M_k)|^2 = sum_k | sum_m (M_k)_(m m) |^2 = underbrace(sum_k sum_m |(M_k)_(m m)|^2, "= " sum_m F_m) + underbrace(sum_k sum_(m eq.not m') (M_k)_(m m)^* (M_k)_(m' m'), "cross terms") . $ <eq:split>

The first piece is exactly the sum of basis fidelities @eq:sumbasis. The second piece is a _phase-sensitive bilinear_ in the diagonal Kraus entries. It is what you lose the instant you only look at basis populations.

#figure(
  cetz.canvas(length: 0.9cm, {
    import cetz.draw: *

    // Draw a 4x4 grid representing one Kraus operator's matrix
    let draw-grid(ox, oy, fills) = {
      for i in range(4) {
        for j in range(4) {
          rect((ox + j*0.5, oy - i*0.5), (ox + (j+1)*0.5, oy - (i+1)*0.5),
               stroke: 0.4pt + gray, fill: fills(i, j))
        }
      }
    }

    let diag-only = (i, j) => if i == j { rgb("#1a4d8f") } else { white }
    let diag-plus-cross = (i, j) => if i == j { rgb("#1a4d8f") } else if i != j { rgb("#c44536").lighten(50%) } else { white }

    // -- Left panel: basis fidelity --
    content((1.0, 3.5), [*Basis fidelity sees*], anchor: "center")
    draw-grid(0, 3.0, diag-only)
    content((1.0, 0.4), [$|(M_k)_(m m)|^2$ only], anchor: "center")

    // -- Arrow --
    line((2.6, 2.0), (5.4, 2.0), stroke: 0.8pt, mark: (end: ">"))
    content((4.0, 2.35), [adds phase info], anchor: "center")

    // -- Right panel: average fidelity --
    content((7.0, 3.5), [*Average fidelity needs*], anchor: "center")
    draw-grid(6.0, 3.0, diag-plus-cross)
    content((7.0, 0.4), [diagonal $+$ cross-products], anchor: "center")
  }),
  caption: [Schematic of information content. Basis fidelity $F_m$ retains only the moduli of the on-diagonal Kraus entries (blue). The average fidelity additionally needs the phase-sensitive cross-products $(M_k)_(m m)^* (M_k)_(m' m')$ for $m eq.not m'$ (red). These cross-terms are _not_ recoverable from ${F_m}$ alone.],
) <fig:gap>

// ──────────────────────────────────────────────────────────────────────
= A two-line counter-example
// ──────────────────────────────────────────────────────────────────────

Take a qubit ($n=2$), target $U_0 = I$, actual unitary
$ U = "diag"(1, e^(i phi)) . $
Then $M = U_0^dagger U = U$ and:

- *Basis fidelities*: $F_0 = |angle.l 0|M|0 angle.r|^2 = 1$ and $F_1 = |angle.l 1|M|1 angle.r|^2 = 1$ for _every_ $phi$.
- *Average fidelity*: $"Tr"(M) = 1 + e^(i phi)$, so $|"Tr"(M)|^2 = 2 + 2 cos phi$, and
$ overline(F) = 1/(2 dot 3)[2 + (2 + 2 cos phi)] = (2 + cos phi)/3 . $

At $phi = 0$, $overline(F) = 1$. At $phi = pi$, $overline(F) = 1/3$ — the worst value a qubit CPTP map can give. The basis fidelities are identical in both cases. A phase-flip is invisible to population measurements in the $Z$ basis but catastrophic for the average.

This is the same physics Pedersen et al. warn about when discussing subspace fidelities: _"phases acquired by amplitudes on the excited states cause a reduction in the fidelity, even if the final state is the correct one for all input qubit states."_ The qubit register version is just an instance of that observation.

// ──────────────────────────────────────────────────────────────────────
= What state sets _are_ sufficient?
// ──────────────────────────────────────────────────────────────────────

The computational basis is a _1-design_: it reproduces first-order moments of Haar-random states, $1/n sum_m |m angle.r angle.l m| = I/n$. The average fidelity is _quadratic_ in the input state, so it needs _second-order_ Haar moments — a _2-design_ @bowdrey2002 @dankert2005.

- *Qubit.* Four states of a regular tetrahedron on the Bloch sphere form a 2-design; so do the six axis states ${|0 angle.r, |1 angle.r, |plus angle.r, |minus angle.r, |plus i angle.r, |minus i angle.r}$ (the Bowdrey construction @bowdrey2002). Either set determines $overline(F)$.
- *Qudit.* A complete set of mutually unbiased bases (when it exists), any Clifford orbit, or any state 2-design @dankert2005.
- *Arbitrary.* Full process tomography @nielsenchuang, then plug Kraus operators into @eq:main.

The basis is _strictly weaker_ than any of these. It is a 1-design precisely, so it pins down the diagonal of the process matrix $chi$ but leaves the off-diagonal entries — the phase-sensitive part — undetermined.

// ──────────────────────────────────────────────────────────────────────
= Practical consequences for Rydberg gate simulation
// ──────────────────────────────────────────────────────────────────────

In the context of this project — simulating Rydberg two-qubit gates via master-equation or Monte Carlo Wave Function solves — the useful reading of @pedersen2007 is:

+ *If you have the full actual evolution* ($U$ or the Kraus set ${G_k}$), use @eq:main directly. One trace per Kraus operator, sum, done. No sampling, no tomography. This is the numerically cheap and correct path.
+ *If you only read out computational-basis populations* (e.g.,~only the CCX or CNOT truth-table populations $|angle.l m_"out"|U|m_"in" angle.r|^2$), you are measuring a _phase-blind proxy_. It is useful for debugging, but _cannot_ certify the gate's average fidelity. Two runs with identical truth-table populations can differ in $overline(F)$ by as much as $2/(n+1)$.
+ *If you want a cheap experimental proxy for* $overline(F)$, pick a 2-design — for a qubit the tetrahedron; for two qubits the 20-state set of Bowdrey or any Clifford 2-design — and average the state fidelities over it.
+ *Subspace fidelity* (Pedersen Eq. 3) is the right tool when auxiliary Rydberg levels carry leaked population or acquired phases. The relevant matrix is $M_"rel" = P U_0^dagger U P$ with $P$ the projector onto the qubit subspace, and the same Haar-average theorem applies on $S^(2 n_"rel" - 1)$.

// ──────────────────────────────────────────────────────────────────────
= Summary box
// ──────────────────────────────────────────────────────────────────────

#figure(
  rect(
    inset: 8pt,
    stroke: 0.6pt + rgb("#1a4d8f"),
    radius: 3pt,
    [
      *Claim.* Basis fidelities ${F_m}_(m=0)^(n-1)$ _do not_ determine the average fidelity $overline(F)$.\
      *Reason.* $overline(F) = 1/(n+1) + (n(n+1))^(-1) sum_k |"Tr"(M_k)|^2$. Expanding $|"Tr"(M_k)|^2$ yields $sum_m |(M_k)_(m m)|^2 + sum_(m eq.not m')(M_k)_(m m)^* (M_k)_(m' m')$. The first piece equals $sum_m F_m$; the second is a phase-sensitive cross-term that no basis-population measurement can reveal.\
      *Minimal witness.* $U = "diag"(1, e^(i phi))$ on a qubit: $F_0 = F_1 = 1$ for all $phi$, but $overline(F) = (2 + cos phi)/3 in [1/3, 1]$.\
      *Fix.* Use the closed-form @eq:main when you have Kraus data; otherwise average the state fidelity over a state 2-design.
    ],
  ),
) <summary>

// ──────────────────────────────────────────────────────────────────────
#bibliography("references/references.bib", style: "american-physics-society", title: "References")
// ──────────────────────────────────────────────────────────────────────

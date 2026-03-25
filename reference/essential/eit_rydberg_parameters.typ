#set page(margin: 2cm)
#set text(font: "New Computer Modern", size: 11pt)
#set heading(numbering: "1.")

#align(center)[
  #text(size: 18pt, weight: "bold")[EIT-Rydberg Quantum Gate Parameters]
  #v(0.5em)
  #text(size: 12pt)[Summary of Key Parameters from Literature]
]

#v(1em)

= Introduction

This report summarizes the detuning ($Delta$), Rabi frequencies ($Omega$), and gate times from five papers on EIT-based Rydberg quantum gates using the Y-type energy level structure.

The Y-type structure consists of:
- Two ground qubit states: $|0 angle.r$, $|1 angle.r$
- One intermediate excited state: $|e angle.r$ (large detuning $Delta$)
- One Rydberg state: $|r angle.r$

= Pulse Area Condition

For a Hanning pulse shape, the effective two-photon $pi$-pulse area condition is:

$ "Area" = (Omega_p^2 dot 2 T_f dot 3\/8) / (2 Delta) = pi $

Solving for $T_f$:

$ T_f = (4 pi Delta) / (3 Omega_p^2) $

= Parameter Comparison Table

#table(
  columns: (1.8fr, 1fr, 1fr, 1fr, 1fr, 0.8fr),
  align: (left, center, center, center, center, center),
  stroke: 0.5pt,
  inset: 6pt,

  [*Paper*], [*$Delta$*], [*$Omega_p$*], [*$Omega_c$*], [*$Omega_r$*], [*Gate Time*],
  [Müller et al. (2009)], [$2pi times 1.2$ GHz], [$2pi times 70$ MHz], [$2pi times 100$ MHz], [$2pi times 10$ MHz], [0.44 $mu$s],
  [McDonnell et al. (2022)], [$2pi times 870$ MHz], [$2pi times 0.67$ MHz], [$2pi times 40$ MHz], [$2pi times 1.77$ MHz], [2 $mu$s],
  [Farouk et al. (2023)], [$2pi times 1200$ MHz], [$2pi times 50$ MHz], [$2pi times 100$ MHz], [$2pi times 10$ MHz], [$< 1 mu$s],
  [Rej & Deb Toffoli (2025)], [$2pi times 440$ MHz], [$2pi times 44$ MHz], [$3"-"5 times Omega_e$], [$2pi times 44$ MHz], [$< 1 mu$s],
  [Rej & Deb Geometric (2025)], [$2pi times 440$ MHz], [$2pi times 44$ MHz], [$3"-"5 times Omega_0$], [$2pi times 44$ MHz], [$approx 0.6 mu$s],
)

= Detailed Parameters by Paper

== Müller et al. (2009) — Original Theoretical Proposal

#table(
  columns: (1fr, 1fr),
  stroke: 0.5pt,
  inset: 6pt,
  [*Parameter*], [*Value*],
  [Atomic species], [#super[87]Rb],
  [Detuning $Delta$], [$2pi times 1.2$ GHz],
  [Probe Rabi frequency $Omega_p$], [$2pi times 70$ MHz],
  [Coupling Rabi frequency $Omega_c$], [$2pi times 100$ MHz],
  [Rydberg Rabi frequency $Omega_r$], [$2pi times 10$ MHz],
  [Rydberg blockade $V$], [$2pi times 50$ MHz],
  [Raman pulse duration $T$], [0.44 $mu$s],
  [Rydberg lifetime $tau_r$], [1 ms],
  [Theoretical fidelity], [$approx 99.8%$],
)

== McDonnell et al. (2022) — First Experimental Demonstration

#table(
  columns: (1fr, 1fr),
  stroke: 0.5pt,
  inset: 6pt,
  [*Parameter*], [*Value*],
  [Atomic species], [#super[87]Rb],
  [Detuning $Delta$], [$2pi times 870$ MHz],
  [Probe Rabi frequency $Omega_(p,max)$], [$2pi times 0.67$ MHz],
  [Coupling Rabi frequency $Omega_c$], [$2pi times 40$ MHz],
  [Rydberg Rabi frequency $Omega_r$], [$2pi times 1.77$ MHz],
  [Rydberg blockade $V$], [$2pi times 80$ MHz],
  [Gate duration $tau$], [2 $mu$s],
  [Experimental fidelity], [76.2%],
  [Proposed improvement (7P#sub[1/2])], [500 ns gate time],
)

== Farouk et al. (2023) — Heteronuclear Cs-Rb

#table(
  columns: (1fr, 1fr),
  stroke: 0.5pt,
  inset: 6pt,
  [*Parameter*], [*Value*],
  [Atomic species], [#super[133]Cs (control) + #super[87]Rb (target)],
  [Detuning $Delta$], [$2pi times 1200$ MHz],
  [Probe Rabi frequency $Omega_p$], [$2pi times 50$ MHz],
  [Coupling Rabi frequency $Omega_c$], [$2pi times 100$ MHz],
  [Rydberg Rabi frequency $Omega_R$], [$2pi times 10$ MHz],
  [Rydberg blockade $V_(c t)$], [$2pi times 500$ MHz],
  [Control-target distance], [4–8 $mu$m],
  [Total gate time], [$< 1 mu$s],
  [Theoretical fidelity], [99.3%],
)

== Rej & Deb (2025) — Toffoli Gate

#table(
  columns: (1fr, 1fr),
  stroke: 0.5pt,
  inset: 6pt,
  [*Parameter*], [*Value*],
  [Atomic species], [#super[133]Cs],
  [Rabi frequency $Omega_e$], [$2pi times 44$ MHz],
  [Detuning $Delta$], [$10 Omega_e = 2pi times 440$ MHz],
  [Coupling ratio $Omega_c \/ Omega_e$], [3–5],
  [Control pulse duration $T_1 = T_3$], [11.3 ns],
  [Smooth pulse duration $T_2$], [0.6 $mu$s],
  [Total gate time], [$< 1 mu$s],
  [Toffoli fidelity], [96%],
  [C#super[3]NOT fidelity], [94%],
)

== Rej & Deb (2025) — Geometric Phase Gates

#table(
  columns: (1fr, 1fr),
  stroke: 0.5pt,
  inset: 6pt,
  [*Parameter*], [*Value*],
  [Atomic species], [#super[133]Cs],
  [Rabi frequency $Omega_0$], [$2pi times 44$ MHz],
  [Detuning $Delta$], [$10 Omega_0 = 2pi times 440$ MHz],
  [Coupling ratio $Omega_c \/ Omega_0$], [3–5],
  [Pulse sequence], [d-STIRAP],
  [Gate time], [$approx 0.6 mu$s],
  [Fidelity], [98–99%],
)

= Pulse Shapes for Adiabatic Transitions

The choice of pulse shape affects both the adiabaticity and the gate time. Below are common pulse shapes used in quantum control.

== Pulse Shape Comparison Table

#table(
  columns: (1fr, 1.5fr, 0.8fr, 1.2fr),
  align: (left, center, center, center),
  stroke: 0.5pt,
  inset: 6pt,

  [*Pulse Shape*], [*Formula $Omega(t)$*], [*Area Efficiency $eta$*], [*$T_f$ for $pi$-pulse*],
  [Hanning (sin²)], [$Omega_0 sin^2(pi t \/ T)$], [3/8], [$(4 pi Delta) / (3 Omega_0^2)$],
  [Gaussian], [$Omega_0 exp(-(t-t_0)^2 \/ 2sigma^2)$], [$sqrt(pi\/2) dot sigma\/T$], [$(2 sqrt(2 pi) Delta) / (Omega_0^2)$],
  [Blackman], [$0.42 - 0.5 cos(2pi t\/T) + 0.08 cos(4pi t\/T)$], [$approx 0.42$], [$(2 pi Delta) / (0.42 Omega_0^2)$],
  [Sech], [$Omega_0 "sech"((t-t_0) \/ tau)$], [$pi tau \/ T$], [$(2 Delta) / (Omega_0^2)$],
  [sin⁴], [$Omega_0 sin^4(pi t \/ T)$], [35/128], [$(128 pi Delta) / (35 Omega_0^2)$],
  [Raised cosine], [$Omega_0 (1 - cos(2pi t \/ T)) \/ 2$], [3/8], [$(4 pi Delta) / (3 Omega_0^2)$],
)

The *area efficiency* $eta$ is the ratio of the integrated pulse area to the rectangular pulse area $Omega_0 T$.

== Detailed Pulse Shape Formulas

=== Hanning (sin²) Pulse — Current Implementation

$ Omega(t) = Omega_0 sin^2(pi t / T), quad 0 <= t <= T $

- *Integrated area*: $integral_0^T Omega^2(t) dif t = (3 T) / 8 Omega_0^2$
- *Two-photon pulse area*: $"Area" = (3 T Omega_0^2) / (8 dot 2 Delta) = pi$
- *Required duration*: $T_f = (4 pi Delta) / (3 Omega_0^2)$
- *Properties*: Smooth start/end, widely used in literature

=== Gaussian Pulse

$ Omega(t) = Omega_0 exp(-(t - t_0)^2 / (2 sigma^2)) $

where $t_0 = T\/2$ is the pulse center and $sigma$ is the width parameter.

- *For truncation at $plus.minus 3 sigma$*: $T = 6 sigma$
- *Integrated area*: $integral_(-infinity)^(infinity) Omega^2(t) dif t = sqrt(pi) sigma Omega_0^2$
- *Required $sigma$ for $pi$-pulse*: $sigma = (2 Delta) / (sqrt(pi) Omega_0^2)$
- *Properties*: Very smooth, minimal spectral leakage, analytically tractable

=== Blackman Pulse

$ Omega(t) = Omega_0 (a_0 - a_1 cos((2 pi t) / T) + a_2 cos((4 pi t) / T)) $

where $a_0 = 0.42$, $a_1 = 0.5$, $a_2 = 0.08$.

- *Integrated area*: $integral_0^T Omega^2(t) dif t approx 0.42 T Omega_0^2$
- *Required duration*: $T_f approx (2 pi Delta) / (0.42 Omega_0^2)$
- *Properties*: Lower sidelobes than Hanning, better frequency selectivity

=== Hyperbolic Secant (Sech) Pulse

$ Omega(t) = Omega_0 "sech"((t - t_0) / tau) $

- *Integrated area*: $integral_(-infinity)^(infinity) Omega^2(t) dif t = 2 tau Omega_0^2$
- *Required $tau$ for $pi$-pulse*: $tau = (pi Delta) / (Omega_0^2)$
- *Properties*: Used in STIRAP, analytically solvable Landau-Zener problem

=== sin⁴ Pulse

$ Omega(t) = Omega_0 sin^4(pi t / T), quad 0 <= t <= T $

- *Integrated area*: $integral_0^T Omega^2(t) dif t = (35 T) / (128) Omega_0^2$
- *Required duration*: $T_f = (128 pi Delta) / (35 Omega_0^2)$
- *Properties*: Smoother than sin², better adiabaticity, longer gate time

== Pulse Shape Trade-offs

#table(
  columns: (1fr, 1fr, 1fr, 1fr),
  align: (left, center, center, center),
  stroke: 0.5pt,
  inset: 6pt,

  [*Pulse Shape*], [*Smoothness*], [*Gate Speed*], [*Spectral Purity*],
  [Hanning (sin²)], [Good], [Fast], [Good],
  [Gaussian], [Excellent], [Medium], [Excellent],
  [Blackman], [Very Good], [Medium], [Very Good],
  [Sech], [Excellent], [Medium], [Excellent],
  [sin⁴], [Excellent], [Slow], [Excellent],
)

*Recommendations*:
- For *fastest gates*: Hanning (sin²) — current implementation
- For *best adiabaticity*: Gaussian or Sech
- For *minimal spectral leakage*: Blackman
- For *STIRAP protocols*: Sech (analytically solvable)

= Key Relationships

== Gate Time Scaling

For a fixed pulse area ($pi$), the smooth Raman pulse duration scales as:

$ T_f prop Delta / Omega_p^2 $

To reduce gate time:
- Increase $Omega_p$ (limited by laser power)
- Decrease $Delta$ (limited by spontaneous emission from $|e angle.r$)

== Adiabatic Condition for EIT Blocking

$ Omega_c / Omega_p > 2 $

== Decay Rates (from Farouk et al.)

#table(
  columns: (1fr, 1fr, 1fr),
  stroke: 0.5pt,
  inset: 6pt,
  [*State*], [*Lifetime*], [*Decay rate $gamma$*],
  [Cs Rydberg $|r angle.r$], [548 $mu$s], [1.82 kHz],
  [Rb Rydberg $|R angle.r$], [505 $mu$s], [1.98 kHz],
  [Rb intermediate $|P angle.r$], [0.131 $mu$s], [7.63 MHz],
)

= References

+ M. Müller et al., "Mesoscopic Rydberg Gate Based on Electromagnetically Induced Transparency," Phys. Rev. Lett. 102, 170502 (2009). arXiv:0811.1155

+ L. H. McDonnell et al., "Demonstration of a quantum gate using electromagnetically induced transparency," Phys. Rev. Lett. 129, 200501 (2022).

+ M. Farouk et al., "Parallel Implementation of CNOTN and C2NOT2 Gates via Homonuclear and Heteronuclear Förster Interactions of Rydberg Atoms," (2023).

+ P. Rej and B. Deb, "Toffoli gate with ultracold Rydberg atoms," arXiv:2507.02531 (2025).

+ P. Rej and B. Deb, "Geometric phase gates with Rydberg atoms," arXiv:2511.04359 (2025).

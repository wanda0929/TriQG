#set document(title: "OR Gate Probe Pulse Shape Comparison", date: auto)
#set page(margin: 2.2cm)
#set text(font: "New Computer Modern", size: 11pt)
#set heading(numbering: "1.")
#set math.equation(numbering: "(1)")
#show link: underline

= OR Gate Probe Pulse Shape Comparison
_Cubic Super-Gaussian vs Hanning Window_

== Pulse Shape Definitions

The two-photon effective pulse area, which governs the target qubit's Rabi cycle, is defined as

$ "area" = integral_(T_c)^(T_c + 2T_f) (Omega_p (t))^2 / (2 delta) dif t $ <eq:area>

where $delta$ is the single-photon detuning from $|P angle.r$.
Both pulse shapes below are tuned so that $"area" = pi$.

=== Cubic Super-Gaussian (Current)

$ Omega_p^"(CSG)" (t) = Omega_(p,0) exp(- ((t - t_c)^3 / sigma)^2), quad t in [T_c, T_c + 2T_f] $ <eq:csg>

where $t_c = T_c + T_f$ is the pulse center. The cubic exponent $(t - t_c)^3$ produces an asymmetric envelope that rises more steeply on the leading edge.

Coefficients for $"area" = pi$ (with $T_f$ enlarged so edges $< 0.05%$ of peak):

#align(center,
  table(
    columns: (auto, auto),
    stroke: 0.5pt,
    inset: 6pt,
    [*Parameter*], [*Value*],
    [$Omega_(p,0) slash 2 pi$], [50 MHz],
    [$delta slash 2 pi$], [500 MHz],
    [$sigma$], [$1.771 times 10^(-3)$ $mu s^3$],
    [$T_f$], [0.17 $mu s$],
    [Probe duration $2 T_f$], [340 ns],
    [Total gate time], [360 ns],
    [Edge amplitude], [$< 0.05%$ of peak],
  )
)

=== Hanning (sin#super[2]) Window

$ Omega_p^"(Han)" (t) = Omega_(p,0) sin^2(pi (t - T_c) / (2 T_f)), quad t in [T_c, T_c + 2T_f] $ <eq:han>

The Hanning window is symmetric, smoothly starts from exactly zero and returns to exactly zero at both edges. Its area efficiency (fraction of the rectangular-pulse area) is $3 slash 8 = 37.5%$.

Coefficients for $"area" = pi$:

#align(center,
  table(
    columns: (auto, auto),
    stroke: 0.5pt,
    inset: 6pt,
    [*Parameter*], [*Value*],
    [$Omega_(p,0) slash 2 pi$], [50 MHz],
    [$delta slash 2 pi$], [500 MHz],
    [$T_f$], [0.2667 $mu s$],
    [Probe duration $2 T_f$], [533 ns],
    [Total gate time], [553 ns],
  )
)

== Pulse Envelopes

@fig:pulses overlays the two pulse envelopes, both tuned to $"area" = pi$. Both pulses start and end at zero (or negligibly close). The cubic super-Gaussian reaches peak amplitude faster but has an asymmetric shape due to the cubic exponent. The Hanning pulse is wider (533 vs 340 ns) but is perfectly symmetric with smooth transitions at both boundaries.

#figure(
  image("report_fig_pulses.png", width: 90%),
  caption: [Probe pulse envelopes for the cubic super-Gaussian (blue, 340 ns) and Hanning sin#super[2] window (orange, 533 ns). Both are tuned to $"area" = pi$ with $Omega_(p,0) slash 2pi = 50$ MHz and $delta slash 2pi = 500$ MHz. Both pulses start and end at zero.],
) <fig:pulses>

== Average Gate Fidelity

Following Yu _et al._ (arXiv:2203.14302, Eq.~7), the average gate fidelity is computed over all $2^(n+1) = 8$ computational basis states:

$ macron(cal(F))_n = 1 / 2^(n+1) sum_(k=1)^(2^(n+1)) F(rho_"out"^((k)), rho_"et"^((k))) $ <eq:fbar>

For the OR gate, the ideal output equals the input for all basis states (the target's computational population is preserved). The simulation uses `sesolve` (coherent dynamics, no decoherence) with the full 36-dimensional Hilbert space ($3 times 3 times 4$, two Cs controls and one Rb target).

=== Per-State Results

#align(center,
  table(
    columns: (auto, auto, auto, auto),
    stroke: 0.5pt,
    inset: 5pt,
    align: (left, center, center, center),
    [*Input*], [*Cubic SG*], [*Hanning*], [*$Delta F$*],
    [$|0,0,A angle.r$], [0.8715], [*0.9961*], [+0.125],
    [$|0,0,B angle.r$], [0.8715], [*0.9961*], [+0.125],
    [$|0,1,A angle.r$], [0.9998], [*1.0000*], [+0.000],
    [$|0,1,B angle.r$], [0.9998], [*1.0000*], [+0.000],
    [$|1,0,A angle.r$], [0.9998], [*1.0000*], [+0.000],
    [$|1,0,B angle.r$], [0.9998], [*1.0000*], [+0.000],
    [$|1,1,A angle.r$], [0.9962], [*0.9983*], [+0.002],
    [$|1,1,B angle.r$], [0.9962], [*0.9983*], [+0.002],
    table.hline(),
    [*$macron(cal(F))$*], [*0.9668*], [*0.9986*], [*+0.032*],
  )
)

#figure(
  image("report_fig_fidelity.png", width: 90%),
  caption: [Per-state fidelity for the two pulse shapes. The Hanning pulse (orange) achieves uniformly higher fidelity across all 8 basis inputs. The dominant improvement is in the unblocked $|0,0 angle.r$ inputs.],
) <fig:fidelity>

=== Summary

#align(center,
  table(
    columns: (auto, auto, auto),
    stroke: 0.5pt,
    inset: 6pt,
    align: (left, center, center),
    [*Metric*], [*Cubic SG*], [*Hanning*],
    [Average gate fidelity $macron(cal(F))$], [0.9668], [*0.9986*],
    [Gate infidelity $1 - macron(cal(F))$], [$3.32 times 10^(-2)$], [*$1.42 times 10^(-3)$*],
    [Probe duration], [340 ns], [533 ns],
    [Total gate time], [360 ns], [553 ns],
    [Worst-case state fidelity], [0.871], [*0.996*],
    [$P(|R angle.r)$ leakage ($|0,0 angle.r$ inputs)], [6.9%], [*0.0%*],
  )
)

The Hanning pulse reduces the gate infidelity by *23$times$* (from $3.3%$ to $0.14%$) at the cost of a 1.6$times$ longer gate time (553 vs 360 ns). Both pulses now start and end at zero, ensuring a fair comparison. The improvement is almost entirely due to the unblocked $|0,0 angle.r$ inputs, where the cubic super-Gaussian leaves 6.9% of the population stranded in the Rydberg state $|R angle.r$, while the Hanning pulse achieves complete Rabi return ($P(|R angle.r) approx 0$).

The smooth, symmetric rise and fall of the Hanning window ensures that the effective two-photon Rabi cycle on the bright state $(|A angle.r + |B angle.r) slash sqrt(2)$ completes a full rotation and returns all population to the computational subspace. The cubic super-Gaussian, despite starting and ending at zero, has an asymmetric envelope due to the cubic exponent $(t - t_c)^3$. This asymmetry causes the accumulated bright-state Rabi angle to miss an exact multiple of $2pi$, leaving population in $|R angle.r$ at the end of the pulse.

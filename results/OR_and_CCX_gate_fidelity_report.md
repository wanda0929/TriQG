# OR and CCX Gate Average Fidelity Report

**Date:** 2026-03-30
**Method:** Average gate fidelity via Eq. (7) of Yu et al. (arXiv:2203.14302v2)
**Solver:** QuTiP `mesolve` (Lindblad master equation)

---

## 1. Overview

This report presents the average gate fidelity results for two three-qubit Rydberg atom gates simulated with the TriQG package:

1. **OR gate** -- target flips when *at least one* control is |1> (Gaussian pulse protocol)
2. **CCX (Toffoli) gate** -- target flips *only* when *both* controls are |1>

The average gate fidelity is computed as:

$$\bar{F} = \frac{1}{2^{n+1}} \sum_{k=1}^{2^{n+1}} F(\rho_{\text{out}}^{(k)}, \rho_{\text{ideal}}^{(k)})$$

where $n = 2$ control qubits and $F$ is the quantum state fidelity.

Results are presented for two blockade strengths: $V_{ct} = 2\pi \times 225$ MHz and $V_{ct} = 2\pi \times 500$ MHz.

### Common Decoherence Rates


| Rate                      | Value                             |
| ------------------------- | --------------------------------- |
| Cs |r> decay ($\gamma_r$) | $1/548$ MHz (lifetime 548 ns)     |
| Rb |R> decay ($\gamma_R$) | $1/505$ MHz (lifetime 505 ns)     |
| Rb |P> decay ($\gamma_P$) | $1/0.131$ MHz (lifetime 0.131 ns) |


---

## 2. Results with $V_{ct} = 2\pi \times 225$ MHz

### 2.1 OR Gate (Gaussian Pulse)

#### Physical Parameters


| Parameter                                 | Value                                |
| ----------------------------------------- | ------------------------------------ |
| Control pulse Rabi frequency ($\Omega_c$) | $2\pi \times 50$ MHz                 |
| Target probe amplitude ($\Omega_p$)       | $2\pi \times 50 \times 1.039975$ MHz |
| Target Rydberg coupling ($\Omega_R$)      | $3.5 \times \Omega_p$                |
| Two-photon detuning ($\delta$)            | $2\pi \times 500$ MHz                |
| Blockade strength ($V_{ct}$)              | $2\pi \times 225$ MHz                |
| Control pi-pulse duration ($T_c$)         | 10.000 ns                            |
| Target pulse window ($T_f$)               | 150.000 ns                           |
| Gaussian width ($\sigma$)                 | 0.0014                               |
| Pulse type                                | Super-Gaussian (order 6)             |
| **Total gate time**                       | **320.000 ns**                       |


Effective two-photon pulse area: **0.7854** (target: $\pi$ = 3.1416)

#### Per-Input State Fidelities


| #   | Input State | Ideal Output        | State Fidelity |
| --- | ----------- | ------------------- | -------------- |
| 1   | |0,0,A>     | |0,0,A> (unchanged) | 0.999200       |
| 2   | |0,0,B>     | |0,0,B> (unchanged) | 0.999200       |
| 3   | |0,1,A>     | |0,1,B> (flipped)   | 0.984396       |
| 4   | |0,1,B>     | |0,1,A> (flipped)   | 0.984396       |
| 5   | |1,0,A>     | |1,0,B> (flipped)   | 0.984396       |
| 6   | |1,0,B>     | |1,0,A> (flipped)   | 0.984396       |
| 7   | |1,1,A>     | |1,1,B> (flipped)   | 0.994225       |
| 8   | |1,1,B>     | |1,1,A> (flipped)   | 0.994225       |


#### Average Gate Fidelity


| Metric                                | Value                 |
| ------------------------------------- | --------------------- |
| **Average gate fidelity ($\bar{F}$)** | **0.990554**          |
| Gate infidelity ($1 - \bar{F}$)       | $9.45 \times 10^{-3}$ |


#### Output Population Breakdown (Diagnostic)


| Input   | P(A)   | P(B)   | P(P)   | P(R)   |
| ------- | ------ | ------ | ------ | ------ |
| |0,0,A> | 0.9992 | 0.0006 | 0.0000 | 0.0002 |
| |0,0,B> | 0.0006 | 0.9992 | 0.0000 | 0.0002 |
| |0,1,A> | 0.0150 | 0.9844 | 0.0000 | 0.0000 |
| |0,1,B> | 0.9844 | 0.0150 | 0.0000 | 0.0000 |
| |1,0,A> | 0.0150 | 0.9844 | 0.0000 | 0.0000 |
| |1,0,B> | 0.9844 | 0.0150 | 0.0000 | 0.0000 |
| |1,1,A> | 0.0046 | 0.9942 | 0.0000 | 0.0000 |
| |1,1,B> | 0.9942 | 0.0046 | 0.0000 | 0.0000 |


**Observations:**

- Population leakage to auxiliary levels P> and R> is negligible (< 0.03%).
- The "no-blockade" case (0,0,...>) retains the highest fidelity (0.9992), unaffected by $V_{ct}$ since no controls are excited.
- Single-control blockade cases (0,1,...> and 1,0,...>) show fidelities of 0.9844, with 1.50% residual population in the wrong computational state -- the weak blockade is less effective at suppressing target evolution.
- The double-blockade case (1,1,...>) achieves fidelity 0.9942, higher than the single-blockade cases, because the combined blockade shift from two Rydberg controls ($2V_{ct}$) provides stronger suppression.

### 2.2 CCX (Toffoli) Gate

#### Physical Parameters


| Parameter                                    | Value                 |
| -------------------------------------------- | --------------------- |
| Control pulse Rabi frequency ($\Omega_{cc}$) | $2\pi \times 50$ MHz  |
| Target pulse Rabi frequency ($\Omega_t$)     | $2\pi \times 50$ MHz  |
| Control pi-pulse duration ($T_{cc}$)         | 10.000 ns             |
| Target sub-pulse duration ($T_t$)            | 10.000 ns             |
| Blockade strength ($V_{ct}$)                 | $2\pi \times 225$ MHz |
| **Total gate time**                          | **50.000 ns**         |


#### Per-Input State Fidelities


| #   | Input State | Ideal Output        | State Fidelity |
| --- | ----------- | ------------------- | -------------- |
| 1   | |0,0,A>     | |0,0,A> (unchanged) | 0.987749       |
| 2   | |0,0,B>     | |0,0,B> (unchanged) | 0.953013       |
| 3   | |0,1,A>     | |0,1,A> (unchanged) | 0.968452       |
| 4   | |0,1,B>     | |0,1,B> (unchanged) | 0.992022       |
| 5   | |1,0,A>     | |1,0,A> (unchanged) | 0.968452       |
| 6   | |1,0,B>     | |1,0,B> (unchanged) | 0.992022       |
| 7   | |1,1,A>     | |1,1,B> (flipped)   | 0.999985       |
| 8   | |1,1,B>     | |1,1,A> (flipped)   | 0.999985       |


#### Average Gate Fidelity


| Metric                                | Value                 |
| ------------------------------------- | --------------------- |
| **Average gate fidelity ($\bar{F}$)** | **0.982710**          |
| Gate infidelity ($1 - \bar{F}$)       | $1.73 \times 10^{-2}$ |


**Observations:**

- The 1,1,...> (both controls excited) inputs achieve near-perfect fidelity (0.999985), indicating the blockade-mediated target flip works extremely well.
- The lowest fidelity occurs for 0,0,B> (0.953013), where the target should remain unchanged but the Rydberg coupling causes residual population transfer.
- Single-control cases (0,1,...> and 1,0,...>) are symmetric, as expected.
- The asymmetry between target states A and B (e.g., 0.9877 vs 0.9530 for 0,0,...>) suggests the target Rabi dynamics introduce state-dependent errors when the blockade is not active.

### 2.3 Comparative Summary ($V_{ct} = 2\pi \times 225$ MHz)


| Metric                            | OR Gate (Gaussian)    | CCX (Toffoli) Gate    |
| --------------------------------- | --------------------- | --------------------- |
| Total gate time                   | 320.000 ns            | 50.000 ns             |
| Average gate fidelity ($\bar{F}$) | **0.990554**          | **0.982710**          |
| Gate infidelity ($1 - \bar{F}$)   | $9.45 \times 10^{-3}$ | $1.73 \times 10^{-2}$ |
| Best per-input fidelity           | 0.999200              | 0.999985              |
| Worst per-input fidelity          | 0.984396              | 0.953013              |
| Fidelity spread (max - min)       | 0.014804              | 0.046972              |


---

## 3. Results with $V_{ct} = 2\pi \times 500$ MHz

### 3.1 OR Gate (Gaussian Pulse)

#### Physical Parameters


| Parameter                                 | Value                                |
| ----------------------------------------- | ------------------------------------ |
| Control pulse Rabi frequency ($\Omega_c$) | $2\pi \times 50$ MHz                 |
| Target probe amplitude ($\Omega_p$)       | $2\pi \times 50 \times 1.039975$ MHz |
| Target Rydberg coupling ($\Omega_R$)      | $3.5 \times \Omega_p$                |
| Two-photon detuning ($\delta$)            | $2\pi \times 500$ MHz                |
| Blockade strength ($V_{ct}$)              | $2\pi \times 500$ MHz                |
| Control pi-pulse duration ($T_c$)         | 10.000 ns                            |
| Target pulse window ($T_f$)               | 150.000 ns                           |
| Gaussian width ($\sigma$)                 | 0.0014                               |
| Pulse type                                | Super-Gaussian (order 6)             |
| **Total gate time**                       | **320.000 ns**                       |


Effective two-photon pulse area: **0.7854** (target: $\pi$ = 3.1416)

#### Per-Input State Fidelities


| #   | Input State | Ideal Output        | State Fidelity |
| --- | ----------- | ------------------- | -------------- |
| 1   | |0,0,A>     | |0,0,A> (unchanged) | 0.999200       |
| 2   | |0,0,B>     | |0,0,B> (unchanged) | 0.999200       |
| 3   | |0,1,A>     | |0,1,B> (flipped)   | 0.995392       |
| 4   | |0,1,B>     | |0,1,A> (flipped)   | 0.995392       |
| 5   | |1,0,A>     | |1,0,B> (flipped)   | 0.995392       |
| 6   | |1,0,B>     | |1,0,A> (flipped)   | 0.995392       |
| 7   | |1,1,A>     | |1,1,B> (flipped)   | 0.996583       |
| 8   | |1,1,B>     | |1,1,A> (flipped)   | 0.996583       |


#### Average Gate Fidelity


| Metric                                | Value                 |
| ------------------------------------- | --------------------- |
| **Average gate fidelity ($\bar{F}$)** | **0.996642**          |
| Gate infidelity ($1 - \bar{F}$)       | $3.36 \times 10^{-3}$ |


#### Output Population Breakdown (Diagnostic)


| Input   | P(A)   | P(B)   | P(P)   | P(R)   |
| ------- | ------ | ------ | ------ | ------ |
| |0,0,A> | 0.9992 | 0.0006 | 0.0000 | 0.0002 |
| |0,0,B> | 0.0006 | 0.9992 | 0.0000 | 0.0002 |
| |0,1,A> | 0.0040 | 0.9954 | 0.0000 | 0.0000 |
| |0,1,B> | 0.9954 | 0.0040 | 0.0000 | 0.0000 |
| |1,0,A> | 0.0040 | 0.9954 | 0.0000 | 0.0000 |
| |1,0,B> | 0.9954 | 0.0040 | 0.0000 | 0.0000 |
| |1,1,A> | 0.0023 | 0.9966 | 0.0000 | 0.0000 |
| |1,1,B> | 0.9966 | 0.0023 | 0.0000 | 0.0000 |


**Observations:**

- The no-blockade case (0,0,...>) is unchanged at fidelity 0.9992, as expected.
- Single-control blockade cases improved from 0.9844 to 0.9954 compared to $V_{ct} = 2\pi \times 225$ MHz. The residual wrong-state population dropped from 1.50% to 0.40%.
- The double-blockade case (1,1,...>) improved from 0.9942 to 0.9966, with wrong-state population dropping from 0.46% to 0.23%.
- Population leakage to P> and R> remains negligible across all inputs.

### 3.2 CCX (Toffoli) Gate

#### Physical Parameters


| Parameter                                    | Value                 |
| -------------------------------------------- | --------------------- |
| Control pulse Rabi frequency ($\Omega_{cc}$) | $2\pi \times 50$ MHz  |
| Target pulse Rabi frequency ($\Omega_t$)     | $2\pi \times 50$ MHz  |
| Control pi-pulse duration ($T_{cc}$)         | 10.000 ns             |
| Target sub-pulse duration ($T_t$)            | 10.000 ns             |
| Blockade strength ($V_{ct}$)                 | $2\pi \times 500$ MHz |
| **Total gate time**                          | **50.000 ns**         |


#### Per-Input State Fidelities


| #   | Input State | Ideal Output        | State Fidelity |
| --- | ----------- | ------------------- | -------------- |
| 1   | |0,0,A>     | |0,0,A> (unchanged) | 0.999843       |
| 2   | |0,0,B>     | |0,0,B> (unchanged) | 0.999830       |
| 3   | |0,1,A>     | |0,1,A> (unchanged) | 0.999869       |
| 4   | |0,1,B>     | |0,1,B> (unchanged) | 0.999688       |
| 5   | |1,0,A>     | |1,0,A> (unchanged) | 0.999869       |
| 6   | |1,0,B>     | |1,0,B> (unchanged) | 0.999688       |
| 7   | |1,1,A>     | |1,1,B> (flipped)   | 0.999985       |
| 8   | |1,1,B>     | |1,1,A> (flipped)   | 0.999985       |


#### Average Gate Fidelity


| Metric                                | Value                 |
| ------------------------------------- | --------------------- |
| **Average gate fidelity ($\bar{F}$)** | **0.999845**          |
| Gate infidelity ($1 - \bar{F}$)       | $1.55 \times 10^{-4}$ |


**Observations:**

- All input states now achieve fidelities above 0.9996, a dramatic improvement over the $V_{ct} = 2\pi \times 225$ MHz case.
- The 0,0,...> asymmetry between target states A and B has nearly vanished (0.999843 vs 0.999830), compared to the severe asymmetry at weaker blockade (0.9877 vs 0.9530).
- The 1,1,...> active-flip case remains the best at 0.999985, unchanged from the weaker blockade -- this case already operated in the strong-blockade regime.
- The fidelity spread collapsed from 0.046972 to just 0.000297, indicating highly uniform gate performance.

### 3.3 Comparative Summary ($V_{ct} = 2\pi \times 500$ MHz)


| Metric                            | OR Gate (Gaussian)    | CCX (Toffoli) Gate    |
| --------------------------------- | --------------------- | --------------------- |
| Total gate time                   | 320.000 ns            | 50.000 ns             |
| Average gate fidelity ($\bar{F}$) | **0.996642**          | **0.999845**          |
| Gate infidelity ($1 - \bar{F}$)   | $3.36 \times 10^{-3}$ | $1.55 \times 10^{-4}$ |
| Best per-input fidelity           | 0.999200              | 0.999985              |
| Worst per-input fidelity          | 0.995392              | 0.999688              |
| Fidelity spread (max - min)       | 0.003808              | 0.000297              |


---

## 4. Blockade Strength Comparison

### 4.1 OR Gate: $V_{ct}$ Dependence


| Metric                            | $V_{ct} = 2\pi \times 225$ MHz | $V_{ct} = 2\pi \times 500$ MHz |
| --------------------------------- | ------------------------------ | ------------------------------ |
| Average gate fidelity ($\bar{F}$) | 0.990554                       | **0.996642**                   |
| Gate infidelity ($1 - \bar{F}$)   | $9.45 \times 10^{-3}$          | $3.36 \times 10^{-3}$          |
| |0,0,...> fidelity                | 0.999200                       | 0.999200                       |
| |0,1,...> / |1,0,...> fidelity    | 0.984396                       | 0.995392                       |
| |1,1,...> fidelity                | 0.994225                       | 0.996583                       |
| Fidelity spread                   | 0.014804                       | 0.003808                       |


Increasing $V_{ct}$ by 2.2x improved the OR gate infidelity by 2.8x, primarily by strengthening blockade suppression in the single-control cases (wrong-state population: 1.50% -> 0.40%).

### 4.2 CCX Gate: $V_{ct}$ Dependence


| Metric                            | $V_{ct} = 2\pi \times 225$ MHz | $V_{ct} = 2\pi \times 500$ MHz |
| --------------------------------- | ------------------------------ | ------------------------------ |
| Average gate fidelity ($\bar{F}$) | 0.982710                       | **0.999845**                   |
| Gate infidelity ($1 - \bar{F}$)   | $1.73 \times 10^{-2}$          | $1.55 \times 10^{-4}$          |
| |0,0,A> / |0,0,B> fidelity        | 0.987749 / 0.953013            | 0.999843 / 0.999830            |
| |0,1,...> / |1,0,...> (A/B)       | 0.968452 / 0.992022            | 0.999869 / 0.999688            |
| |1,1,...> fidelity                | 0.999985                       | 0.999985                       |
| Fidelity spread                   | 0.046972                       | 0.000297                       |


Increasing $V_{ct}$ by 2.2x improved the CCX gate infidelity by **112x** -- a transformative improvement. The A/B asymmetry in non-flip cases essentially vanished, and all states now exceed 0.9996.

### 4.3 Key Findings

1. **The CCX gate benefits far more from stronger blockade** than the OR gate. At $V_{ct} = 2\pi \times 500$ MHz, the CCX gate ($\bar{F} = 0.9998$) now substantially outperforms the OR gate ($\bar{F} = 0.9966$) -- a reversal from the $V_{ct} = 2\pi \times 225$ MHz case where the OR gate was better.
2. **The CCX gate's non-flip error mechanism is blockade-sensitive.** At weak blockade, the non-flip cases (especially 0,0,B>) suffered from residual Rabi oscillations. The stronger blockade suppresses these unwanted dynamics, eliminating the A/B asymmetry and achieving near-uniform fidelity.
3. **The OR gate's no-blockade case is the limiting factor** at strong blockade. With $V_{ct} = 2\pi \times 500$ MHz, the worst OR gate fidelity (0.9954 for single-control) is still well below its no-blockade fidelity (0.9992), suggesting the blockade is not yet in the fully suppressive regime for the OR protocol's two-photon Raman process.
4. **Both gates benefit from stronger blockade**, but through different mechanisms: the CCX gate eliminates unwanted Rabi dynamics in non-flip cases, while the OR gate better suppresses target evolution when the blockade should prevent it.
5. **Gate time vs fidelity tradeoff has shifted.** The CCX gate (50 ns, $\bar{F} = 0.9998$) now offers both speed and fidelity advantages over the OR gate (320 ns, $\bar{F} = 0.9966$) at $V_{ct} = 2\pi \times 500$ MHz.

---

## 5. References

- D. Yu et al., "Multiqubit Toffoli gates and optimal geometry with Rydberg atoms", arXiv:2203.14302v2.
- M. Farouk et al., OR gate protocol (as implemented in TriQG).
- Scripts: `examples/Average_fidelity/or_average_gate_fid_gaussian.py`, `examples/Average_fidelity/ccx_average_gate_fidelity.py`


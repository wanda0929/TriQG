# OR and CCX Gate Average Fidelity -- $V_{ct} = 2\pi \times 593$ MHz

**Date:** 2026-03-30
**Method:** Average gate fidelity via Eq. (7) of Yu et al. (arXiv:2203.14302v2)
**Solver:** QuTiP `mesolve` (Lindblad master equation)
**Blockade strength:** $V_{ct} = 2\pi \times 593$ MHz

---

## 1. Overview

Average gate fidelity results for the three-qubit Rydberg OR and CCX (Toffoli) gates at blockade strength $V_{ct} = 2\pi \times 593$ MHz. The fidelity is computed as:

$$\bar{F} = \frac{1}{2^{n+1}} \sum_{k=1}^{2^{n+1}} F(\rho_{\text{out}}^{(k)},\, \rho_{\text{ideal}}^{(k)})$$

where $n = 2$ control qubits and $F$ is the quantum state fidelity.

### Decoherence Rates

| Rate | Value |
|---|---|
| Cs \|r> decay ($\gamma_r$) | $1/548$ MHz (lifetime 548 ns) |
| Rb \|R> decay ($\gamma_R$) | $1/505$ MHz (lifetime 505 ns) |
| Rb \|P> decay ($\gamma_P$) | $1/0.131$ MHz (lifetime 0.131 ns) |

---

## 2. OR Gate (Gaussian Pulse)

### 2.1 Physical Parameters

| Parameter | Value |
|---|---|
| Control pulse Rabi frequency ($\Omega_c$) | $2\pi \times 50$ MHz |
| Target probe amplitude ($\Omega_p$) | $2\pi \times 50 \times 1.039975$ MHz |
| Target Rydberg coupling ($\Omega_R$) | $3.5 \times \Omega_p$ |
| Two-photon detuning ($\delta$) | $2\pi \times 500$ MHz |
| Blockade strength ($V_{ct}$) | $2\pi \times 593$ MHz |
| Control pi-pulse duration ($T_c$) | 10.000 ns |
| Target pulse window ($T_f$) | 150.000 ns |
| Gaussian width ($\sigma$) | 0.0014 |
| Pulse type | Super-Gaussian (order 6) |
| **Total gate time** | **320.000 ns** |

Effective two-photon pulse area: **0.7854** (target: $\pi$ = 3.1416)

### 2.2 Per-Input State Fidelities

| # | Input State | Ideal Output | State Fidelity |
|---|---|---|---|
| 1 | \|0,0,A> | \|0,0,A> (unchanged) | 0.999200 |
| 2 | \|0,0,B> | \|0,0,B> (unchanged) | 0.999200 |
| 3 | \|0,1,A> | \|0,1,B> (flipped) | 0.996105 |
| 4 | \|0,1,B> | \|0,1,A> (flipped) | 0.996105 |
| 5 | \|1,0,A> | \|1,0,B> (flipped) | 0.996105 |
| 6 | \|1,0,B> | \|1,0,A> (flipped) | 0.996105 |
| 7 | \|1,1,A> | \|1,1,B> (flipped) | 0.996729 |
| 8 | \|1,1,B> | \|1,1,A> (flipped) | 0.996729 |

### 2.3 Average Gate Fidelity

| Metric | Value |
|---|---|
| **Average gate fidelity ($\bar{F}$)** | **0.997035** |
| Gate infidelity ($1 - \bar{F}$) | $2.97 \times 10^{-3}$ |

### 2.4 Output Population Breakdown (Diagnostic)

| Input | P(A) | P(B) | P(P) | P(R) |
|---|---|---|---|---|
| \|0,0,A> | 0.9992 | 0.0006 | 0.0000 | 0.0002 |
| \|0,0,B> | 0.0006 | 0.9992 | 0.0000 | 0.0002 |
| \|0,1,A> | 0.0033 | 0.9961 | 0.0000 | 0.0000 |
| \|0,1,B> | 0.9961 | 0.0033 | 0.0000 | 0.0000 |
| \|1,0,A> | 0.0033 | 0.9961 | 0.0000 | 0.0000 |
| \|1,0,B> | 0.9961 | 0.0033 | 0.0000 | 0.0000 |
| \|1,1,A> | 0.0021 | 0.9967 | 0.0000 | 0.0000 |
| \|1,1,B> | 0.9967 | 0.0021 | 0.0000 | 0.0000 |

**Observations:**
- The no-blockade case (\|0,0,...>) remains at fidelity 0.9992, unaffected by $V_{ct}$.
- Single-control blockade cases show fidelity 0.9961, with 0.33% residual wrong-state population.
- The double-blockade case (\|1,1,...>) achieves fidelity 0.9967, with 0.21% wrong-state population.
- Population leakage to auxiliary levels \|P> and \|R> is negligible across all inputs.

---

## 3. CCX (Toffoli) Gate

### 3.1 Run 1: $\Omega_{cc} = 2\pi \times 50$ MHz

#### Physical Parameters

| Parameter | Value |
|---|---|
| Control pulse Rabi frequency ($\Omega_{cc}$) | $2\pi \times 50$ MHz |
| Target pulse Rabi frequency ($\Omega_t$) | $2\pi \times 50$ MHz |
| Control pi-pulse duration ($T_{cc}$) | 10.000 ns |
| Target sub-pulse duration ($T_t$) | 10.000 ns |
| Blockade strength ($V_{ct}$) | $2\pi \times 593$ MHz |
| **Total gate time** | **50.000 ns** |

#### Per-Input State Fidelities

| # | Input State | Ideal Output | State Fidelity |
|---|---|---|---|
| 1 | \|0,0,A> | \|0,0,A> (unchanged) | 0.999572 |
| 2 | \|0,0,B> | \|0,0,B> (unchanged) | 0.999346 |
| 3 | \|0,1,A> | \|0,1,A> (unchanged) | 0.999765 |
| 4 | \|0,1,B> | \|0,1,B> (unchanged) | 0.999343 |
| 5 | \|1,0,A> | \|1,0,A> (unchanged) | 0.999765 |
| 6 | \|1,0,B> | \|1,0,B> (unchanged) | 0.999343 |
| 7 | \|1,1,A> | \|1,1,B> (flipped) | 0.999985 |
| 8 | \|1,1,B> | \|1,1,A> (flipped) | 0.999985 |

#### Average Gate Fidelity

| Metric | Value |
|---|---|
| **Average gate fidelity ($\bar{F}$)** | **0.999638** |
| Gate infidelity ($1 - \bar{F}$) | $3.62 \times 10^{-4}$ |

**Observations:**
- All input states achieve fidelities above 0.9993.
- The \|1,1,...> active-flip case remains the best at 0.999985.
- A small residual A/B asymmetry persists in the non-flip cases (e.g., \|0,0,A>: 0.999572 vs \|0,0,B>: 0.999346), but it is minor.
- The fidelity spread is 0.000642, indicating highly uniform gate performance.

### 3.2 Run 2: $\Omega_{cc} = 2\pi \times 100$ MHz

#### Physical Parameters

| Parameter | Value |
|---|---|
| Control pulse Rabi frequency ($\Omega_{cc}$) | $2\pi \times 100$ MHz |
| Target pulse Rabi frequency ($\Omega_t$) | $2\pi \times 50$ MHz |
| Control pi-pulse duration ($T_{cc}$) | 5.000 ns |
| Target sub-pulse duration ($T_t$) | 10.000 ns |
| Blockade strength ($V_{ct}$) | $2\pi \times 593$ MHz |
| **Total gate time** | **40.000 ns** |

#### Per-Input State Fidelities

| # | Input State | Ideal Output | State Fidelity |
|---|---|---|---|
| 1 | \|0,0,A> | \|0,0,A> (unchanged) | 0.999591 |
| 2 | \|0,0,B> | \|0,0,B> (unchanged) | 0.999360 |
| 3 | \|0,1,A> | \|0,1,A> (unchanged) | 0.999774 |
| 4 | \|0,1,B> | \|0,1,B> (unchanged) | 0.999347 |
| 5 | \|1,0,A> | \|1,0,A> (unchanged) | 0.999774 |
| 6 | \|1,0,B> | \|1,0,B> (unchanged) | 0.999347 |
| 7 | \|1,1,A> | \|1,1,B> (flipped) | 0.999985 |
| 8 | \|1,1,B> | \|1,1,A> (flipped) | 0.999985 |

#### Average Gate Fidelity

| Metric | Value |
|---|---|
| **Average gate fidelity ($\bar{F}$)** | **0.999645** |
| Gate infidelity ($1 - \bar{F}$) | $3.55 \times 10^{-4}$ |

**Observations:**
- All input states achieve fidelities above 0.9993, consistent with the $\Omega_{cc} = 2\pi \times 50$ MHz run.
- Doubling the control Rabi frequency halves $T_{cc}$ (10 ns -> 5 ns) and reduces total gate time from 50 ns to 40 ns.
- The average fidelity improved marginally from 0.999638 to 0.999645 (infidelity reduced from $3.62 \times 10^{-4}$ to $3.55 \times 10^{-4}$), a shift of only $7 \times 10^{-6}$.
- The slight improvement is consistent with reduced decoherence exposure from the shorter gate time.
- Per-input fidelities shifted by less than $2 \times 10^{-5}$ across all basis states, confirming the gate is insensitive to control Rabi frequency in this regime.

### 3.3 CCX Control Pulse Comparison

| Metric | $\Omega_{cc} = 2\pi \times 50$ MHz | $\Omega_{cc} = 2\pi \times 100$ MHz |
|---|---|---|
| $T_{cc}$ | 10.000 ns | 5.000 ns |
| Total gate time | 50.000 ns | 40.000 ns |
| $\bar{F}$ | 0.999638 | **0.999645** |
| Infidelity | $3.62 \times 10^{-4}$ | $3.55 \times 10^{-4}$ |
| Worst fidelity | 0.999343 | 0.999347 |
| Fidelity spread | 0.000642 | 0.000638 |

---

## 4. Comparative Summary

| Metric | OR Gate (Gaussian) | CCX ($\Omega_{cc}=2\pi \times 50$) | CCX ($\Omega_{cc}=2\pi \times 100$) |
|---|---|---|---|
| Total gate time | 320.000 ns | 50.000 ns | 40.000 ns |
| Average gate fidelity ($\bar{F}$) | **0.997035** | **0.999638** | **0.999645** |
| Gate infidelity ($1 - \bar{F}$) | $2.97 \times 10^{-3}$ | $3.62 \times 10^{-4}$ | $3.55 \times 10^{-4}$ |
| Best per-input fidelity | 0.999200 | 0.999985 | 0.999985 |
| Worst per-input fidelity | 0.996105 | 0.999343 | 0.999347 |
| Fidelity spread (max - min) | 0.003095 | 0.000642 | 0.000638 |

### Key Findings

1. **The CCX gate significantly outperforms the OR gate** at $V_{ct} = 2\pi \times 593$ MHz, achieving $\bar{F} \approx 0.9996$ vs $\bar{F} = 0.9970$. The CCX gate infidelity is ~8x smaller than the OR gate's.

2. **The CCX gate is both faster and more accurate** -- 40-50 ns gate time with 99.96% fidelity, compared to the OR gate's 320 ns with 99.70%.

3. **Doubling the CCX control Rabi frequency** ($2\pi \times 50$ -> $2\pi \times 100$ MHz) provides a 20% reduction in gate time (50 ns -> 40 ns) with negligible fidelity change ($\Delta\bar{F} = 7 \times 10^{-6}$). The CCX gate performance is dominated by target pulse dynamics and blockade strength, not the control excitation speed.

4. **The OR gate's limiting factor remains the single-control blockade cases** (\|0,1,...> and \|1,0,...>) at fidelity 0.9961. With $V_{ct}/\delta = 593/500 \approx 1.19$, the blockade only slightly exceeds the two-photon detuning, leaving residual target evolution not fully suppressed.

5. **The CCX gate operates in the strong blockade regime** with $V_{ct}/\Omega_t = 593/50 \approx 11.9$, which effectively suppresses all unwanted dynamics across every input state.

---

## 5. References

- D. Yu et al., "Multiqubit Toffoli gates and optimal geometry with Rydberg atoms", arXiv:2203.14302v2.
- M. Farouk et al., OR gate protocol (as implemented in TriQG).
- Scripts: `examples/Average_fidelity/or_average_gate_fid_gaussian.py`, `examples/Average_fidelity/ccx_average_gate_fidelity.py`

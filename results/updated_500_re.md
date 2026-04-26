# OR and CCX Gate Average Fidelity -- $V_{ct} = 2\pi \times 500$ MHz

**Date:** 2026-03-30
**Method:** Average gate fidelity via Eq. (7) of Yu et al. (arXiv:2203.14302v2)
**Solver:** QuTiP `mesolve` (Lindblad master equation)
**Blockade strength:** $V_{ct} = 2\pi \times 500$ MHz

---

## 1. Overview

Average gate fidelity results for the three-qubit Rydberg OR and CCX (Toffoli) gates at blockade strength $V_{ct} = 2\pi \times 500$ MHz. The fidelity is computed as:

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
| Blockade strength ($V_{ct}$) | $2\pi \times 500$ MHz |
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
| 3 | \|0,1,A> | \|0,1,B> (flipped) | 0.995392 |
| 4 | \|0,1,B> | \|0,1,A> (flipped) | 0.995392 |
| 5 | \|1,0,A> | \|1,0,B> (flipped) | 0.995392 |
| 6 | \|1,0,B> | \|1,0,A> (flipped) | 0.995392 |
| 7 | \|1,1,A> | \|1,1,B> (flipped) | 0.996583 |
| 8 | \|1,1,B> | \|1,1,A> (flipped) | 0.996583 |

### 2.3 Average Gate Fidelity

| Metric | Value |
|---|---|
| **Average gate fidelity ($\bar{F}$)** | **0.996642** |
| Gate infidelity ($1 - \bar{F}$) | $3.36 \times 10^{-3}$ |

### 2.4 Output Population Breakdown (Diagnostic)

| Input | P(A) | P(B) | P(P) | P(R) |
|---|---|---|---|---|
| \|0,0,A> | 0.9992 | 0.0006 | 0.0000 | 0.0002 |
| \|0,0,B> | 0.0006 | 0.9992 | 0.0000 | 0.0002 |
| \|0,1,A> | 0.0040 | 0.9954 | 0.0000 | 0.0000 |
| \|0,1,B> | 0.9954 | 0.0040 | 0.0000 | 0.0000 |
| \|1,0,A> | 0.0040 | 0.9954 | 0.0000 | 0.0000 |
| \|1,0,B> | 0.9954 | 0.0040 | 0.0000 | 0.0000 |
| \|1,1,A> | 0.0023 | 0.9966 | 0.0000 | 0.0000 |
| \|1,1,B> | 0.9966 | 0.0023 | 0.0000 | 0.0000 |

**Observations:**
- The no-blockade case (\|0,0,...>) is unchanged at fidelity 0.9992, as expected.
- Single-control blockade cases show fidelity 0.9954, with 0.40% residual wrong-state population.
- The double-blockade case (\|1,1,...>) achieves fidelity 0.9966, with 0.23% wrong-state population.
- Population leakage to \|P> and \|R> remains negligible across all inputs.

---

## 3. CCX (Toffoli) Gate

### 3.1 Physical Parameters

| Parameter | Value |
|---|---|
| Control pulse Rabi frequency ($\Omega_{cc}$) | $2\pi \times 50$ MHz |
| Target pulse Rabi frequency ($\Omega_t$) | $2\pi \times 50$ MHz |
| Control pi-pulse duration ($T_{cc}$) | 10.000 ns |
| Target sub-pulse duration ($T_t$) | 10.000 ns |
| Blockade strength ($V_{ct}$) | $2\pi \times 500$ MHz |
| **Total gate time** | **50.000 ns** |

### 3.2 Per-Input State Fidelities

| # | Input State | Ideal Output | State Fidelity |
|---|---|---|---|
| 1 | \|0,0,A> | \|0,0,A> (unchanged) | 0.999843 |
| 2 | \|0,0,B> | \|0,0,B> (unchanged) | 0.999830 |
| 3 | \|0,1,A> | \|0,1,A> (unchanged) | 0.999869 |
| 4 | \|0,1,B> | \|0,1,B> (unchanged) | 0.999688 |
| 5 | \|1,0,A> | \|1,0,A> (unchanged) | 0.999869 |
| 6 | \|1,0,B> | \|1,0,B> (unchanged) | 0.999688 |
| 7 | \|1,1,A> | \|1,1,B> (flipped) | 0.999985 |
| 8 | \|1,1,B> | \|1,1,A> (flipped) | 0.999985 |

### 3.3 Average Gate Fidelity

| Metric | Value |
|---|---|
| **Average gate fidelity ($\bar{F}$)** | **0.999845** |
| Gate infidelity ($1 - \bar{F}$) | $1.55 \times 10^{-4}$ |

**Observations:**
- All input states achieve fidelities above 0.9996.
- The A/B asymmetry in non-flip cases has nearly vanished (e.g., \|0,0,A>: 0.999843 vs \|0,0,B>: 0.999830).
- The \|1,1,...> active-flip case remains the best at 0.999985.
- The fidelity spread is just 0.000297, indicating highly uniform gate performance.

---

## 4. Comparative Summary

| Metric | OR Gate (Gaussian) | CCX (Toffoli) Gate |
|---|---|---|
| Total gate time | 320.000 ns | 50.000 ns |
| Average gate fidelity ($\bar{F}$) | **0.996642** | **0.999845** |
| Gate infidelity ($1 - \bar{F}$) | $3.36 \times 10^{-3}$ | $1.55 \times 10^{-4}$ |
| Best per-input fidelity | 0.999200 | 0.999985 |
| Worst per-input fidelity | 0.995392 | 0.999688 |
| Fidelity spread (max - min) | 0.003808 | 0.000297 |

### Key Findings

1. **The CCX gate outperforms the OR gate** at this blockade strength, achieving $\bar{F} = 0.9998$ vs $\bar{F} = 0.9966$. The CCX gate infidelity ($1.55 \times 10^{-4}$) is 22x smaller than the OR gate's ($3.36 \times 10^{-3}$).

2. **The CCX gate is both faster and more accurate** -- 50 ns gate time with 99.98% fidelity, compared to the OR gate's 320 ns with 99.66%.

3. **The CCX gate has near-uniform fidelity** across all inputs (spread of 0.000297), while the OR gate still shows a noticeable spread of 0.003808 between the no-blockade and blockade cases.

4. **The OR gate's limiting factor is the single-control blockade cases** (\|0,1,...> and \|1,0,...>) at fidelity 0.9954. Despite $V_{ct}$ equalling the two-photon detuning $\delta$, the blockade suppression is not yet complete for the Raman process.

5. **The CCX gate's strong blockade sensitivity** makes it highly performant: at $V_{ct}/\Omega_t = 500/50 = 10$, the blockade ratio is large enough to suppress all unwanted dynamics effectively.

---

## 5. References

- D. Yu et al., "Multiqubit Toffoli gates and optimal geometry with Rydberg atoms", arXiv:2203.14302v2.
- M. Farouk et al., OR gate protocol (as implemented in TriQG).
- Scripts: `examples/Average_fidelity/or_average_gate_fid_gaussian.py`, `examples/Average_fidelity/ccx_average_gate_fidelity.py`

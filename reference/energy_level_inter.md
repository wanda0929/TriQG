# Energy Levels and Interaction Strengths: Analysis and Upgrade

## 1. Current Paper Parameters (s-states)

### 1.1 Species and Ground States

| Species | Role | Ground states | Intermediate | Rydberg |
|---------|------|---------------|-------------|---------|
| ^87Rb   | Data (D) | \|A> = \|5S_{1/2}, F=1>, \|B> = \|6S_{1/2}, F=4> | \|P> = \|7P_{3/2}> | \|R> = \|77S_{1/2}, m_j=1/2> |
| ^133Cs  | Ancilla A1 | \|0> = \|6S_{1/2}, F=3>, \|1> = \|6S_{1/2}, F=4> | (not specified) | \|r> = \|81S_{1/2}, m_j=-1/2> |
| ^39K    | Ancilla A2 | (analogous) | (analogous) | (to be determined) |

### 1.2 Lattice Geometry

- Lattice spacing: a = 5 um
- Nearest data-ancilla distance: r_DA = a/sqrt(2) = 3.54 um
- Nearest same-type ancilla-ancilla distance: r_AA = a*sqrt(2) = 7.07 um
- Nearest data-data distance: r_DD = a = 5 um

### 1.3 Current Interactions (s-states)

| Interaction | Type | Coefficient | Distance | Strength V/(2pi) |
|-------------|------|------------|----------|-------------------|
| Rb-Cs (inter) | Dipole-dipole (Forster) | C3/(2pi) = 10.0 GHz um^3 | 3.54 um | **225 MHz** |
| Cs-Cs (intra) | van der Waals | C6^(Cs)/(2pi) = 2364 GHz um^6 | 7.07 um | **18.9 MHz** |
| Rb-Rb (intra) | van der Waals | C6^(Rb)/(2pi) = 2036 GHz um^6 | 5.0 um | ~130 MHz (irrelevant*) |

*Rb-Rb interaction is negligible in practice because data qubits are never simultaneously excited to Rydberg states during gate operations.

### 1.4 Current Selectivity

- Ratio V_inter / V_Cs-Cs = 225 / 18.9 = **~12**
- Forster defect: delta_F/(2pi) = 2 MHz
- Forster channel: |81S_{1/2}; 77S_{1/2}> -> |80P_{1/2}; 77P_{3/2}>
- Source: Alkali Rydberg Calculator (ARC), Ref [photonics10111280]

### 1.5 Problem

The inter-species interaction V/(2pi) = 225 MHz is too weak. We want V/(2pi) > 500 MHz while keeping Cs-Cs interaction weak.

---

## 2. Why Switch from s-states to d-states?

Key reference: **Ireland, Walker & Pritchard, Phys. Rev. Research 6, 013293 (2024)** [arXiv:2401.02308]

### 2.1 Advantages of d-states over s-states

1. **Stronger interspecies Forster resonances**: d-state pairs have C3 coefficients up to 5x larger than s-state pairs at comparable n.
2. **Suppressed intraspecies interactions**: d-states naturally have weaker same-species couplings due to larger pair defects in the intraspecies channels.
3. **Lower principal quantum number**: comparable interaction strengths are achieved at lower n, giving longer Rydberg lifetimes and stronger Rabi frequencies.
4. **Larger dipole matrix elements**: the P->D transition matrix element is ~2x larger than P->S, allowing stronger two-photon Rabi frequencies at the same laser power.

### 2.2 Excitation Pathway (selection rule Delta_l = +/- 1)

```
Current (s-state):   S --(photon 1)--> P --(photon 2)--> S_Rydberg
Proposed (d-state):  S --(photon 1)--> P --(photon 2)--> D_Rydberg
```

Both pathways use the same intermediate P state. Only the frequency/polarization of the second photon changes. The ground states and intermediate states are preserved.

Specifically (from Ireland et al. Appendix B):
- Rb: 5S_{1/2} -> 6P_{3/2} (or 7P_{3/2}) -> n'D_{5/2}  (sigma+ polarization)
- Cs: 6S_{1/2} -> 7P_{1/2} (or 7P_{3/2}) -> n'D_{3/2}  (sigma- polarization)

### 2.3 Comparison: s-state vs d-state Forster resonances for Rb-Cs

From Ireland et al. Table VI (s-states) vs Table I (d-states), at R=6 um, theta=90 deg:

**Best s-state**: Rb 72s_{1/2} - Cs 70s_{1/2}
- C3_eff = 9.82 GHz um^3
- Cs-Cs C6 = -605.3 GHz um^6
- P_1r(t_pi) = 0.9992 (good blockade)
- BUT Rb-Rb C6 = -1106.6 GHz um^6 (very strong -- irrelevant for our protocol but shows s-states have strong intraspecies)

**Best d-state for multi-qubit gates**: Rb 59d_{5/2} - Cs 68d_{3/2}
- C3_eff = -14.36 GHz um^3 (1.5x stronger)
- Cs-Cs C6 = -206.1 GHz um^6 (3x weaker Cs-Cs!)
- P_1r(t_pi) = 0.9992

---

## 3. Candidate d-state Pairs with V > 500 x 2pi MHz

### 3.1 Screening Criteria

1. V_inter/(2pi) > 500 MHz at r_DA = 3.54 um, i.e., C3_eff > 22.2 GHz um^3
2. V_Cs-Cs weak at r_AA = 7.07 um (ideally < current 18.9 MHz)
3. High blockade fidelity P_1r(t_pi) > 0.999
4. Rb-Rb interaction: **irrelevant** (data qubits not simultaneously in Rydberg)

### 3.2 All Candidates Exceeding 500 MHz

Data from Ireland et al. (2024) Tables I-V, evaluated at theta=90 deg.
V_inter calculated at r=3.54 um; V_Cs-Cs calculated at r=7.07 um.

| # | Rb state | Cs state | Forster channel | C3_eff (GHz um^3) | V_inter/(2pi) (MHz) | Cs-Cs C6 (GHz um^6) | V_Cs-Cs/(2pi) (MHz) | Ratio | P_1r |
|---|----------|----------|-----------------|-------------------|---------------------|---------------------|---------------------|-------|------|
| 1 | 66d_{5/2} | 76d_{3/2} | 67p_{3/2}-74f_{5/2} | 22.84 | **514** | -692.9 | 5.5 | **93** | 0.9997 |
| 2 | 69d_{5/2} | 79d_{5/2} | 70p_{3/2}-77f_{7/2} | 26.34 | **593** | -1448.9 | 11.6 | **51** | 0.9998 |
| 3 | 70d_{5/2} | 80d_{5/2} | 71p_{3/2}-78f_{7/2} | 29.03 | **654** | -1812.1 | 14.5 | **45** | 0.9990 |
| 4 | 73d_{3/2} | 83d_{3/2} | 74p_{1/2}-81f_{5/2} | 28.54 | **643** | -2156.1 | 17.2 | **37** | 0.9995 |
| 5 | 76d_{5/2} | 87d_{5/2} | 77p_{3/2}-85f_{5/2} | 40.58 | **914** | -4482.5 | 35.8 | **26** | 0.9997 |
| 6 | 79d_{5/2} | 91d_{3/2} | 80p_{3/2}-89f_{5/2} | 49.07 | **1105** | -3596.5 | 28.7 | **38** | 0.9999 |
| 7 | 80d_{5/2} | 92d_{3/2} | 81p_{3/2}-90f_{5/2} | 51.12 | **1151** | -5756.9 | 46.0 | **25** | 0.9999 |

### 3.3 Comparison with Current Design

| Metric | Current (s-state) | Candidate #1 | Candidate #2 |
|--------|-------------------|-------------|-------------|
| Rb Rydberg | 77S_{1/2} | 66D_{5/2} | 69D_{5/2} |
| Cs Rydberg | 81S_{1/2} | 76D_{3/2} | 79D_{5/2} |
| V_inter/(2pi) | 225 MHz | 514 MHz | 593 MHz |
| V_Cs-Cs/(2pi) | 18.9 MHz | 5.5 MHz | 11.6 MHz |
| Selectivity ratio | 12 | 93 | 51 |
| Blockade fidelity | -- | 0.9997 | 0.9998 |
| Improvement factor (V) | 1x | 2.3x | 2.6x |
| Improvement factor (ratio) | 1x | 7.8x | 4.3x |

---

## 4. Recommended Pair

### Primary recommendation: Rb 69D_{5/2} - Cs 79D_{5/2}

This pair offers the best balance of strong interspecies interaction solidly above the 500 MHz threshold, weak Cs-Cs coupling, and excellent blockade fidelity.

### 4.1 Complete Level Structure

**^133Cs ancilla controls:**
- Ground states: |0>_c = |6S_{1/2}, F=3>, |1>_c = |6S_{1/2}, F=4> (unchanged)
- Intermediate: via P state (unchanged, e.g., 7P_{3/2} at 455 nm or 6P_{3/2} at 852 nm)
- Rydberg: |r>_c = |79D_{5/2}, m_j=-5/2>    **(changed from 81S_{1/2})**

**^87Rb data targets:**
- Ground states: |A>_t = |5S_{1/2}, F=1>, |B>_t = |5S_{1/2}, F=2> (unchanged)
- Intermediate: |P>_t = |7P_{3/2}> (unchanged)
- Rydberg: |R>_t = |69D_{5/2}, m_j=+5/2>    **(changed from 77S_{1/2})**

**Two-photon excitation pathways (selection rule S -> P -> D, Delta_l = +1 each step):**
- Rb: 5S_{1/2} -> 7P_{3/2} (existing laser) -> 69D_{5/2} (retune second laser, sigma+ polarization)
- Cs: 6S_{1/2} -> P intermediate (existing laser) -> 79D_{5/2} (retune second laser, sigma- polarization)

**Magnetic sublevel choice:** (m_j^Rb, m_j^Cs) = (+5/2, -5/2) for theta=90 deg operation (quantization axis perpendicular to array plane, as recommended by Ireland et al.).

### 4.2 Interaction Parameters

**Inter-species (Rb-Cs, data-ancilla):**
- Forster resonance channel: |Rb 69D_{5/2}> |Cs 79D_{5/2}> <-> |Rb 70P_{3/2}> |Cs 77F_{7/2}>
- Forster defect: delta_F/(2pi) = -19.14 MHz
- Bare channel coefficient: C_{3,k}/(2pi) = -21.97 GHz um^3
- Effective (fitted) coefficient: C3_eff/(2pi) = 26.34 GHz um^3
- At r_DA = 3.54 um: **V_dd/(2pi) = 593 MHz**

**Intra-species (Cs-Cs, ancilla-ancilla):**
- Interaction type: van der Waals
- C6^(Cs)/(2pi) = -1448.9 GHz um^6
- At r_AA = 7.07 um: **V_vdW/(2pi) = 11.6 MHz**

**Selectivity:**
- V_inter / V_Cs-Cs = 593 / 11.6 = **51** (vs 12 for current s-state design)
- Blockade leakage: P_1r(t_pi) = 0.9998 (leakage error ~2 x 10^{-4})

### 4.3 Derived Quantities (at a = 5 um, Omega_p/(2pi) = 50 MHz)

- Blockade ratio: V_dd / Omega_p = 593/50 = **11.9** (strong blockade)
- Blockade radius: r_b = (C3_eff / Omega)^{1/3} = (26340/50)^{1/3} = 8.1 um
  (comfortably encompasses nearest data-ancilla pairs at 3.54 um)
- Next-nearest data-ancilla distance: a*sqrt(5)/2 = 5.59 um (inside blockade radius --
  may need to verify if this causes issues, same as current design)

### 4.4 Rydberg State Lifetimes (estimates)

d-state lifetimes scale as ~n^3 (radiative) but are reduced by blackbody radiation at 300 K.
For n ~ 70-80, typical d-state lifetimes are:
- Rb 69D_{5/2}: tau ~ 200-400 us (comparable to current Rb 77S at 505 us)
- Cs 79D_{5/2}: tau ~ 200-400 us (comparable to current Cs 81S at 548 us)

Per-gate decay probability remains ~10^{-3} for sub-microsecond pulses, consistent with current error budget.

### 4.5 Impact on Gate Parameters

| Parameter | Current (s-state) | Proposed (d-state) |
|-----------|-------------------|--------------------|
| V_ct/(2pi) | 225 MHz | 593 MHz |
| V_ct / Omega_p | 4.5 | 11.9 |
| Blockade leakage | ~10^{-2} | ~2 x 10^{-4} |
| EIT condition (V >> Delta_EIT) | marginal | strong |
| Cs-Cs perturbation during global pulse | 18.9 MHz | 11.6 MHz |

The 2.6x stronger blockade directly improves:
1. OR gate fidelity (EIT condition more cleanly broken under blockade)
2. CCX gate fidelity (resonant transitions more effectively suppressed)
3. Reduced blockade leakage errors

---

## 5. Alternative Candidates

### 5.1 If Weakest Cs-Cs Is Priority: Rb 66D_{5/2} - Cs 76D_{3/2}

- V_inter/(2pi) = 514 MHz (just above threshold)
- V_Cs-Cs/(2pi) = 5.5 MHz (3.4x weaker than current!)
- Ratio = 93
- Forster channel: 67P_{3/2} - 74F_{5/2}
- Trade-off: marginally above 500 MHz threshold

### 5.2 If Strongest V Is Priority: Rb 79D_{5/2} - Cs 91D_{3/2}

- V_inter/(2pi) = 1105 MHz (5x current)
- V_Cs-Cs/(2pi) = 28.7 MHz (1.5x current -- worse)
- Ratio = 38
- Forster channel: 80P_{3/2} - 89F_{5/2}
- Trade-off: very strong blockade but Cs-Cs becomes comparable to current s-state

### 5.3 Additional Approach: RF-tuned Forster Resonance

Palm et al. (arXiv:2603.07958, 2026) demonstrated using AC Stark shifts from a microwave drive to tune into Forster resonance, converting 1/R^6 vdW into 1/R^3 dipolar interaction. This technique could enhance any s-state pair without changing Rydberg state choice, at the cost of requiring an additional microwave field.

### 5.4 K-Rb Forster Resonances (for A2-D interaction)

Otto, Kjaergaard & Deb, Phys. Rev. Research 2, 033474 (2020) [arXiv:2009.02004] identified an "ultrastrong" zero-field K-Rb Forster resonance with crossover distance > 100 um. This is relevant for the A2 (K) - D (Rb) interaction channel and could provide strong interactions without electric field tuning.

---

## 6. References

1. Ireland, Walker & Pritchard, "Interspecies Forster resonances of Rb-Cs Rydberg d-states for enhanced multi-qubit gate fidelities", Phys. Rev. Research 6, 013293 (2024). arXiv:2401.02308.
2. Beterov & Saffman, "Rydberg blockade, Forster resonances, and quantum state measurements with different atomic species", Phys. Rev. A 92, 042710 (2015). arXiv:1508.07111.
3. Otto, Kjaergaard & Deb, "Strong zero-field Forster resonances in K-Rb Rydberg systems", Phys. Rev. Research 2, 033474 (2020). arXiv:2009.02004.
4. Palm et al., "Enhanced Rydberg Blockade through RF-tuned Forster Resonance", arXiv:2603.07958 (2026).

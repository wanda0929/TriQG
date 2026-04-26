# Farouk et al. 2023 — Coefficient Summary

**Paper:** *Parallel Implementation of CNOTⁿ and C²NOT² Gates via Homonuclear and Heteronuclear Förster Interactions of Rydberg Atoms*, Photonics **10**, 1280 (2023).

> Core trick — and the thing you noticed — is that they deliberately **unbalance** the interactions: the control–target coupling is put on a near-resonant **dipole–dipole Förster channel** (C₃/R³, strong), while the target–target coupling stays in the **van der Waals regime** (C₆/R⁶, weak). They never call it "unbalanced"; their word is **"asymmetric"** interaction. Same idea.

---

## 1. The unbalance, in one line


| Coupling | Regime | Strength | Where it comes from |
|---|---|---|---|
| V_CT (control ↔ target) | **dipole–dipole** (Förster) | **strong** | Small energy defect δ_F between |r;R⟩ and |r';R'⟩ → near-resonant d–d |
| V_TT (target ↔ target) | **van der Waals** | weak | Same-species targets have large defects (~15 GHz) → d–d suppressed |
| V_CC (control ↔ control, only in C²NOT²) | **van der Waals** | **deceptively large** | Same-species controls, vdW, but **R_CC < R_CT** in the rhombus → V_CC can exceed V_CT |


Because V_CT scales as 1/R³ and V_TT as 1/R⁶, the ratio V_CT / V_TT scales as **R³** — so the farther you space the atoms, the more unbalanced it gets, which is exactly why they can push to "large" R_CT ∼ 6–10 µm and still have Rydberg blockade on the control↔target channel.

---

## 2. Atomic states used

| Role | Species | Ground | Intermediate |P⟩ | Rydberg (|r⟩ or |R⟩) | Lifetime |
|---|---|---|---|---|---|
| Control (heteronuclear) | ¹³³Cs | |6S₁/₂, F=3⟩, |6S₁/₂, F=4⟩ | |7P₃/₂⟩ (2nd resonance) | **|81S₁/₂, m_j = −1/2⟩** | τ_c = 548 µs |
| Target (heteronuclear) | ⁸⁷Rb | |5S₁/₂, F=1⟩, |5S₁/₂, F=2⟩ | **|6P₃/₂, m_j = 3/2⟩** (2nd resonance) | **|77S₁/₂, m_j = +1/2⟩** | — |
| Homonuclear variant | Rb / Cs | — | same | all atoms to **|77S₁/₂, +1/2⟩** (Rb) or **|81S₁/₂, −1/2⟩** (Cs) | — |

Intermediate-state lifetimes:

- Rb |5P₃/₂⟩ = 26.4 ns — 1st resonance, too lossy
- Rb |6P₃/₂⟩ = **0.131 µs** — their choice
- Cs |6P₃/₂⟩ = 30.5 ns
- Cs |7P₃/₂⟩ = **0.118 µs** — their choice

They pick the **2nd resonance** intermediate state on purpose to get a ~5× longer P-state lifetime, buying ~0.3 % fidelity (e.g. CNOT² goes from 99.43 % → 99.75 % at R_CT = 8.33 µm).

---

## 3. Interaction coefficients (the numbers you asked about)

### 3.1 Heteronuclear Cs(81S) + Rb(77S) — the headline case

Dominant Förster channel:
$$|81S_{1/2},-\tfrac12;77S_{1/2},+\tfrac12\rangle \longrightarrow |80P_{1/2},+\tfrac12;77P_{3/2},+\tfrac12\rangle$$


| Quantity                               | Symbol             | Value                     | Notes                                                     |
| -------------------------------------- | ------------------ | ------------------------- | --------------------------------------------------------- |
| Control–target dipole–dipole           | **C₃ / 2π**        | **10 GHz·µm³**            | Used in V_CT = C₃(1 − 3cos²θ)/R³, with θ = π/2            |
| Energy defect of dominant channel      | **δ_F / 2π**       | **2 MHz**                 | Tiny → near Förster resonance, boosts d–d                 |
| Target–target (Rb–Rb, both in 77S) vdW | **C₆ / 2π**        | **2036 GHz·µm⁶**          | From ARC fit, r ∈ [R_LR, 20 µm], minStateContribution = 0 |
| Le Roy radius (Rb–Rb)                  | R_LR               | 2 µm                      | Lower bound for validity                                  |
| Target–target vdW radius               | R_vdW              | 4.5 µm                    | d–d → vdW crossover                                       |
| Dimensionless coupling                 | χ_3 = C₃/(R̄³·δ_F) | **36.56** at R̄_CT = 5 µm | χ ≫ 1 → strongly d–d                                      |


### 3.2 Homonuclear Cs–Cs (comparison case)

- Cs |81S₁/₂, −1/2⟩ ↔ |81S₁/₂, −1/2⟩: **C₆ / 2π = 2364 GHz·µm⁶**, pure vdW (no Förster).

### 3.3 Table A1 — full channel list (Appendix A, for completeness)

R̄_CT = 5 µm, R̄_TT = √2·R̄_CT, χ_α = C₃^(α) / (R̄³ · δ_Fα). χ ≫ 1 means **dipole–dipole**; χ ≪ 1 means **vdW**.

**Heteronuclear Cs(81S)–Rb(77S) → |(81±n̄)P; (77±m̄)P⟩, with ∆n ≤ 2, |δ_F|/2π ≤ 2 GHz**


| α                | Channel (Rb → ; Cs →) | C₃/2π (GHz·µm³) | δ_F/2π (GHz) | χ_α         | C₆/2π (GHz·µm⁶) |
| ---------------- | --------------------- | --------------- | ------------ | ----------- | --------------- |
| 1                | 76P₃/₂; 81P₁/₂        | 10.9            | +1.87        | 4.65 × 10⁻² | 5.812           |
| 2                | 77P₁/₂; 80P₁/₂        | 5.88            | +0.214       | 2.20 × 10⁻¹ | 27.52           |
| **3 (dominant)** | **77P₃/₂; 80P₁/₂**    | **10.0**        | **+0.002**   | **36.56**   | **4570**        |
| 4                | 78P₁/₂; 79P₁/₂        | 0.105           | −0.437       | 1.93 × 10⁻³ | −0.241          |
| 5                | 78P₃/₂; 79P₁/₂        | 0.190           | −0.640       | 2.37 × 10⁻³ | −0.297          |


**Homonuclear Rb(77S)–Rb(77S) → |(77±n̄)P; (77±n̄)P⟩** (same ∆n cuts)


| α   | Channel        | C₃/2π (GHz·µm³) | δ_F/2π (GHz) | χ_α         | C₆/2π (GHz·µm⁶) |
| --- | -------------- | --------------- | ------------ | ----------- | --------------- |
| 1   | 76P₃/₂; 76P₃/₂ | 11.2            | +16.84       | 6.66 × 10⁻⁴ | 67.26           |
| 2   | 77P₁/₂; 77P₁/₂ | 4.29            | −15.40       | 2.78 × 10⁻⁴ | −96.78          |
| 3   | 77P₃/₂; 77P₃/₂ | 12.5            | −15.82       | 7.87 × 10⁻⁴ | −88.39          |
| 4   | 78P₁/₂; 78P₁/₂ | 0.060           | −46.79       | 1.29 × 10⁻⁶ | −0.006          |
| 5   | 78P₃/₂; 78P₃/₂ | 0.197           | −47.20       | 4.18 × 10⁻⁶ | −0.007          |


**The punchline:** heteronuclear χ₃ ≈ **3.7 × 10¹**, homonuclear max χ ≈ **7.9 × 10⁻⁴**. That is a **~4.6 × 10⁴** gap in the effective dipole coupling — that's the "unbalance" in one number.

### 3.4 Why heteronuclear has the tiny δ_F

Cs and Rb have *independent* level spacings, so you can tune the principal quantum numbers n₁(Cs)=81, n₂(Rb)=77 (difference = 4) to make |81S, 77S⟩ almost degenerate with |80P, 77P⟩. With two species, the two-atom Förster defect δ_F/2π drops from ~15 GHz (same species) to ~2 MHz (Cs–Rb). Four orders of magnitude — for free.

---

## 4. Laser / pulse coefficients (the dynamics side)


| Parameter | Symbol | Value | What it drives |
|---|---|---|---|
| Peak Raman Rabi on target | **Ω_p,max / 2π** | **50 MHz** | |A⟩ ↔ |B⟩ via |P⟩ on target |
| Intermediate-state detuning | **Δ / 2π** | **1200 MHz** | Single-photon detuning from |P⟩ |
| Coupling Rabi (EIT coupling) | **Ω_c / 2π** | **125 MHz** ≡ 2.5 Ω_p (sweet spot) | |P⟩ ↔ |R⟩ on target; sets EIT blocking |
| Ω_c range scanned | — | 0.15 Ω_p – 8 Ω_p (7.5–400 MHz) | Scan for EIT regime |
| Control π-pulse duration | **T_c** | **1 µs** | Excites |1⟩ → |r⟩ |
| Raman π-pulse shape | Ω_p(t) | √(16πΔ/(3T_p)) · sin²(πt/T_p) | Smooth, pulse-area condition ∫Ω_p² dt = 2πΔ |
| Raman pulse duration | **T_p** | 16πΔ / (3·Ω_p,max²) ≈ **1.28 µs** | Follows from Ω_p,max and Δ |
| Total CNOTⁿ gate time | τ | **3.28 µs** (= 2 T_c + T_p) | Sequence: π_r , Ω_p , π_r |
| Total C²NOT² gate time | τ | **5.28 µs** | Sequence: π₁, π₂, Ω_p, π₂, π₁ |


EIT regime (Müller et al. condition): Ω_c > 2 Ω_p is required to **block** the transfer when the control is in |0⟩.

---

## 5. Geometries and interatomic distances

All configurations have the quantization axis ẑ **perpendicular** to the interatomic plane (θ = π/2). Because V_CT ∝ (1 − 3cos²θ), this kills the angular factor for the CT bond direction… but they still get strong d–d because they picked m_j = −1/2 on Cs (use negative projection) — this gives the second-largest |C₃| (see Fig. 5a).

### CNOTⁿ (single control, N targets)


| N   | Spatial layout     | T-T distance R_TT in terms of R_CT | Optimal R_CT       | Best fidelity                                            |
| --- | ------------------ | ---------------------------------- | ------------------ | -------------------------------------------------------- |
| 1   | 1 target collinear | —                                  | R_LR – 5.5 µm      | **99.8 %** (Ω_c ≥ 2.5 Ω_p)                               |
| 2   | 2 targets, linear  | 2 R_CT                             | ~6 µm              | **99.75 %** (at R = 8.33 µm, 2nd-resonance intermediate) |
| 3   | isoceles triangle  | √2 R_CT                            | ~3–4 µm or 6–10 µm | **99.68 %**                                              |
| 4   | square             | √2 R_CT                            | 6–10 µm            | **99.3 %**                                               |


Homonuclear CNOT⁴ (symmetric, Rb) saturates at ~96.5 % even at Ω_c > 3.5 Ω_p and R_CT = 5 µm — limited by the T-T coupling.

### C²NOT² (two controls, two targets) — the V_CT vs V_CC problem

This is where the three-way interaction hierarchy gets interesting and **V_CC becomes the dominant error source**.

#### Geometry (Fig. 8a)

Rhombus with perpendicular diagonals. All four atoms in the xy-plane, ẑ quantization axis perpendicular to it (θ = π/2 for all bonds).

$$R_{TT} = 2\,R_{CC}, \qquad R_{CT} = \frac{\sqrt{5}}{2}\,R_{CC}$$

Inverting to use R_CT as the free parameter:

$$R_{CC} = \frac{2}{\sqrt{5}}\,R_{CT} \approx 0.894\,R_{CT}, \qquad R_{TT} = \frac{4}{\sqrt{5}}\,R_{CT} \approx 1.789\,R_{CT}$$

**Key geometric fact:** R_CC < R_CT — the controls are *closer together* than any control–target pair. This amplifies V_CC relative to what you'd naively expect.

#### Three interaction formulas (all in terms of R_CT)

| Interaction | Species pair | Regime | Formula | Coefficient |
|---|---|---|---|---|
| V_CT | Cs ↔ Rb | **dipole–dipole** | C₃ / R_CT³ | C₃/2π = 10 GHz·µm³ |
| V_CC | Cs ↔ Cs | van der Waals | C₆^(CC) · (√5/2)⁶ / R_CT⁶ = **4617 / R_CT⁶** | C₆^(CC)/2π = 2364 GHz·µm⁶ |
| V_TT | Rb ↔ Rb | van der Waals | C₆^(TT) · (√5/4)⁶ / R_CT⁶ = **62.1 / R_CT⁶** | C₆^(TT)/2π = 2036 GHz·µm⁶ |

(All V values are /2π when substituting the coefficients above.)

#### Numerical comparison at each distance

| R_CT (µm) | V_CT/2π (MHz) | V_CC/2π (MHz) | V_TT/2π (MHz) | V_CT / V_CC | V_CC / V_TT |
|---|---|---|---|---|---|
| 5 | **80** | **296** | 4.0 | 0.27 | 74 |
| **6** (sweet spot) | **46** | **99** | 1.3 | **0.47** | 74 |
| 7 | 29 | 39 | 0.53 | 0.74 | 74 |
| **7.7** (crossover) | **22** | **22** | 0.30 | **1.0** | 74 |
| 8 | 20 | 18 | 0.24 | 1.1 | 74 |
| 10 | 10 | 4.6 | 0.062 | 2.2 | 74 |

Three things to notice:

**① V_CC / V_TT = 74, always.** Both scale as 1/R_CT⁶, so the ratio is pure geometry × coefficients:
$$\frac{V_{CC}}{V_{TT}} = \frac{C_6^{(CC)}}{C_6^{(TT)}} \times \left(\frac{R_{TT}}{R_{CC}}\right)^6 = \frac{2364}{2036} \times 2^6 = 1.161 \times 64 = 74.3$$
This is a constant — independent of R_CT.

**② V_CT / V_CC grows as R_CT³.** Because V_CT ∝ 1/R³ and V_CC ∝ 1/R⁶:
$$\frac{V_{CT}}{V_{CC}} = \frac{64\,C_3}{125\,C_6^{(CC)}} \cdot R_{CT}^3 \approx 0.00217 \cdot R_{CT}^3$$

They cross at:
$$R_{CT}^* = \left(\frac{125\,C_6^{(CC)}}{64\,C_3}\right)^{1/3} = \left(\frac{125 \times 2364}{64 \times 10}\right)^{1/3} \approx \mathbf{7.7\;\mu m}$$

**③ At the paper's sweet spot R_CT ≈ 6 µm: V_CC is 2× larger than V_CT.** The control–control vdW *dominates* the control–target Förster at the actual operating point.

#### The hierarchy at the operating point (R_CT ≈ 6 µm)

$$V_{CC}\;(99\;\text{MHz}) \;>\; V_{CT}\;(46\;\text{MHz}) \;\gg\; V_{TT}\;(1.3\;\text{MHz})$$

#### Why V_CC > V_CT doesn't kill the gate (but does limit it)

V_CT and V_CC play **completely different roles** in the gate dynamics:

- **V_CT breaks EIT on targets** → enables the Raman transfer. Needs V_CT ≳ Ω_p²/Ω_c ≈ 20 MHz. At 6 µm, V_CT = 46 MHz — comfortably above threshold. ✅

- **V_CC adds phase to |11⟩ branch.** When both controls are in |r⟩ (during T_p ≈ 1.28 µs), V_CC shifts the |r,r⟩ energy:
$$\varphi_{CC} = V_{CC} \times T_p \approx 2\pi \times 99\;\text{MHz} \times 1.28\;\mu\text{s} \approx 2\pi \times 127$$
That's **~127 full phase wraps**. A tiny change in R_CT shifts φ_CC by 2π, flipping fidelity between constructive and destructive. This is why Fig. 9a shows a **sharp oscillatory dip at 8–12 µm**.

- **V_TT is negligible.** At 1.3 MHz × 1.28 µs ≈ 0.01 rad — irrelevant. The heteronuclear Förster trick does its job here.

#### Paper's evidence: Fig. 9a vs 9b

| Figure | V_CC setting | Result |
|---|---|---|
| **9a** (realistic) | V_CC ≠ 0 (99 MHz at R_CT = 6 µm) | F = 99.7 % at R_CT ≈ 6 µm; **sharp dip at 8–12 µm** |
| **9b** (hypothetical) | V_CC = 0 | F = 99.7 % over a **wide, smooth** range of R_CT |

The paper states: *"This case proves that the destructive pattern in system dynamics is a direct result of V_CC."*

#### Why can't they just go to larger R_CT where V_CT > V_CC?

At R_CT > 7.7 µm, V_CT does exceed V_CC. But V_CT also drops below the EIT-breaking threshold Ω_p²/Ω_c ≈ 20 MHz, so the gate transfer stops working. They're caught between:
- **Small R_CT:** V_CT is strong enough, but V_CC is even stronger and wraps uncontrolled phases.
- **Large R_CT:** V_CC is weak, but V_CT is also too weak to break EIT.
- **Sweet spot (R_CT ≈ 6 µm):** V_CT safely above threshold, and V_CC happens to land on a constructive phase point.

#### Implication: V_CC is the bottleneck for C²NOT²

If R_CC could be tuned independently of R_CT (i.e., break the rigid rhombus geometry), one would push R_CC larger to suppress V_CC while keeping R_CT at 6–8 µm. But the rhombus couples them: R_CC = 2R_CT/√5. This geometric constraint is what limits C²NOT² fidelity to 99.7 % vs 99.8 % for CNOTⁿ.

Alternative paths (not explored in this paper):
- **Asymmetric geometries** that decouple R_CC from R_CT.
- **Spin-echo / dynamical decoupling** sequences to refocus the V_CC phase.
- **Different Rydberg states** for the two control atoms that reduce C₆^(CC) while keeping C₃(CT) large.

---

## 6. How the paper actually uses these numbers

The full control–target Hamiltonian with **one** dominant channel (Eq. 11):

$$
\hat H_{CT} = \sum_{j=1}^{N_T} \frac{C_3(1-3\cos^2\theta_{CTj})}{R_{CTj}^3}|r\rangle\langle r'|\otimes|R\rangle_j\langle R'|*j + \sum*{j} \delta_F|r'\rangle\langle r'|\otimes|R'\rangle_j\langle R'|_j
$$

Target–target (Eq. 13), pure vdW:

$$
\hat H_{TT} = \sum_{l<k}\frac{C_6}{R_{T_lT_k}^6}|R\rangle_l|R\rangle_k\langle R|_l\langle R|_k
$$

So the "unbalanced" / asymmetric design enters through **three separate coefficients**:

1. **C₃ ≈ 2π × 10 GHz·µm³** on CT — large, resonant.
2. **δ_F ≈ 2π × 2 MHz** on CT — small enough that the channel is effectively resonant.
3. **C₆ ≈ 2π × 2036 GHz·µm⁶** on TT — pure vdW, no Förster boost.

### 6.1 What does "Rydberg blockade" mean in this paper?

Rydberg blockade is an umbrella term for "use V between Rydberg atoms to prevent some process." There are (at least) two distinct schemes both called Rydberg blockade in the literature; Farouk et al. implement the second one:


| Scheme | Who | Mechanism | Condition on V |
|---|---|---|---|
| **Direct** | Jaksch–Cirac–Zoller–Lukin 2000 | Both atoms driven to \|r⟩; V shifts \|rr⟩ out of resonance → only one can be excited | **V ≫ Ω_r** |
| **EIT-based** (this paper) | Müller–Lesanovsky–Weimer–Büchler–Zoller 2009 | Control: bare π-pulse \|1⟩→\|r⟩. Target: EIT dark state (\|P⟩↔\|R⟩ dressed by Ω_c) blocks Raman. V_CT shifts \|R⟩, breaks EIT, unblocks Raman. | **V_CT ≳ Ω_p²/Ω_c** (much looser) |


Both are Rydberg blockade; the interaction that blocks is still V_CT between Rydberg atoms. But the EIT variant has a **dramatically looser** requirement on V_CT, and that's what makes the "push R_CT out to ~8–10 µm" idea actually work.

### 6.2 The numbers at the operating point

At the sweet spot R_CT = 8.33 µm (where they report F = 99.75 % for CNOT²):
$$V_{CT} = \frac{C_3}{R^3} = \frac{2\pi \cdot 10\ \text{GHz}\cdot\mu\text{m}^3}{(8.33\mu\text{m})^3} \approx 2\pi \cdot 17\ \text{MHz}$$


| Candidate blockade condition        | Required V_CT  | Actual V_CT | Satisfied?       |
| ----------------------------------- | -------------- | ----------- | ---------------- |
| Direct Rydberg blockade: V_CT ≫ Ω_c | ≫ 2π · 125 MHz | 2π · 17 MHz | ❌ (off by ~7×)   |
| EIT blockade: V_CT ≳ Ω_p²/Ω_c       | ≳ 2π · 20 MHz  | 2π · 17 MHz | ✅ (at threshold) |


The EIT condition is *marginally* met at R_CT = 8.33 µm, which matches why they lose some fidelity there (99.75 % vs 99.8 % closer in). At R_CT = 5 µm, V_CT ≈ 2π · 80 MHz, which is 4× above the EIT threshold — deep in the blockade regime.

### 6.3 Where V_CT > Ω_p²/Ω_c comes from (sketch)

In the target's rotating frame, |A⟩, |B⟩, |R⟩ are at 0 and |P⟩ is at −Δ. The bright combination (|A⟩+|B⟩)/√2 forms a Λ with |R⟩ via |P⟩. Its dark state is

$$|D\rangle \propto \Omega_c(|A\rangle+|B\rangle)/\sqrt 2 -\sqrt 2\Omega_p|R\rangle$$

and as long as the two-photon resonance (E_{bright} = E_{|R⟩}) holds, |D⟩ is completely decoupled from both lasers — the Raman π-pulse does nothing. V_CT shifts |R⟩ by exactly V_CT, breaking two-photon resonance. The threshold for this to actually break the dark state is set by the dark state's own AC-Stark protection scale, which is ~Ω_p²/Ω_c. So:

$$\boxed{V_{CT}\gtrsim\Omega_p^2 / \Omega_c}$$

with Ω_p = 2π·50 MHz, Ω_c = 2π·125 MHz → threshold ≈ 2π · 20 MHz. That's the real operating point of the blockade, **not** V_CT vs Ω_c.

Figure 4 of the paper is the empirical signature: at Ω_c = 0.15 Ω_p the EIT regime isn't reached and blocking fails; at Ω_c = 2 Ω_p blocking is marginal; at Ω_c = 8 Ω_p blocking is clean. It's the Ω_c/Ω_p *ratio* — i.e., the EIT depth — that gates the blocking, which is the fingerprint of scheme (B).

---

## 7. TL;DR for your TriQG notes

### Interaction coefficients
- **C₃ / 2π = 10 GHz·µm³** — Cs(81S)–Rb(77S) control–target dipole–dipole (dominant Förster channel, θ = π/2, m_j(Cs) = −1/2)
- **δ_F / 2π = 2 MHz** — Förster energy defect (tiny → near-resonant → strong d–d)
- **C₆^(TT) / 2π = 2036 GHz·µm⁶** — Rb(77S)–Rb(77S) target–target vdW
- **C₆^(CC) / 2π = 2364 GHz·µm⁶** — Cs(81S)–Cs(81S) control–control vdW
- χ₃ ≈ **37** (heteronuclear) vs **< 10⁻³** (homonuclear) — ~4–5 orders-of-magnitude gap in effective d–d coupling

### Laser / pulse parameters
- **Ω_p,max / 2π = 50 MHz**, **Ω_c / 2π = 125 MHz** (= 2.5 Ω_p), **Δ / 2π = 1200 MHz**
- **T_c = 1 µs**, **T_p ≈ 1.28 µs**
- Total CNOTⁿ gate ≈ **3.28 µs**, C²NOT² ≈ **5.28 µs**

### Blockade condition
- This is the **Müller 2009 EIT-based** Rydberg blockade, not Jaksch–Lukin 2000 direct blockade.
- Operating condition: **V_CT ≳ Ω_p²/Ω_c ≈ 2π · 20 MHz** (not V_CT ≫ Ω_c = 125 MHz).
- Sweet spot: **R_CT ≈ 6–10 µm** for CNOTⁿ.

### Three-way hierarchy for C²NOT² (at R_CT ≈ 6 µm)

$$V_{CC}\;(99)\;>\;V_{CT}\;(46)\;\gg\;V_{TT}\;(1.3) \qquad [\text{MHz}]$$

- **V_TT suppressed** ← heteronuclear Förster trick ✅
- **V_CT above EIT threshold** (46 > 20 MHz) ✅
- **V_CC is the dominant error source** — wraps ~127 phase cycles during gate, creates oscillatory fidelity dips at 8–12 µm ⚠️
- V_CT = V_CC crossover at **R_CT* ≈ 7.7 µm**
- V_CC / V_TT = **74** (constant, independent of distance)
- C²NOT² fidelity capped at **99.7 %** by V_CC; removing V_CC gives smooth, wide 99.7 % over all R_CT (Fig. 9b vs 9a)


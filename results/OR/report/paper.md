# Multiqubit Toffoli gates and optimal geometry with Rydberg atoms 

Dongmin Yu[1] ,[∗] Han Wang[1] ,[∗] Jin-ming Liu[1] , Shi-Lei Su[2][,][‡] , Jing Qian[1][,][†] , and Weiping Zhang[3][,][4][,][5] 

> 1State Key Laboratory of Precision Spectroscopy, Department of Physics, School of Physics and Electronic Science, East China Normal University, Shanghai 200062, China 

> 2School of Physics, Zhengzhou University, Zhengzhou 450001, China 

> 3School of Physics and Astronomy, Tsung-Dao Lee Institute, Shanghai Jiao Tong University, Shanghai 200240, China 

> 4Collaborative Innovation Center of Extreme Optics, Shanxi University, Taiyuan, Shanxi 030006, China and 

> 5Shanghai Research Center for Quantum Sciences, Shanghai 201315, China† 

Due to its potential for implementing a scalable quantum computer, multiqubit Toffoli gate lies in the heart of quantum information processing. In this article, we demonstrate a multiqubit blockade gate with atoms arranged in a three-dimension spheroidal array. The gate performance is greatly improved by the method of optimizing control-qubit distributions on the spherical surface via evolutionary algorithm, which leads to an enhanced asymmetric Rydberg blockade. This spheroidal configuration, not only arises a well preservation for the dipole blockade energy between arbitrary control-target pairs, which keeps the asymmetric blockade error at a very low level; but also manifests an unprecedented robustness to the spatial position variations, leading to a negligible position error. Taking account of intrinsic errors and with typical experimental parameters, we numerically show that a C6NOT Rydberg gate can be created with a fidelity of 0.992 which is only limited by the Rydberg state decays. Our protocol opens up a new platform of higher-dimensional atomic arrays for achieving multiqubit neutral-atom quantum computation. 

## I. INTRODUCTION 

Rydberg atoms serve as a reliable platform for studying quantum computing and quantum simulation because of their strong and tunable interactions, which can block the excitation of surrounding atoms in the vicinity of a preexcited atom [1–3]. Via this so-called Rydberg blockade mechanism, versatile quantum gates can be created [4–8] which manifest as basic logic-calculation units for universal quantum computation [9, 10]. Among existing Rydberg-mediated quantum gates, a multiqubit Toffoli (CnNOT) gate is an important family member, which can offer an efficient implementation of Grover quantum search algorithm to speedup the searches on a programmable quantum computer [11] or to extend into any dimensional quantum systems [12]. A conventional threequbit Toffoli (C2NOT) gate can be implemented in a 1D array [13, 14] where two outer control atoms constrain the behavior of middle target atom with strong control-target interactions. However such a linear model is unsuited for engineering a multiqubit gate because one target atom can not be simultaneously manipulated by two nearest-neighbor control qubits due to the blockade mechanism [15]. Therefore previous contributions to a multiqubit Toffoli gate often rely on the assembly of several elementary gates [16–20] or the parallel operation on some clusters of atoms in a 1D array of optical tweezers [21]. Direct execution of multiqubit Toffoli gates (n ≥ 3) 

> ∗ First Author and Second Author contribute equally to this work. 

> † jqian1982@gmail.com 

remains a big challenge both in theory and experiment. 

To date, several studies for multiqubit gates use the way of adiabatic passages in which the evolution of states can be performed by obeying a multiqubit dark eigenstate with complex optimal pulses [22, 23]. An alternative way for this target depends on nonadiabatic holonomic quantum computation showing a C3NOT gate with an error of 0.0018 [24]. Another prominent idea to the realization of multiqubit Rydberg gates adopts asymmetric blockade as proposed in [25], in which there exists a large separation of scales between different types of Rydberg interactions [26–29]. However we note that, the asymmetric interaction condition breaks easily when the number of qubits is enlarged, especially for atoms arranged in 1D or 2D arrays where distant control-target interaction suffers from a dramatic decrease. Recently J. Young and coworkers propose a 2D multiqubit gate by placing many control and many target atoms at the same time, in which the strong control-control and controltarget interactions can be engineered via extra microwave fields, leading to a perfect asymmetric blockade [29]. But this scheme is still unsuitable for implementing individual multicontrol [30–32] or multitarget [33] gates due to the absence of strong and tunable interactions between distant atoms. 

In the present work, inspired by the development of defect-free atom arrays from 2D to 3D platforms where arbitrary atoms can be arranged expectantly in space [34–39] (a recent work has reported mixed-species atom arrays with arbitrary geometry [40]), we propose a scheme for implementing CnNOT gates with atoms individually arranged in a 3D spheroidal atomic array. As illustrated in Fig.1, we consider a single target atom(in 

2 

green) located at the center and n control atoms(in red) on the surface. Such a 3D atomic array can be treated as an assembling of multi-layer 2D lattices and easily achieve single-site Rydberg addressing [41]. Compared to the existing asymmetric-blockade-based protocols [25– 29], our scheme benefits from an optimal 3D configuration to maximize the asymmetry of blockade, representing an unprecedented robustness to the 3D atomic position variations. The synthetic interplay between interatomic Rydberg-Rydberg interactions and the optimal geometry results in a huge asymmetric blockade, making the gate imperfection dominated by a intrinsic decay error, and the asymmetric blockade error can be suppressed to a negligible level. Our results show that a simple estimate of decay errors gives rise to an acceptable fidelity of 0.9537 for a multiqubit C12NOT gate. This 3D Rydberg quantum gates can serve as a new gate-unit for parallel operation in 3D optical tweezers, promising for scalable quantum computation with more flexibility. 

## II. MAXIMIZING ASYMMETRIC BLOCKADE VIA OPTIMIZATION 

To achieve desirable asymmetric interactions we adopt optimization with evolutionary algorithm. Since each atom contains two Rydberg states |pj,t⟩ and |sj,t⟩ with subscript j(t) for the control(target) atom, we consider the interaction between atoms in product states |pjst⟩ and |sjpt⟩ is of resonant dipole-dipole feature [42] 

**==> picture [238 x 26] intentionally omitted <==**

with C3[sp][(][θ][c] j[t][)] = C3(1 − 3 cos[2] θcjt) and C3 = |µsp|[2] /(8πǫ0) (see Appendix A for more details). Note that the spheroidal structure can preserve the controltarget distance Rct unchanged so Ucj t only depends on the polarizing angle θcjt. On the other hand there is also a vdWs interaction between two control atoms, given by [43] 

**==> picture [209 x 27] intentionally omitted <==**

with Rcj cj′ = |rj −rj[′] | the control-control distance, which is a vdWs energy shift of the pair state |pjpj′ ⟩ considered arising from a second-order approximation of the nonresonant dipole-dipole interaction. 

Below we focus on how to achieve best asymmetric blockade using[87] Rb Rydberg states: |st,j⟩ = |(m + 1)S1/2, mj = 12[⟩][,][|][p][t,j][⟩][=][|][mP][3][/][2][, m][j][=] 23[⟩][as][in][[26].] The C3 coefficient for m = 60 is about C3/2π = 4.194 GHz·µm[3] and C6/2π = −12.0 GHz·µm[6] calculated by the ARC open source library [44]. Since Ucjt ∝ (1 − 3 cos[2] θcjt) and Ucjcj′ ∝ 1/Rc[6] jcj′[,][the][condition][of] strongly asymmetric interactions Ucjt ≫ Ucj cj′ i.e. any control atom can block the excitation of the target atom 

**==> picture [245 x 207] intentionally omitted <==**

FIG. 1. Realization of spheroidal multiqubit Toffoli gates in a 3D atomic array. Upper panels: optimized distributions of n control atoms on the surface accompanied by the optimal geometries explicitly shown below. Main panel: Amplification of the C6NOT gate with four control atoms c1∼4 at the equatorial plane and two control atoms c5∼6 at the south and north poles. The target atom is placed at the center. (a) Distribution of c2 after 10[4] optimizations via evolutionary algorithm. (b-c) The atomic energy levels as well as the atomlight interactions. In the presence of a static electric field, we fix the quantization axis along +ˆz to simplify the optimization, which arises an angular-dependent dipole-dipole interaction Ucj t(θcj t) for each control-target atomic pair, while the control-control interaction Ucj c′j[(][|][r][j][−][r][j][′][|][)][is][of][vdWs][-type] which depends on the intraspecies distance Rcj c′j[.] 

without blocking other control atoms, can be readily met if Rcj c[′] j[is][appropriate.] In Appendix B we verify the establishment of strong asymmetric interactions by calculating the leakage error due to nonresonant Rydberg couplings nearby. We also note that the asymmetry increases for small principal quantum number because the coefficient C3(C6) scales as ∼ m[4] (∼ m[11] ) [25]. Lowering m can realize a 3D quantum gate with more control qubits (details in Sec. V). 

Here, in order to maximize asymmetric blockade, we have to optimize the spatial positions of all control atoms accompanied by a dipole-angle optimization. For arbitrary control atom cj a factor characterizing asymmetry is as 

**==> picture [161 x 13] intentionally omitted <==**

which must be maximized. Intuitively, as increasing n the vdWs interaction is enhanced so as to easily break the asymmetry. To determine maximal nmax permitted for a chosen radius Rct(the scale of 3D array), we perform a global optimization to the atomic positions via evolutionary algorithm [45]. A detailed description of optimization algorithm can be found in Appendix 

3 

C. For achieving strong asymmetric blockade, we set χj > 100 which means the minimal value of χj should satisfy min(χj ) > 100 for any cj. This limitation leads to nmax = 8 when Rct = 5.0 µm and m = 60. Several optimal geometries are shown in Fig. 1(upper panels) where the positions of control atoms denoted as red dots, are precisely obtained by sufficient optimization. This optimal structure does not depend on the coefficients C3 or C6 chosen and is stably existing. For an even n value, the structure looks more regular. Physics behind these optimal geometries can be understood by seeking for a maximal asymmetry between dipole-dipole interaction and vdWs interaction, where the potential energy of system reaches its global minimum, corresponding to a maximal magnitude of dipole-dipole interaction. This specific geometry is formed by a competition between attractive vdWs interactions and inhomogeneous dipole-dipole interactions which is discussed in Appendix D. Fig.1(a) represents an amplified position distribution of c2 under 10[4] optimization. They are extremely condensed in space, confirming the accuracy of algorithm. 

## III. GATE PERFORMANCE AND DECAY ERROR 

As examples we investigate the gate performance of an optimal 3D C6NOT gate. The effective non-Hermitian Hamiltonian including the dissipative dynamics, is expressed as 

**==> picture [192 x 26] intentionally omitted <==**

with k the indices of Rydberg levels |pj,j′,t⟩ and |sj,t⟩, and the Hamiltonians 

**==> picture [226 x 72] intentionally omitted <==**

represent the atom-light couplings and the atom-atom Rydberg interactions. Ωc(t) is the Rabi frequency for the control(target) atoms. To characterize the gate performance we calculate the average gate fidelity 

**==> picture [215 x 21] intentionally omitted <==**

by solving the stochastic Schr¨odinger equation subject to arbitrary computational basis |Ψ⟩ [46]: 

**==> picture [167 x 13] intentionally omitted <==**

During each time interval δt, one generates a random number δp and compares it with the instantaneous population on Rydberg states. If δp is larger, the system 

**==> picture [245 x 111] intentionally omitted <==**

**----- Start of picture text -----**<br>
1 0.9990 0.9985 theoretical prediction<br>0.9990 0.9983 0.9958 0.9945 MC<br>0.995 0.9935<br>0.9956 0.9920<br>0.9943 0.9906<br>0.9932 0.9892<br>0.99 0.9919<br>0.9904<br>0.9891<br>0.985<br>1 2 3 4 5 6 7 8<br>n<br>¯ Fn<br>Fidelity<br>**----- End of picture text -----**<br>


FIG. 2. Average gate fidelity F[¯] n as a function of the control atomic number n, estimated by the numerical Monte Carlo method(red triangles) and the theoretical expression(blue stars). 

will evolve by obeying the Schr¨odinger equation (8); otherwise one generates a random Rydberg excitation via a quantum jump [47]. The total random number is tdet/δt where tdet = 2π/Ωc + 3π/Ωt is the gate duration. By initializing 2[n][+1] input states, Ψ[¯] out denotes the average output at t = tdet after 500 stochastic evolutions and ρet is an etalon matrix. In addition, the operator 

**==> picture [181 x 13] intentionally omitted <==**

indicates the spontaneous population decay of Rydberg levels, in which the decay rates are Γk = Γp for k = pj,j′,t and Γk = Γs for k = sj,t. 

By performing further calculations, for n = 6 we find (min(Ucjt ), max(Ucj c[′] j[))][/][2][π][= (33][.][552][,][ 0][.][096) MHz,][lead-] ing to the asymmetry: min(χj) = 349.5 > 100. Such a huge asymmetry can keep the intrinsic asymmetric error originating from imperfect control-target(control) (anti)blockade at a very low level < 10[−][5] . In turn we extend this asymmetric condition to a more generalized form, as 

**==> picture [200 x 13] intentionally omitted <==**

which is also related to relevant pulse strengths [25]. Based on Eq.(10) we assume 

**==> picture [203 x 12] intentionally omitted <==**

throughout the paper. The decay rate is Γs = 5.0 kHz and Γp = 3.4 kHz in a cryogenic environment [48], we find a gate fidelity of F[¯] 6 = 0.9920 which is mainly constrained by the decay error from Rydberg levels(the decay error is about 8 × 10[−][3] estimated by Eq.(12)). The overall gate time is tdet ≈ 3π/Ωt = 894 ns. Detailed description of the gate operation can be found in Appendix E. In Fig.2 we show that the average gate fidelity F[¯] n(red triangles and texts, estimated by MC) decreases with the control atom number n. For comparison, it is instructive to recall the decay-error expression [26] 

**==> picture [198 x 23] intentionally omitted <==**

4 

by which the gate fidelity can be analytically obtained according to F[¯] n ≈ 1 −En,se(blue stars and texts). A good agreement is observable between the theoretical and numerical predictions which confirms that other intrinsic asymmetric errors including blockade error and antiblockade error, are both negligible due to the huge asymmetry in our scheme. 

## IV. RESILIENCE TO POSITION VARIATIONS 

Owing to the finite temperature which leads to atomic position variations in the optical trap, the interatomic interaction strength is slightly different for each measurement. This so-called position error could catastrophically break the implementation of Rydberg antiblockade(RAB)-based gates which depend on a severely modified RAB condition [49–52]. Although the excitation annihilation as reviewed in [53] or transition slow-down effect [54] makes blockade gates benefited from a robustness against interaction fluctuations, most current achievements are still constrained to fewer-qubit gates [55, 56] because the blockade strength decreases significantly for two distant atoms. Here we express the control(target) atom position as 

**==> picture [249 x 200] intentionally omitted <==**

FIG. 3. Imperfection of the gate fidelity based on two different geometries vs the position variations along (a-b) xˆ and (c-d) zˆ directions. The standard deviations are σy,z = 0.27 µm in (a-b) and σx,y = 0.27 µm in (c-d). Each point denotes an average of 500 measurements. (a,c)(or (b,d)) are obtained from an optimal 3D C6NOT gate(a 2D honeycomb-type (6+1) CNOT gate). The shadings indicate a maximal position error during the calculation. 

**==> picture [171 x 12] intentionally omitted <==**

where r0,j = (Rct, θcjt, φj ) is obtained by optimization and r0,t = (0, 0, 0). The displacements δrj(t) originating from thermal motion of atoms, can be modeled as a 3D Gaussian function with widths σx,y,(z) = �kBTa/mwx,y,[2] (z)[for][radial(axial)][localizations.] Inspired by the experimental data in Ref. [57] we consider two cases: σx ∈ [0, 2.0] µm, σy,z = 0.27 µm and σz ∈ [0, 2.0] µm, σx,y = 0.27 µm. For Rb atoms held at a low temperature Ta = 10 µK, the optical trap with frequencies 2π × 18.22 kHz and 2π × 2.46 kHz arise a position uncertainty of 0.27 µm and 2.0 µm, respectively. To estimate the errors from 3D position variations we also use the way of stochastic Schr¨odinger equation and obtain the numerical solution by averaging over 500 independent trajectories. 

The numerical solutions in Fig.3(a) and (c) indicate that the 3D gate protocol can show an unprecedented robustness to the fluctuated interactions in all directions. Because in a 3D optimal configuration the position variations of atoms can be partially overcome keeping the infidelity at a small level of 10[−][4] . In contrast, arranging (6 + 1) atoms in a 2D honeycomb lattice will lead to a clear enhancement of the infidelity as shown in Fig.3(b) and (d). Especially for the radial fluctuation the imperfection dramatically increases with σx, agreeing with previous results [57–59]. Note that a 1D chain model can not preserve the asymmetric blockade condition so as to be unable to engineer a multiqubit quantum gate. Other technical imperfections such as the sensitivity to motional dephasing, laser intensity noise and laser phase 

noise would be discussed in Appendix F. A specific discussion for the leakage error due to off-resonantly coupled Rydberg pair states, will be given in Appendix B. 

## V. LARGE-SCALE MULTIQUBIT GATE 

This multiqubit Toffoli gates can be treated as a new calculation unit for large-scale quantum information processor [60]. Compared to traditional fewer-qubit gates [61–65] our protocol benefits from an optimal 3D geometry with arbitrary n(< nmax) control atoms. To fully determine maximal number nmax we graphically study it by tuning the principal quantum number m and the spherical radius Rct simultaneously. A shown in Fig.4 it is clear that nmax increases by lowering m since the asymmetry of interactions increases then, yet at the expense of the gate fidelity. Because the decay rates Γs(p) grow at the same time, leading to a larger decay error. On the other hand, we find nmax has a dramatic increase with the radius Rct. Because for a same m the absolute values Ωt and Ωc which depend on the control-target interaction Ucjt[Eq. (11)], would strongly decrease if Rct is enhanced. A smaller laser Rabi frequency will elongate the gate operation time, making the decay error dominant. So for Rct = 7 µm a C16NOT gate suffers from a very low fidelity of 0.2821 when m = 45. A high-fidelity multicontrol Rydberg-blockade gate can be accomplished by a dual consideration of both the asymmetry and the Rydberg-state decay. Our theoretical estimation shows 

5 

**==> picture [247 x 136] intentionally omitted <==**

FIG. 4. Estimates of maximal control-atom number fidelitynmax(top,F¯n(blackcolor blocks)texts, estimatedtogether withby 1the−Eaveragen,se) in gatethe space of (Rct, m). For different principal quantum number m = (45, 50, 55, 60, 65, 70, 75, 80) and at a temperature of 77 K, the calculated coefficients are C6/2π = (−1, −2, −6, −12, −24, −57, −137, −288) GHz·µm[6] , C3/2π = (1.370, 1.950, 2.912, 4.194, 5.859, 7.976, 10.620, 13.873) GHz·µm[3] , the Rydberg decay rates are Γs = (11.70, 8.55, 6.53, 5.00, 4.09, 3.33, 2.76, 2.33) kHz and Γp = (6.85, 5.50, 4.26, 3.40, 2.73, 2.25, 1.88, 1.60) kHz. All parameters are taken from the ARC open source library [44]. 

that a C12NOT(12 control qubits and 1 target qubit) gate with a fidelity of 0.9537 is possible when Rct = 7 µm and m = 60. 

## VI. CONCLUSION AND OUTLOOK 

We have studied a protocol of multicontrol-qubit Toffoli gates in which all control atoms are precisely arranged on a 3D spherical surface via optimization, which ensures a best asymmetric Rydberg blockade. These optimal geometries(see Fig. 1) are obtained by performing sufficient optimization based on evolutionary algorithm. Such a spheroidal gate has many advantages. First, it allows for a perfect preservation of strong Rydberg blockade between any control-target atom pairs, avoiding the effect of dramatic reduction in blockade strength due to distant control-target atoms as in 1D or 2D arrays. Second, an efficient optimization can ensure best asymmetric Rydberg blockade leading to a negligible asymmetric blockade error. Finally and most importantly, an unprecedented insensitivity ∼ 10[−][4] to the position variations can be observed within the 3D gate due to the compensation of three-dimensional spacial fluctuations. In comparison, Ref.[57] reports a position error of 2.5 × 10[−][3] for a radius deviation of 0.16 µm via excitation annihilation mechanism. Our work shows a minimal error of 0.0463 when 12 control atoms monitor one target atom. 

The scheme for an arbitrary (n + 1)-qubit Toffoli gate can offer a direct route to multiqubit quantum computation. Upon the basis of one-step implementation to fewer-qubit quantum gates by our group recently [56], 

this 3D blockade-gate scheme can be used to reduce the number of fewer-qubit gates, greatly lowering the complexity of quantum device design [19, 20]. Other straightforward applications with multiqubit gates refer to the production of Rydberg-mediated entanglement between two atom qubits [57] or within a mesoscopic ensemble of atoms [66], and of fast quantum computation with neutral Rydberg qubits [67] . 

## ACKNOWLEDGMENTS 

This work is supported by the National Key Research and Development Program of China under Grant No. 2016YFA0302001; by the National Natural Science Foundation of China under Grants Nos. 12174106, 11474094, 11104076, 11804308, 91950112, 11174081; the Science and Technology Commission of Shanghai Municipality under Grant No. 18ZR1412800, and the ECNU Academic Innovation Promotion Program for Excellent Doctoral Students under Grant No. YBNLTS2019-023. 

## APPENDIX A: ASYMMETRIC INTERACTIONS 

Taking account of the scheme feasibility we present details on the Rydberg pair states and their interactions in order to show the establishment of asymmetric interactions. We assume the two Rydberg states of each atom which are |st,j⟩ = |61S1/2, 2[1][⟩][,][|][p][t,j][⟩][=][ |][60][P][3][/][2][,][3] 2[⟩][.][In][the] presence of a static electric field, when two atoms(control and target) are prepared in two different dipole-coupled Rydberg states such as |pj⟩ and |st⟩, the pair state |pjst⟩ is directly coupled to the same-energy state |sjpt⟩ by a resonant dipole-dipole exchange interaction. Typically this dipole-dipole interaction B between a pair of Rydberg atoms can be given by [68] 

**==> picture [211 x 24] intentionally omitted <==**

where µ1,2 stands for the electric dipole transition operators and ε0 the permittivity of a vacuum. R is the internuclear distance and R = |R|. Moreover an external static electric field E[⃗] defines the quantization axis zˆ which controls the orientation of the dipole moments relative to the separation vector R, yielding an anisotropic dipole-dipole interaction, 

**==> picture [41 x 12] intentionally omitted <==**

**==> picture [259 x 88] intentionally omitted <==**

6 

with θcj t the polarizing angle betweenˆ the internuclear axisprojectionsand theofquantizationdipole matrixaxiselementz. µv,µ(x,y,zv onto) denotesaxis x,ˆ ˆy,the ˆz and µv,± = µv,x ± iµv,y with v ∈ (1, 2). Accounting for the use of σ-polarized transition between |60P3/2, mj = 3/2⟩ and |61S1/2, mj = 1/2⟩ with respect to ∆mj = ±1, we can ignore the term µ1zµ2z which requires ∆mj = 0. Then Eq.(15) can be reorganized as 

**==> picture [245 x 99] intentionally omitted <==**

Apparently, there are three types of angular dependence in Eq.(16) while only the first term ∝ (1 − 3 cos[2] θcjt) is appropriate. This corresponds to a resonant exchange energy between states |sjpt⟩ ⇆ |pjst⟩ where ∆mj = +1 for one atom and ∆mj = −1 for the other. Other possible transitions connecting with same combinations of ∆mj = ±1, are off-resonantly coupled due to a big Stark shift via the electric field [68]. An estimation of the leakage error to the gate fidelity from these nonresonant Rydberg levels is illustrated in Appendix B. 

In the main text we expert a resonant dipole-dipole interaction strength that only varies as 

**==> picture [220 x 26] intentionally omitted <==**

where Rct = R means the two-atom separation and the interaction coefficient C3[sp] scaling as m[4] takes a complex form of 

**==> picture [198 x 25] intentionally omitted <==**

and the transition matrix element is 

**==> picture [199 x 12] intentionally omitted <==**

Finally, we can obtain the electric dipole-dipole Hamiltonian between a pair of control and target atoms, which is 

**==> picture [240 x 25] intentionally omitted <==**

On the other hand, as for two control atoms which are prepared in same Rydberg level such as |60P3/2, mj = 3/2⟩ the electric dipole-dipole interaction B only plays roles at the second-order in perturbation theory since an atomic state has a vanishing average electric dipole moment to the first-order of perturbation [69]. As a result via B the pair state |pjpj′ ⟩ = |60P3/2, 60P3/2⟩ 

is coupled to other nearby pair states of opposite parity where the energy of those states differs from that of |60P3/2, 60P3/2⟩ by a big quantity. The average effect gives rise to a second-order vdWs shift of the considered pair state |pjpj′ ⟩ scaling as ∝ C6/Rc[6] jc[′] j[where][the][coeffi-] cient C6 roughly scales as m[11] (m is the principal quantum number). Details about the influence from original nonresonant dipole-dipole coupled states would be discussed in Appendix B. Therefore, the reduced vdWs-type interaction Hamiltonian can be described by 

**==> picture [186 x 27] intentionally omitted <==**

From Eqs.(19) and (20) it is apparent that both the dipole-dipole and vdWs interactions between two Rydberg atoms are separation-dependent. To reach a huge asymmetry in the interaction i.e. C3[sp] R[(][3] ct[θ][c] j[t][)] ≫ R[6] cj cjC6 ′[,][we] have to seek for optimal distributions of all control atoms on the spherical surface, see more details in Sec.II. 

## APPENDIX B: LEAKAGE ERROR ESTIMATION 

Leakage error based on two-atom states. As illustrated in Fig.5a we consider a resonant dipole-dipole interaction between one control atom and one target atom for the |pjst⟩ ⇆ |sjpt⟩ transition. In a real implementation these two-atom pair states might still experience nonresonant dipole-dipole couplings to other undesired Rydberg pair states, resulting in a leakage error to the gate fidelity. Here the nonresonant coupling strength and the F¨orster energy defect are denoted as Bκ and δκ respectively. Our task is to find out the influence of these nonresonant couplings to the gate fidelity estimated in our protocol. In principle we should sum over all selection-rule permitted transitions over a wide range of principal quantum numbers and calculate the leakage error. Here we have checked all possible transitions from |pjst⟩ and |sjpt⟩ to other leakage states and find that the influence of a farther state can be almost negligible due to its weaker coupling strength Bκ or a larger energy defect δκ. In the calculation the factor Bκ/δκ is used to characterize the leakage strength that is proportional to the leakage error. If Bκ/δκ ≪ 1 the leakage from the resonantly-coupled state can be suppressed [63]. 

In Fig.5a and Table I(left) we show the possible transitions with the change of principal quantum number up to ±2 from the resonant pair states |pjst⟩ ⇆ |sjpt⟩. To estimate the gate error due to the leakage of population from these states, we solve the stochastic Schr¨odinger equation (8) with respect to the Hamiltonian 

**==> picture [246 x 63] intentionally omitted <==**


In this repo, I want to perform pulse simulation through the Qutip package 'https://github.com/qutip/qutip'

I will design pulse sequence to perform three qubit OR gate:
1. The gate contains two control atoms and one target atom. The control atom contains the ground state (g) and Rydberg state (r) and the system is resonantly coupled from ground state to Rydberg state through a resonant pulse "Omega_c". The target atom contains two ground states (A,B) and intermediate excited state (P) and Rydberg state (R). Two ground states coupled to the intermediate excited state through an off-resonant laser (Omega_p) with detuning (Delta). The intermediate state P is contineously coupled to the Rydberg state (R) through a constant laser (Omega_R) resonantly.
2. First, consider the coefficients of the system: I take the paper (M. Farouk et al., “Parallel Implementation of CNOTN and C2NOT2 Gates via Homonuclear and Heteronuclear Förster Interactions of Rydberg Atoms.”) as a reference source of the coefficient source. The target qubit is 87-Rb atom and the control atom is the 133-Cs atom. First, the pulse $Omega_p_amp = 2pi * 50 MHz$, the pulse $Omega_c_amp = 2pi * 50 MHz$ and the pulse $Omega_R_amp = 2.5 * Omega_p_amp$. The detuning of transition from ground state to intermediate can be $Delta = 1200 * 2pi MHz$. The energy levels are: the Cs atoms(ground state to Rydberg state |r> = |81S_(1/2), m_j = - 1/2>); the Rb atoms (ground state to intermediate excited state P = |6P_(3/2), m_j = 3/2 > and then to Rydberg state |R> = |77S_(1/2), m_j = 1/2> ). Through ARC library, we can adjust the quantization angle slightly to make the VDWs interaction strength between control qubits far larger than the dipole interaction between the control qubit and target qubit. The ground state A = |5S_(1/2), F = 1>, B = |5S_(1/2), F = 2>. The control atom |0> = |6S_(1/2), F = 3>, |1> = |6S_(1/2), F = 4>. 
3. I have already performed the pulse shape: 
   1. the control atom pulse is from 0 to $T_c$; from $T_c + T_f$ to $2 * T_c + T_f$. where T_c is the time of pulse applied on control qubit and T_c = pi / Omega_c. 
   (function Omega_c(t)
        if 0 <= t < T_c
            return omega_c_amp / 2
        elseif T_c + T_f <= t < 2 * T_c + T_f
            return -omega_c_amp / 2
        else
            return 0.0  
        end
    end)
   2. The target atom pulse is in a shape : 
   function Omega_p(t)
            if T_c <= t < T_c + T_f
                return omega_p_amp * exp(-((t - T_c - T_f/2)^3 / sigma)^2)/2
            else
                return 0.0
            end
        end
   I need to set suitable sigma and T_f value to make that the area of F = Omega_p(t)^2/(2*Delta) equals to pi. 
   3. Lifetime: the lifetime of Rydberg state R is 505 $mu s$ and the lifetime of Rydberg state r is 548 $mu s$. The lifetime of the intermediate excited state is 0.131 microseconds. The decay rate can be derived from : gamma = 1/lifetime. And we need to add all decay channels in our simulation.
   4. The Hamiltonian includs the Rydberg interaction(between control atom and target atom); the P-level detuning and the transition Hamiltonian under all the pulses. 
4. Finally, I want the mesolve function in the Qutip package to simulate. I want: 1. the pulse shape. 2. the population transfer among computation basis. 3. the fidelity of the three qubit gate. 

5. In the next step, I also want to probe into the optimal pulse sequence....

I also want to constuct a three qubit CCX gate, too. 
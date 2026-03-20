In this repo, I want to perform pulse simulation through the Qutip package 'https://github.com/qutip/qutip'

I will design pulse sequence to perform three qubit ccx gate:

The two control qubit(with 0,1,r energy levels) and one target qubit (with A,B,P,R energy levels) do not change.
The only change is another pulse sequence for  the new ccx gate:

1. We apply the Omega_cc_amp = 2pi * 100 MHz which coupled the state 0 and r. and Omega_t_amp = 2pi * 50 MHz. 

2. I have already performed the pulse shape: 
   1. the control atom pulse is from 0 to $T_cc$; from $T_cc + 3 * T_t$ to $2 * T_c + 3 * T_t$. where T_cc is the time of pulse applied on control qubit and T_cc = pi / Omega_cc_amp and the T_t is the time of pulse applied on the target qubit with T_t = pi / Omega_t_amp 
3. function Omega_cc(p, t)
        if 0 <= t < p.T_cc
            return Omega_cc_amp / 2
        elseif p.T_cc + 3 * p.T_t <= t <= tmax
            return -Omega_cc_amp / 2
        else
            return 0.0
        end
    end

    # Target: |1> <-> |R> (first and third sub-pulse)
    function Omega_t1(p, t)
        if p.T_cc <= t < p.T_t + p.T_cc
            return Omega_t_amp / 2
        elseif 2 * p.T_t + p.T_cc <= t <= 3 * p.T_t + p.T_cc
            return Omega_t_amp / 2
        else
            return 0.0
        end
    end

    # Target: |0> <-> |R> (second sub-pulse)
    function Omega_t2(p, t)
        if p.T_t + p.T_cc <= t < 2 * p.T_t + p.T_cc
            return Omega_t_amp / 2
        else
            return 0.0
        end
    end

   1. The Hamiltonian includs the Rydberg interaction(between control atom and target atom); and the transition Hamiltonian under all the pulses. 
4. Finally, I want the mesolve function in the Qutip package to simulate. I want: 1. the pulse shape. 2. the population transfer among computation basis. 3. the fidelity of the three qubit gate. 

5. In the next step, I also want to probe into the optimal pulse sequence....


# OR Gate Average Gate Fidelity Record

This record keeps both the previously recorded Hanning-pulse fidelity data and the current Gaussian-pulse fidelity data from the terminal.

## Previous Run: Hanning Pulse

Source script: `examples/or_average_gate_fidelity.py`

### Hanning Summary

- Average gate fidelity: `F_bar = 0.997716`
- Gate infidelity: `1 - F_bar = 2.28e-03`
- Effective two-photon pulse area: `0.7854`
- Target pulse area reference: `pi = 3.1416`
- Number of control qubits: `2`
- Number of computational basis inputs: `8`
- Total gate time: `673.060 ns`

### Hanning Parameters From `examples/or_average_gate_fidelity.py`

### Drive amplitudes

- `omega_c_amp = 2 * pi * 50` MHz
- `omega_p_amp = 2 * pi * 70` MHz
- `omega_R_amp = 2.5 * omega_p_amp = 2 * pi * 175` MHz

### Interaction parameters

- `delta = 2 * pi * 1200` MHz
- `V_ct = 2 * pi * 500` MHz

### Pulse timing

- `T_c = pi / omega_c_amp = 0.010000 us = 10.000 ns`
- `T_f = 0.32653 us = 326.530 ns`
- `t_total = 2 * T_c + 2 * T_f = 0.673060 us = 673.060 ns`

### Decay rates

- `gamma_r = 1 / 548.0 = 0.001825 MHz`
- `gamma_R = 1 / 505.0 = 0.001980 MHz`
- `gamma_P = 1 / 0.131 = 7.633588 MHz`

### Solver setup

- Time grid: `500` points from `0` to `t_total`
- Solver options: `{"store_final_state": True, "nsteps": 10000}`
- Method: `mesolve`

### Script args dictionary

```python
args = {
    "omega_c_amp": omega_c_amp,
    "omega_p_amp": omega_p_amp,
    "omega_R_amp": omega_R_amp,
    "T_c": T_c,
    "T_f": T_f,
}
```

### Hanning Per-Input State Fidelities

| Input | Fidelity |
| --- | ---: |
| `0,0,A` | `0.998841` |
| `0,0,B` | `0.998841` |
| `0,1,A` | `0.997655` |
| `0,1,B` | `0.997655` |
| `1,0,A` | `0.997655` |
| `1,0,B` | `0.997655` |
| `1,1,A` | `0.996714` |
| `1,1,B` | `0.996714` |

### Hanning Output Population Breakdown

| Input | P(A) | P(B) | P(P) | P(R) |
| --- | ---: | ---: | ---: | ---: |
| `0,0,A` | `0.9988` | `0.0007` | `0.9988` | `0.0005` |
| `0,0,B` | `0.0007` | `0.9988` | `0.9988` | `0.0005` |
| `0,1,A` | `0.0011` | `0.9977` | `0.9977` | `0.0000` |
| `0,1,B` | `0.9977` | `0.0011` | `0.9977` | `0.0000` |
| `1,0,A` | `0.0011` | `0.9977` | `0.9977` | `0.0000` |
| `1,0,B` | `0.9977` | `0.0011` | `0.9977` | `0.0000` |
| `1,1,A` | `0.0009` | `0.9967` | `0.9967` | `0.0000` |
| `1,1,B` | `0.9967` | `0.0009` | `0.9967` | `0.0000` |

### Hanning Note

These values were transcribed from the active terminal output and matched against the parameter definitions in `examples/or_average_gate_fidelity.py`; the script was not rerun for this record.

## Current Run: Gaussian Pulse

Source script: `examples/or_average_gate_fid_gaussian.py`

### Gaussian Summary

- Average gate fidelity: `F_bar = 0.997811`
- Gate infidelity: `1 - F_bar = 2.19e-03`
- Effective two-photon pulse area: `0.7854`
- Target pulse area reference: `pi = 3.1416`
- Number of control qubits: `2`
- Number of computational basis inputs: `8`
- Total gate time: `320.000 ns`
- Pulse type: `super-Gaussian`
- `sigma = 0.0014`

### Gaussian Parameters From `examples/or_average_gate_fid_gaussian.py`

#### Gaussian drive amplitudes

- `omega_c_amp = 2 * pi * 50` MHz
- `omega_p_amp = 2 * pi * 50.0 * 1.039975` MHz
- `omega_R_amp = 3.5 * omega_p_amp`

#### Gaussian interaction parameters

- `delta = 2 * pi * 500` MHz
- `V_ct = 2 * pi * 5000` MHz

#### Gaussian pulse timing

- `T_c = pi / omega_c_amp = 0.010000 us = 10.000 ns`
- `T_f = 0.15 us = 150.000 ns`
- `t_total = 2 * T_c + 2 * T_f = 0.320000 us = 320.000 ns`
- `sigma = 0.0014`

#### Gaussian decay rates

- `gamma_r = 1 / 548.0 = 0.001825 MHz`
- `gamma_R = 1 / 505.0 = 0.001980 MHz`
- `gamma_P = 1 / 0.131 = 7.633588 MHz`

#### Gaussian solver setup

- Time grid: `500` points from `0` to `t_total`
- Solver options: `{"store_final_state": True, "nsteps": 10000}`
- Method: `mesolve`

#### Gaussian script args dictionary

```python
args = {
    "omega_c_amp": omega_c_amp,
    "omega_p_amp": omega_p_amp,
    "omega_R_amp": omega_R_amp,
    "T_c": T_c,
    "T_f": T_f,
    "sigma": sigma,
}
```

### Gaussian Per-Input State Fidelities

| Input | Fidelity |
| --- | ---: |
| `0,0,A` | `0.999200` |
| `0,0,B` | `0.999200` |
| `0,1,A` | `0.997543` |
| `0,1,B` | `0.997543` |
| `1,0,A` | `0.997543` |
| `1,0,B` | `0.997543` |
| `1,1,A` | `0.996958` |
| `1,1,B` | `0.996958` |

### Gaussian Output Population Breakdown

| Input | P(A) | P(B) | P(P) | P(R) |
| --- | ---: | ---: | ---: | ---: |
| `0,0,A` | `0.9992` | `0.0006` | `0.0000` | `0.0002` |
| `0,0,B` | `0.0006` | `0.9992` | `0.0000` | `0.0002` |
| `0,1,A` | `0.0019` | `0.9975` | `0.0000` | `0.0000` |
| `0,1,B` | `0.9975` | `0.0019` | `0.0000` | `0.0000` |
| `1,0,A` | `0.0019` | `0.9975` | `0.0000` | `0.0000` |
| `1,0,B` | `0.9975` | `0.0019` | `0.0000` | `0.0000` |
| `1,1,A` | `0.0019` | `0.9970` | `0.0000` | `0.0000` |
| `1,1,B` | `0.9970` | `0.0019` | `0.0000` | `0.0000` |

### Gaussian Note

These values were transcribed from the active terminal output and matched against the parameter definitions in `examples/or_average_gate_fid_gaussian.py`; the script was not rerun for this record.

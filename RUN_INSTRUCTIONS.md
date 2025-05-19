# 🧪 Simulation Instructions

All code is written in MATLAB. Each script corresponds to a specific experiment described below.

---

## 🔁 STN-GPe Weight-Based Model

**File**: `STN_GPe_WeightBased.m`

- Fast simulation using scalar synaptic weights
- Modify `W_GPe2STN` to switch between:
  - Normal: `W_GPe2STN = 0.4`
  - Parkinsonian: `W_GPe2STN = 1.5`
- Add DBS: Enable inverted VGPe feedback injection into STN

---

## 📈 Weight Sweep Analysis

**File**: `STN_GPe_WeightSweep.m`

- Sweeps `W_GPe2STN` from 0.4 to 2.5
- Measures STN–GPe phase difference using:
  - `phase_diff_STN_GPe_range.m`

---

## ⚡ Conductance-Based DBS Model

**File**: `STN_GPe_Conductance.m`

- Includes ion channel and synaptic conductance
- Add cortical input with:
  - `generate_beta_poisson_bursts.m`
- Add DBS with:
  - `generate_biphasic_dbs.m`
- Parkinsonian: `g_GPe_STN = 1.2`
- Normal: `g_GPe_STN = 0.4`

---

## 🎯 Heterogeneous Phase-Shifted Cortical Input

**File**: `STN_GPe_heterogeneous_inputs.m`

- 25 ms phase offset between STN and GPe inputs
- Add jitter (`0.8` recommended) to enhance desync

---

## 🧬 Single Neuron Simulations

- `simulate_STN_HH.m`: STN neuron with Poisson cortical drive
- `simulate_GPe_HH.m`: GPe neuron with noisy input

Use these to explore intrinsic dynamics and tuning of input currents or conductances.

---

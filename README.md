# 🧠 Charge-Balanced Biphasic Stimulation for STN-GPe Desynchronization in Parkinson’s Disease

**MATLAB-based simulation of neuromodulation strategies for suppressing beta-band synchrony in Parkinson’s Disease.**

---

## 🚀 Project Summary

This project models the neural dynamics between the **Subthalamic Nucleus (STN)** and **Globus Pallidus Externus (GPe)** using biophysical equations and proposes stimulation strategies to desynchronize their pathological beta oscillations. These synchronized oscillations (13–30 Hz) are strongly implicated in the motor symptoms of Parkinson's Disease.

---

## 🛠 Techniques Implemented

Three primary Deep Brain Stimulation (DBS) strategies were simulated:

### 1. 🎯 Phase-Shifted Cortical Stimulation
Charge-balanced, biphasic Poisson bursts sent to STN and GPe with a 25 ms delay — disrupting synchronous feedback loops.

### 2. 🔄 Closed-Loop Voltage Follower DBS
Real-time GPe voltage is inverted and injected into STN to act as a negative feedback loop that counters pathological inhibition.

### 3. ⚡ Open-Loop High-Frequency DBS
130 Hz biphasic current injected directly into STN to override synchronized beta activity.

---

## 🧪 Simulation Models

Two modeling frameworks are used:

| Model Type         | Features                                         | Use Case                            |
|--------------------|--------------------------------------------------|-------------------------------------|
| Weight-Based       | Simplified with scalar synaptic weights          | Fast sweep analysis, phase response |
| Conductance-Based  | Full ion channel + synapse dynamics (HH-based)   | High-fidelity DBS simulation        |

---

## 📂 Key Files

All code is under `/models`. Here's a quick map:

| Script | Purpose |
|--------|---------|
| `STN_GPe_WeightBased.m` | Simulates weight-driven feedback dynamics |
| `STN_GPe_Conductance.m` | Biophysical model with DBS control |
| `STN_GPe_heterogeneous_inputs.m` | Simulates phase-shifted cortical input |
| `simulate_STN_HH.m` | Single STN neuron Hodgkin-Huxley sim |
| `simulate_GPe_HH.m` | Single GPe neuron HH model |
| `generate_beta_poisson_bursts.m` | Cortical beta input (Poisson bursts) |
| `generate_biphasic_dbs.m` | DBS waveform generator |
| `STN_Reduced_Step.m` / `GPe_Reduced_Step.m` | Core update functions |

---

## 🧠 Parkinson’s Mechanism Modeled

- **Dopamine depletion** → Reduced GPe inhibition → STN hyperactivity
- **Pathological beta oscillations** emerge in STN-GPe loop
- **Closed-loop entrainment** with cortex locks motor circuits in low-mobility state

---

## 💡 Goal

To **break synchrony** via charge-balanced electrical stimulation that mimics natural desynchronization while minimizing energy use and tissue damage.

---

## 🧾 How to Run

See [`RUN_INSTRUCTIONS.md`](./RUN_INSTRUCTIONS.md) for detailed setup and simulation guidance.

---

## 👨‍🔬 Authors

Manan Bhatt  
Mansi Sharma  
Jacob Varghese

Neurostimulation Systems · 2025

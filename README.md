# Repository for Manuscript: `Compare the Pair: Rotated vs. Unrotated Surface Codes at Equal Logical Error Rates`

## Overview

This Git repository contains:
- The data for the manuscript **`Compare the Pair: Rotated vs. Unrotated Surface Codes at Equal Logical Error Rates`**.
- Code to reproduce the plots in the manuscript using the already-simulated data.
- Code to reproduce the data from scratch.

### Purpose of the Manuscript
While the rotated surface code uses half the number of qubits as the unrotated surface code to achieve the same distance, it was not known what its qubit-saving ratio was in achieving same logical error rate, which is more practically relevant. Additionally, the low-logical-error rate scaling of physical to logical error rate had not been quantified for the surface code using rigorous numerical simulations and circuit-level noise. Our paper investigates these two questions using Monte Carlo sampling of memory experiments performed with the stabilizer simulator **[Stim](https://github.com/quantumlib/Stim)**.

---

## Getting Started

### Main Workflow
If the reader would like to understand the steps used to generate the results and plots in the paper, they are encouraged to start with the notebook `src/example.ipynb`. It goes through the steps:
1. **Generating Stim Circuits**:
    - For **standard depolarizing (SD)** and **superconducting-inspired (SI)** circuit-level noise (using CNOT gates and incorporating idling errors).
    - These circuits get stored in `src/circuits/supplementary_circuits/example_circuits`
2. **Monte Carlo Sampling** of these circuits.
3. **Plotting and Analysis**:
    - Finding thresholds.
    - Analyzing physical-to-logical error rate scaling.
    - Plotting logical error rate vs. qubit count
    - Plotting teraquops
    - Calculating the ratio of the number of qubits used by each code at equal logical error rates.
    - Plotting estimated memory times.


### Availability of Data
- The circuits sampled from to generate the plots in our main text, as well as the appendix plots on SI noise, are available in `circuits/SD` and `circuits/SI`.
- Additional circuits that are unable to be generated using the provided circuit-generating function in `example.ipynb` are included in `circuits/supplementary_circuits`
- our collected data is available in `src/collected_stats/pickled_stats`. 


## Generating our manuscript's plots
If the reader would like to generate the results and plots of our manuscript using our collected data, they can work through the jupyter notebook `src/paper_data_and_plots`.
It produces the same analysis and plots as `example.ipynb` but using our data in `src/collected_stats/pickled_stats`. 
Consequently, the analysis and plots it produces are those which are included in our manuscript.

### Generating Circuits
If the reader would like to regenerate all the circuits sampled from the main text, they can run the notebook: `src/generate_papers_circuits.ipynb`

### Note:
While `example.ipynb` shows how Monte Carlo sampling can be performed within a python environment, to actually collect statistics we ran the bash files `SD_noise_collect.bash` and `SI_noise_collect.bash`

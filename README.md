# Repository for the paper [Compare the Pair: Rotated vs. Unrotated Surface Codes at Equal Logical Error Rates](https://arxiv.org/abs/2409.14765)

## Overview

Overview of this repository's contents:
- The data for the paper [**Compare the Pair: Rotated vs. Unrotated Surface Codes at Equal Logical Error Rates**](https://arxiv.org/abs/2409.14765) in `src/collected_stats`
- Code to reproduce the plots in the paper using the already-simulated data in `src/paper_data_and_plots.ipynb`.
- Code to reproduce the data from scratch in `src/example.ipynb`.
- Bash files that were used to generate data for the paper in `src/bash_files`.

### Purpose of the Paper
While the rotated surface code uses half the number of qubits as the unrotated surface code to achieve the same distance, it was not known what its qubit-saving ratio was in achieving same logical error rate, which is more practically relevant. Additionally, the low-logical-error rate scaling of physical to logical error rate had not been quantified for the surface code using rigorous numerical simulations and circuit-level noise. Our paper investigates these questions using Monte Carlo sampling of memory experiments performed with the stabilizer simulator **[Stim](https://github.com/quantumlib/Stim)**.

---

## Getting Started

### Setup Instructions

We recommend using a **virtual environment** to manage dependencies and avoid conflicts with other projects.

#### 1. Create and Activate a Virtual Environment Called compare_venv
```
python3 -m venv compare_venv
source compare_venv/bin/activate
```
#### 2. Install dependencies
```
pip install -r requirements.txt
```
#### 3. Start JupyterLab
Run the following command then open the link it provides in your browser
```
jupyter lab
```
#### 4. Deactivate compare_venv
After working through the jupyter lab notebooks deactivate the virtual environment with
```
deactivate
```


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


### Availability of Data for the paper
- The circuits sampled from to generate the plots in our main text, as well as the appendix plots on SI noise, are available in `circuits/SD` and `circuits/SI`.
- our collected data is available in `src/collected_stats/pickled_stats`.
- circuits for other appendix figures can be generated using `example.ipynb`.
- Additional circuits that are unable to be generated using the provided circuit-generating function in `example.ipynb` are included in `circuits/supplementary_circuits`.


## Generating our manuscript's plots
If the reader would like to generate the results and plots of our manuscript using our collected data, they can work through the jupyter notebook `src/paper_data_and_plots`.
It produces the same analysis and plots as `example.ipynb` but using our data in `src/collected_stats/pickled_stats`. 
Consequently, the analysis and plots it produces are those which are included in our manuscript.

### Generating Circuits
If the reader would like to regenerate all the circuits we sampled from to produce the results in the main text and the SI noise figures in the appendix, they can run the notebook: `src/generate_papers_circuits.ipynb`. 
This will regenerate the circuits alread contained in `circuits/SD` and `circuits/SI` but saves them instead in `circuits/supplementary_circuits/example_circuits`.

### Note:
While `example.ipynb` shows how Monte Carlo sampling can be performed within a python environment, to actually collect statistics we ran the bash files `SD_noise_collect.bash` and `SI_noise_collect.bash`. Note that within these bash files `--processes` is set to 64.

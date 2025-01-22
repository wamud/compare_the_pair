# Repository for the paper [Compare the Pair: Rotated vs. Unrotated Surface Codes at Equal Logical Error Rates](https://arxiv.org/abs/2409.14765)

## Overview

### Purpose of the Paper
While the rotated surface code uses half the number of qubits as the unrotated surface code to achieve the same distance, it was not known what its qubit-saving ratio was in achieving the same logical error rate, which is more practically relevant. Additionally, the low-logical-error rate scaling of physical to logical error rate had not been found for the surface code using numerical simulations and circuit-level noise. Our paper investigates these questions using Monte Carlo sampling of memory experiments performed with the stabiliser simulator **[Stim](https://github.com/quantumlib/Stim)**.

### This repository's contents:

Data and code:
- The data for the paper [**Compare the Pair: Rotated vs. Unrotated Surface Codes at Equal Logical Error Rates**](https://arxiv.org/abs/2409.14765) in `src/collected_stats/pickled_stats`
- Code to reproduce the plots in the paper using the already-simulated data is in `src/paper_data_and_plots.ipynb`.
- Code to reproduce the data from scratch is in `src/example.ipynb`.
- Bash files that were used to collect data for the paper are in `src/bash_files`.

Circuit files:
- The circuits sampled to generate our data for the plots in the main text, as well as the appendix plots on SI noise, are available in `circuits/SD` and `circuits/SI`.
- circuits for other appendix figures can be generated using `example.ipynb` by varying the CNOT order, noise model and memory type.
- Additional circuits that are unable to be generated using the provided circuit-generating function in `example.ipynb` are included in `circuits/supplementary_circuits`.

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
Run the following command 
```
jupyter lab
```
then open the link it provides in your browser.
#### 4. Work through the jupyter notebooks
Follow the steps under 'Workflow' below

#### 5. Deactivate compare_venv
After working through the jupyter lab notebooks deactivate the virtual environment with
```
deactivate
```


## Workflow

### Main workflow
If the reader would like to redo all the steps used to sample circuits and generate the results and plots in the paper, they are encouraged to start with the notebook `src/example.ipynb`. It goes through the steps:
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

Note that in the sampling step "--processes" can be set to the number of processor cores of your computer. It is currently set to 4.


### Generating our manuscript's plots
If the reader would simply like to analyse our already-collected data and generate the results and plots of our manuscript, they can work through the jupyter notebook `src/paper_data_and_plots`.
This produces the same analysis and plots as `example.ipynb` but using our data in `src/collected_stats/pickled_stats`. 
Consequently, the analysis and plots it produces are those which are included in our manuscript.

### Generating Circuits
If the reader would like to regenerate all the circuits we sampled from to produce the results in the main text and the SI noise figures in the appendix, they can run the notebook: `src/generate_papers_circuits.ipynb`. 
This will regenerate the circuits alread contained in `circuits/SD` and `circuits/SI` but saves them instead in `circuits/supplementary_circuits/example_circuits`.

### Note:
While `example.ipynb` shows how Monte Carlo sampling can be performed within a python environment, to actually collect statistics we ran the bash files `SD_noise_collect.bash` and `SI_noise_collect.bash` in `src/bash_files`. Note that within these bash files `--processes` is set to 64.

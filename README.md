# EV-Deconvolution

This repository contains the code used for the analyses, simulations, and figures presented in **[Paper Title]**.

The project focuses on evaluating and applying multiple **cell-type deconvolution methods**, including simulations with variable sequencing depth and a clinical application to colorectal cancer (CRC) data.

---

## Repository Structure

The repository is organized as follows:

### Figure Generation

These directories contain scripts used to generate the figures presented in the paper.

- `Figure_2/` – Code to generate Figure 2  
- `Figure 5/` – Code to generate Figure 5  
- `Figures_1_6_7/` – Code to generate Figures 1, 6, and 7  
- `Figure_4.R` – R script used to generate Figure 4  
- `FIGURE_S2.R` – R script used to generate Supplementary Figure S2  

---

### Simulations

- `Variable-depth synthetic mixtures/`  
  Shell scripts used to generate **variable sequencing depth synthetic mixtures** for the simulation study.

---

### Clinical Application

- `CRC Abundance estimates/`  
  Code for the **clinical application to colorectal cancer (CRC)** data, including estimation of cell-type abundances.

---

### Deconvolution Methods

The following directories contain the scripts used to run and evaluate different deconvolution methods:

- `CIBERSORT/`
- `MuSiC/`
- `NMF/`
- `NNLS/`
- `NNLS + Quad Prog/`
- `DualSimplex/`
- `HSPE/`

Each directory contains method-specific scripts and settings used in the analyses.

---

## Requirements

- R (version X.X.X or higher)
- Required R packages are listed within individual scripts
- Bash (for simulation scripts)

---

## Usage

Scripts are intended to be run independently by directory, depending on the analysis of interest:

- **Simulations**: Start with `Variable-depth synthetic mixtures/`
- **Method comparisons**: Run scripts in the corresponding deconvolution method directories
- **Figures**: Execute scripts in the relevant `Figure_*` directories

---

## Reproducibility

All scripts were used to generate the results and figures reported in the paper. Random seeds are set where applicable to ensure reproducibility.

---

## Citation

If you use this code, please cite:

> Author(s). *Title*. Journal, Year.

---

## Contact

For questions, please contact:  
**Name** – email@domain

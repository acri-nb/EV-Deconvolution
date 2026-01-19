# EV-Deconvolution

This repository contains the code used for the analyses and figures presented in **ADD OUR PAPER LINK HERE**.

The project focuses on evaluating and applying multiple **cell-type deconvolution methods** on EV-RNA data with the overarching goal of determining the most robust algorithm. Additionally, the top performing method was used in conjonction with our created signature matrix to attempt to deconvolve a clinical EV RNA dataset pertaining to colorectal cancer.

---

## Repository Structure

The repository is organized as follows:

### Figure Generation

These directories contain scripts used to generate the figures presented in the paper.

- `Figures/` – Used to generate Figures 1 through 7
- `Supplemental figures/` – Used to generate supplementary figures 1 and 2 

---

### Mixture creating

- `Variable-depth synthetic mixtures/`  
  Shell script used to generate **variable sequencing depth synthetic mixtures** for the benchmarking study.

---

### Clinical Application

- `CRC Abundance estimates/`  
  Code for the **colorectal cancer (CRC) deconvolution analysis**, including estimation of cell-type abundances.

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


## Citation

If you use this code, please cite:
***What do I put here lol***

---

## Contact

For questions, please contact:  
***What do I put here lol***

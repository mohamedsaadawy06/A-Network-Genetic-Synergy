# A Network-Genetic Synergy Organises Neurodegenerative Atrophy Patterns

This repository contains the official code, statistical models, and supplementary data for the manuscript:
**A Network-Genetic Synergy Organises Neurodegenerative Atrophy Patterns**

## Overview

The progression of neurodegenerative diseases such as Parkinson's disease (PD) and Alzheimer's disease (AD) is driven by the cell-to-cell spread of toxic proteins. This repository provides the analytical pipelines demonstrating that regional vulnerability to atrophy is jointly determined by a **Network-Genetic Synergy**:
1. **Network Centrality:** Proximity to the disease epicenter via the brain's structural connectome.
2. **Genetic Vulnerability:** Local transcriptomic profiles (from the Allen Human Brain Atlas) that govern protein synthesis and clearance.

## Data Cohorts Supported
- **PPMI:** Parkinson's Progression Markers Initiative
- **ADNI:** Alzheimer's Disease Neuroimaging Initiative
- **OASIS-3:** Open Access Series of Imaging Studies

## Repository Structure

* `main_pipeline.py`: Feature extraction and graph metrics for PPMI.
* `adni_pipeline.py`: Feature extraction and region mapping for ADNI cohorts.
* `oasis_pipeline.py`: Feature extraction for OASIS-3 replication cohorts.
* `biophysical_model.py`: Generative ODE network diffusion model simulating longitudinal atrophy.
* `stats.py` & `oasis_stats.py`: Linear Mixed Effects (LME) models and non-parametric permutation testing.
* `tvi.py`: Calculates the Transcriptomic Vulnerability Index (TVI) from AHBA microarray data.
* `graph_metrics.py`: Computes network betweenness centrality and effective distance.
* `results/`: Contains the master feature CSVs and statistical summaries (supplementary data).

## Setup & Requirements

1. Clone the repository.
2. Install dependencies:
   ```bash
   pip install pandas numpy scipy statsmodels scikit-learn networkx
   ```
3. Ensure raw neuroimaging and clinical data are placed in the `/data` directory (not included in this repository due to data use agreements).
4. Run cohort pipelines (e.g., `python main_pipeline.py`) followed by the statistical models (`python stats.py`).

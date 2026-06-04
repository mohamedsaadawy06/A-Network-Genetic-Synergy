# Supplementary Material

This document contains the supplementary files, tables, and methodology details for the manuscript: **"A Network-Genetic Synergy Universally Shapes Neurodegenerative Atrophy Patterns Across Parkinson’s and Alzheimer’s Diseases"**.

---

## Supplementary Methods

### 1. Parcellation Coordinate System
We utilized the **Cammoun033 template parcellation** (comprising 68 cortical regions, 34 per hemisphere). This is an anatomical refinement of the Desikan-Killiany (DK) atlas. 
*   **Adjacency Matrix Construction:** Connectome weights were constructed from high-angular resolution diffusion imaging (HARDI) tractography from healthy subjects, resulting in a 68x68 matrix of consensus connectivity.
*   **Laplacian Transition Probabilities:** Transition matrices ($P_{ij}$) were computed by dividing connectivity weights by row sums, and effective distance was calculated as $D_{eff}(i,j) = -\ln(P_{ij})$.

### 2. Gene Expression Mapping
The Allen Human Brain Atlas (AHBA) was processed using standard pipelines to map post-mortem donor microarray expression values to the 68 parcellated ROIs. SNCA, GBA, and PARK2 expressions were extracted for PD, and MAPT, APP, and APOE were extracted for AD. Individual regional z-scores were calculated relative to the whole-cortex average.

---

## Supplementary Tables

### Supplementary Table 1: Subject Exclusions Log (PPMI Cohort)
To ensure high image quality, strict exclusion criteria were applied to the PPMI Parkinson's cohort:
1. Mean framewise displacement (FD) exceeding $0.5$ mm during fMRI was used as a proxy for excessive in-scanner motion.
2. Missing FreeSurfer v7 reconstructions (or failures during skull-stripping or cortical ribbon rendering).

*A total of 1 subject was excluded from the final LME model analysis:*

| Patient ID (PATNO) | Reason for Exclusion |
| :---: | :--- |
| **3869** | Missing baseline cortical thickness measurement in FreeSurfer `aparc` output. |

---

### Supplementary Table 2: Regional Correlations of Cortical Thickness vs. Disease Duration (PPMI PD)
*Univariate Spearman correlation analysis was performed for each of the 68 cortical regions to correlate baseline cortical thickness with disease duration. False Discovery Rate (FDR) corrections were applied to adjust for multiple comparisons. The top 10 regions showing the strongest correlations are listed below. Note that no individual region reached strict statistical significance after FDR correction ($p_{FDR} < 0.05$), highlighting the necessity of the multivariate mixed-effects network synergy model.*

| Parcellated Region (ROI) | Spearman Correlation ($\rho$) | Raw $p$-value | FDR-Corrected $p$-value | Significant? |
| :--- | :---: | :---: | :---: | :---: |
| **L_paracentral** | -0.2919 | 0.0012 | 0.0827 | False |
| **L_pericalcarine** | -0.2385 | 0.0087 | 0.1812 | False |
| **R_paracentral** | -0.2298 | 0.0116 | 0.1812 | False |
| **L_parsorbitalis** | -0.2285 | 0.0121 | 0.1812 | False |
| **R_superiorfrontal** | -0.2245 | 0.0137 | 0.1812 | False |
| **L_superiorfrontal** | -0.2176 | 0.0170 | 0.1812 | False |
| **R_posteriorcingulate** | -0.2145 | 0.0186 | 0.1812 | False |
| **L_rostralmiddlefrontal** | -0.1865 | 0.0413 | 0.3514 | False |
| **L_isthmuscingulate** | -0.1715 | 0.0611 | 0.4066 | False |
| **L_parstriangularis** | -0.1668 | 0.0687 | 0.4066 | False |

---

### Supplementary Table 3: Template Connectome and Transcriptomic Metrics (Descriptive Top 15 Regions)
*Descriptive metrics of the parcellated templates. Regions are ranked by Betweenness Centrality (BC), showing corresponding z-scored Transcriptomic Vulnerability Index (TVI) for Parkinson's disease and mean baseline patient thickness (mm).*

| ROI Label | Betweenness Centrality (BC) | TVI ($PD$) | Effective Distance ($D_{eff}$) | Mean Thickness (mm) |
| :--- | :---: | :---: | :---: | :---: |
| **L_superiorparietal** | 0.0837 | 4.610 | 4.692 | 2.125 |
| **L_superiorfrontal** | 0.0837 | 6.315 | 2.604 | 2.617 |
| **R_superiorparietal** | 0.0837 | 9.751 | 5.223 | 2.112 |
| **R_superiorfrontal** | 0.0801 | 5.600 | 5.703 | 2.589 |
| **L_insula** | 0.0479 | 2.111 | 1.987 | 2.967 |
| **R_insula** | 0.0574 | 5.420 | 6.282 | 2.997 |
| **L_rostralmiddlefrontal** | 0.0534 | 2.883 | 4.999 | 2.390 |
| **R_rostralmiddlefrontal** | 0.0466 | 2.356 | 5.712 | 2.338 |
| **R_precuneus** | 0.0452 | 2.336 | 5.566 | 2.375 |
| **L_precuneus** | 0.0384 | 1.213 | 2.480 | 2.361 |
| **L_lateraloccipital** | 0.0154 | -0.680 | 4.632 | 2.138 |
| **R_inferiorparietal** | 0.0172 | -0.769 | 6.282 | 2.421 |
| **R_lateralorbitofrontal** | 0.0145 | 2.182 | 7.448 | 2.719 |
| **L_precentral** | 0.0140 | -1.772 | 4.750 | 2.426 |
| **R_precentral** | 0.0140 | -1.573 | 5.565 | 2.377 |

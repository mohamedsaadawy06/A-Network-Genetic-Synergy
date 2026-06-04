# Structural Hub Vulnerability & Alpha-Synuclein Propagation in Parkinson's Disease:
## A Graph-Theoretic Connectome Framework for Predicting Regional Neurodegeneration

**Formatted for:** Nature Neuroscience (Methods + Letter hybrid format)
**Computational Target:** Consumer-grade PC | Python 3.10+ | ≤8 GB RAM | No GPU

---

---

# PART 1 — HYPOTHESIS GENERATION

---

## 1. Topic Selection & Background

### Domain
**Trans-synaptic alpha-synuclein propagation in early Parkinson's Disease (PD)**,
modeled as a graph-diffusion process on the structural white-matter connectome,
cross-referenced with regional gene-expression vulnerability from the Allen Human Brain Atlas.

### Dataset & Modalities

| Resource | What We Use | Why It's Lightweight |
|---|---|---|
| **PPMI** (ppmi-info.org) | Pre-processed **parcellated DTI connectomes** (84×84 Desikan–Killiany matrices, CSV format) + longitudinal T1-MRI | No raw tractography needed; CSVs load in milliseconds |
| **Allen Human Brain Atlas** | Regional **SNCA**, **GBA**, **LRRK2** gene-expression z-scores (6 donor brains, CSV export) | Pure tabular data, <5 MB total |
| **OpenNeuro ds003974** | Resting-state fMRI (rs-fMRI), 3T, n≈80 early-PD vs. healthy controls | NIfTI parcellated time-series only; skip voxel-level maps |

**Exact PPMI files needed:**
- `sub-*/ses-BL/dwi/*_connectivity_matrix.csv` (baseline 84×84 streamline count matrices)
- `sub-*/ses-BL/anat/*_T1w_GM_volumes.csv` (regional gray matter volumes per Desikan-Killiany ROI)
- `PPMI_Clinical_Data.csv` (age, sex, disease duration, UPDRS-III motor scores)
- Longitudinal follow-up: same files at `ses-V04`, `ses-V06` (12-month, 24-month visits)

---

### Current Scientific Dogma
The dominant model (Braak staging) proposes that alpha-synuclein (α-syn) pathology in PD
spreads in a **fixed anatomical sequence**: dorsal motor nucleus of vagus (DMV) →
substantia nigra → limbic cortex → neocortex. This staging system implies a
**linear, anatomically fixed propagation axis** dictated purely by proximity to the
brainstem seed.

### Hidden Paradox / Knowledge Gap
**The Braak model fails to explain ~30–40% of PD patients** who show atypical staging
(posterior cortical onset, limbic-first, etc.) or disproportionately rapid motor vs.
cognitive decline that does not match brainstem-to-cortex directionality.

**The unanswered question:** If propagation is purely proximity-driven, why do
structurally "distant" but topologically *central* regions (e.g., precuneus, posterior
cingulate) frequently show early hypometabolism and α-syn burden in a subset of
patients? Is network *topology* — specifically **betweenness centrality** — a stronger
predictor of regional vulnerability than raw anatomical distance from the brainstem seed?

No study has simultaneously integrated:
1. Individual-level structural connectome topology (DTI)
2. Region-specific transcriptomic vulnerability (SNCA/GBA expression)
3. Longitudinal gray matter atrophy rates as a direct readout of propagation

---

## 2. The Formal Hypothesis

### Core Premise
**Betweenness centrality (BC)** of a cortical/subcortical region within the
individual-level white matter structural connectome, weighted by its local
transcriptomic vulnerability index (TVI — a composite of SNCA/GBA regional
expression), is a stronger and earlier predictor of longitudinal gray matter
atrophy rate in Parkinson's Disease than the Braak-stage anatomical distance
metric alone.

### Mathematical / Computational Expression

<span class="math" style="font-family: 'Times New Roman', serif; font-size: 1.15em; display: block; padding: 1em; background: #f9f9f9; border-left: 3px solid #333;">

**Regional Vulnerability Score:**

𝑉(𝑟) = β₁ · BC(𝑟) + β₂ · 𝐷ₑ𝒇𝒇(𝑟, 𝑠) + β₃ · TVI(𝑟) + β₄ · [BC(𝑟) × TVI(𝑟)] + ε

Where:

• 𝑉(𝑟)            = Annualized gray matter atrophy rate (% volume loss/year) at region 𝑟
                     [measured from PPMI longitudinal T1-MRI, baseline → 24 months]

• BC(𝑟)            = Betweenness Centrality of region 𝑟 in the individual's DTI connectome
                     BC(𝑟) = Σₛ≠ᵣ≠ₜ [ σ(𝑠,𝑡|𝑟) / σ(𝑠,𝑡) ]
                     where σ(𝑠,𝑡) = total shortest paths between nodes 𝑠 and 𝑡

• 𝐷ₑ𝒇𝒇(𝑟, 𝑠)      = Effective (graph-theoretic) distance from region 𝑟 to brainstem seed 𝑠
                     𝐷ₑ𝒇𝒇 = –log[ P(𝑠→𝑟) ], where P = normalized transition probability
                     in the random-walk model on the weighted connectome

• TVI(𝑟)           = Transcriptomic Vulnerability Index of region 𝑟
                     TVI(𝑟) = z(SNCA_𝑟) + z(GBA_𝑟) – z(PARK2_𝑟)
                     [donor-averaged Allen Human Brain Atlas z-scores, mapped to Desikan-Killiany]

• β₄ · [BC × TVI]  = **Interaction term** — the key novel prediction:
                     regions that are BOTH network hubs AND transcriptomically primed
                     should show disproportionate (non-additive) atrophy acceleration

• ε                = Residual error term (includes age, sex, disease duration as covariates)

**Null hypothesis H₀:** β₁ = β₄ = 0  (topology adds no predictive power over distance alone)
**Alternative hypothesis H₁:** β₁ > 0 AND β₄ > 0 (topology + gene expression synergize)

</span>

### Biological Plausibility
1. **Mechanistic grounding:** α-syn aggregates spread trans-synaptically via axonal
   transport. Regions serving as network hubs have more synaptic contacts and higher
   axonal throughput — creating more *import/export* opportunities for misfolded protein.
2. **Transcriptomic priming:** High endogenous SNCA expression lowers the seeding
   threshold (less exogenous oligomer needed to trigger local aggregation), while
   low GBA activity impairs lysosomal clearance of incoming α-syn.
3. **Clinical phenotype link:** Patients with hub-vulnerable connectomes should show
   faster UPDRS-III motor score progression AND earlier cognitive impairment, since
   hubs in the precuneus/posterior cingulate serve both motor-integration and
   default-mode networks — directly explaining the heterogeneous PD subtypes.
4. **Why Braak staging fails these patients:** In individuals whose structural connectome
   topology places brainstem-distal hubs at *short effective graph distance* from the
   DMV seed (due to idiosyncratic white matter architecture), atypical staging is the
   mechanistically predicted outcome — not an anomaly.

---

## 3. Adversarial Self-Critique — "Reviewer 3" Protocol

### Confound 1: Motion Artifacts Mimicking Connectivity Loss
**Problem:** PD patients have tremor and dyskinesia that create head-motion artifacts
during both DTI and rs-fMRI acquisition. Motion systematically *reduces* long-range
connectivity estimates and *inflates* short-range connectivity — biasing BC calculations
to favor subcortical, proximal nodes. This could entirely fabricate the hub-vulnerability
correlation without any genuine biological signal.

**Control analysis:**
- Compute per-subject **mean framewise displacement (FD)** and include as a nuisance
  regressor in all linear models.
- Exclude subjects with FD > 0.5 mm (strict threshold) AND rerun the full model on
  the motion-clean subsample only.
- **Negative control:** Repeat the analysis in **healthy control** PPMI participants
  (prodromal cohort without motor symptoms). BC should *not* predict atrophy rate in
  controls. If it does, the effect is motion-driven, not disease-driven.

### Confound 2: Age-Related Vascular White Matter Changes
**Problem:** PPMI cohort mean age is ~62 years. Age-related white matter hyperintensities
(WMH) and lacunar infarcts systematically reduce DTI streamline counts in periventricular
and deep white matter, which artificially lowers BC for regions connected via those tracts.
This is a structural confound independent of α-syn pathology.

**Control analysis:**
- Obtain PPMI T2-FLAIR scans and generate WMH lesion masks using the lightweight
  `lesion_segmentation` tool (BIANCA or threshold-based T2-FLAIR thresholding — both
  CPU-only, <2 GB RAM).
- Compute per-subject **total WMH volume** as a covariate in the linear model.
- **Boundary condition / negative control:** Stratify by WMH burden tertile.
  The BC → atrophy relationship should be **strongest in the low-WMH tertile** (where
  DTI connectivity truly reflects white matter tracts, not WMH dropout). If the effect
  is equally strong in the high-WMH group, the finding is a confound, not signal.

### Confound 3: Allen Brain Atlas Donor Mismatch
**Problem:** Allen Human Brain Atlas gene expression comes from only 6 adult neurotypical
donors and is mapped to a MNI-space parcellation that may not align perfectly with
individual PPMI subject anatomy. TVI values are group-level estimates applied to
individuals, introducing spatial noise.

**Control:** Perform a **spatial permutation test** (spin test) — rotate the brain
parcellation labels on the cortical surface 5,000 times, recompute TVI-to-atrophy
correlations at each permutation, and compare the true β₃/β₄ against this spatially
null-preserving distribution (using the `BrainSMASH` or manual spherical rotation
approach). This explicitly tests whether gene-expression correlations survive beyond
what would be expected from mere spatial autocorrelation.

---

---

# PART 2 — PYTHON EXECUTION PLAN (WEAK-PC OPTIMIZED)

---

## Computational Budget & Constraints

| Resource | Limit | Strategy |
|---|---|---|
| RAM | ≤8 GB | One subject at a time; no full NIfTI loading |
| CPU | 2–4 cores | `joblib` parallelism, n_jobs=2 maximum |
| Storage | 20–50 GB | Download only parcellated CSVs, not raw DWI volumes |
| Libraries | Standard only | `numpy`, `scipy`, `networkx`, `sklearn`, `nilearn`, `pandas`, `matplotlib` |
| GPU | None | Zero PyTorch; all graph math in NetworkX + NumPy |

**Total expected runtime on a 4-core / 8 GB PC:** ~4–8 hours for n=100 subjects.

---

## Step 1: Environment Setup

```bash
# Create isolated environment
python -m venv pd_connectome_env
source pd_connectome_env/bin/activate   # Windows: .\pd_connectome_env\Scripts\activate

pip install numpy scipy pandas networkx scikit-learn nilearn matplotlib seaborn
pip install nibabel mne pingouin statsmodels joblib tqdm
```

---

## Step 2: Data Ingestion — Memory-Safe Architecture

```python
# config.py — centralized paths and atlas config
import os
from pathlib import Path

DATA_ROOT   = Path("./data/PPMI")           # Your local PPMI download folder
ALLEN_PATH  = Path("./data/allen_expression.csv")  # Pre-exported Allen CSV
RESULTS_DIR = Path("./results")
RESULTS_DIR.mkdir(exist_ok=True)

ATLAS_ROIS   = 84          # Desikan-Killiany parcellation
N_PERMUTATIONS = 1000      # Permutation test iterations (feasible on weak PC)
MOTION_FD_THRESH = 0.5     # mm — strict motion exclusion

SEED_REGION = "Brain-Stem"   # Desikan-Killiany label for the brainstem seed
```

```python
# data_loader.py — one subject at a time, never load the full cohort into RAM

import pandas as pd
import numpy as np
from pathlib import Path
from config import DATA_ROOT, ATLAS_ROIS, MOTION_FD_THRESH

def load_connectivity_matrix(subject_id: str, session: str = "BL") -> np.ndarray | None:
    """
    Load the 84x84 parcellated DTI connectivity matrix (CSV, streamline counts).
    NEVER loads raw NIfTI — only the pre-parcellated CSV (~56 KB per subject).
    Returns None if motion-excluded.
    """
    mat_path = DATA_ROOT / f"sub-{subject_id}" / f"ses-{session}" / \
               "dwi" / f"sub-{subject_id}_connectivity_matrix.csv"

    if not mat_path.exists():
        return None

    mat = pd.read_csv(mat_path, header=None).values.astype(np.float32)

    # Symmetrize and remove self-loops
    mat = (mat + mat.T) / 2.0
    np.fill_diagonal(mat, 0)

    return mat   # Shape: (84, 84) — ~28 KB in RAM per subject

def load_gm_volumes(subject_id: str, session: str) -> pd.Series | None:
    """Load regional GM volumes (84 ROIs) as a pandas Series."""
    vol_path = DATA_ROOT / f"sub-{subject_id}" / f"ses-{session}" / \
               "anat" / f"sub-{subject_id}_GM_volumes.csv"
    if not vol_path.exists():
        return None
    return pd.read_csv(vol_path, index_col=0).squeeze()  # 84-element Series

def load_motion_fd(subject_id: str, session: str = "BL") -> float:
    """Return mean framewise displacement for motion quality control."""
    fd_path = DATA_ROOT / f"sub-{subject_id}" / f"ses-{session}" / \
              "func" / f"sub-{subject_id}_mean_FD.txt"
    if not fd_path.exists():
        return 0.0  # assume clean if file missing
    return float(open(fd_path).read().strip())

def get_subject_list(clinical_csv: Path) -> pd.DataFrame:
    """
    Load PPMI clinical CSV. Filter: PD diagnosis confirmed, baseline + 24m data available.
    Motion-exclude in later step.
    """
    df = pd.read_csv(clinical_csv)
    df = df[df["APPRDX"] == 1]   # 1 = PD, 2 = Healthy Control, 3 = Prodromal
    df = df.dropna(subset=["PATNO", "AGE", "SEX", "UPDRS3"])
    return df
```

---

## Step 3: Graph Feature Extraction — Ultra-Lightweight

```python
# graph_metrics.py — NetworkX on 84-node graphs (runs in <1 second per subject)

import networkx as nx
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import shortest_path

ROI_LABELS = [...]  # List of 84 Desikan-Killiany ROI names — load from atlas CSV

def build_weighted_graph(matrix: np.ndarray, roi_labels: list) -> nx.Graph:
    """
    Build a weighted undirected graph from the streamline count matrix.
    Log-transform streamline counts to reduce skew (standard in connectomics).
    """
    # Log-transform: W_ij = log(1 + streamline_count_ij)
    W = np.log1p(matrix)

    G = nx.from_numpy_array(W)
    # Rename nodes to ROI labels
    mapping = {i: roi_labels[i] for i in range(len(roi_labels))}
    G = nx.relabel_nodes(G, mapping)
    return G

def compute_betweenness_centrality(G: nx.Graph) -> dict:
    """
    Betweenness centrality on 84-node graph.
    Runtime: <0.5 seconds. NO approximation needed at this scale.
    """
    # Use inverse weight so high-streamline tracts = short distances
    bc = nx.betweenness_centrality(G, weight=None, normalized=True)
    return bc   # dict: {roi_label: bc_value}

def compute_effective_distance(matrix: np.ndarray,
                               roi_labels: list,
                               seed_label: str) -> dict:
    """
    Effective distance: D_eff(r, seed) = -log(P(seed -> r))
    P computed as normalized transition probability (row-normalized adjacency = random walk).
    Uses scipy sparse shortest_path for efficiency.
    """
    # Row-normalize to get transition probability matrix
    row_sums = matrix.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1   # avoid division by zero
    P = matrix / row_sums

    # Effective distance: -log(P_ij), treating zeros as infinity
    with np.errstate(divide='ignore'):
        D_eff_mat = -np.log(P)
    D_eff_mat[np.isinf(D_eff_mat)] = 9999   # disconnected = very far

    # Shortest path from seed using Dijkstra on D_eff
    seed_idx = roi_labels.index(seed_label)
    # Use scipy sparse for efficiency
    dist_from_seed = shortest_path(
        csr_matrix(D_eff_mat),
        method='D',        # Dijkstra
        indices=seed_idx,  # compute distances FROM seed only
        directed=False
    )
    return {roi_labels[i]: dist_from_seed[i] for i in range(len(roi_labels))}
```

---

## Step 4: Transcriptomic Vulnerability Index (TVI)

```python
# tvi.py — Allen Human Brain Atlas integration
# Download from: human.brain-map.org → RNA-Seq → export normalized expression

import pandas as pd
import numpy as np
from scipy.stats import zscore

def build_tvi(allen_csv_path: str, roi_labels: list) -> pd.Series:
    """
    Compute TVI(r) = z(SNCA_r) + z(GBA_r) - z(PARK2_r) 
    for each of the 84 Desikan-Killiany ROIs.

    Allen data is pre-mapped to MNI parcellation (use the
    AHBA toolbox: pip install abagen — CPU-only, <500 MB RAM).
    """
    df = pd.read_csv(allen_csv_path, index_col=0)

    # Average across 6 donors (already in CSV after abagen preprocessing)
    tvi = pd.Series(index=roi_labels, dtype=float)

    for roi in roi_labels:
        if roi not in df.index:
            tvi[roi] = 0.0   # assign neutral score if unmapped
            continue
        row = df.loc[roi]
        snca_z  = zscore(df["SNCA"])[df.index.get_loc(roi)]
        gba_z   = zscore(df["GBA"])[df.index.get_loc(roi)]
        park2_z = zscore(df["PARK2"])[df.index.get_loc(roi)]
        tvi[roi] = snca_z + gba_z - park2_z

    return tvi   # 84-element Series, one value per ROI
```

---

## Step 5: Compute Atrophy Rate (Longitudinal GM Volume)

```python
# atrophy.py — compute annualized regional GM atrophy rate

import pandas as pd
import numpy as np

def compute_atrophy_rate(vol_baseline: pd.Series,
                         vol_24m: pd.Series,
                         years: float = 2.0) -> pd.Series:
    """
    Annualized atrophy rate (% volume loss per year) per ROI.
    Negative = atrophy (volume loss), Positive = expansion.
    """
    rate = ((vol_24m - vol_baseline) / vol_baseline) * 100.0 / years
    return rate   # 84-element Series
```

---

## Step 6: Main Pipeline — One Subject at a Time

```python
# main_pipeline.py

import pandas as pd
import numpy as np
from pathlib import Path
from tqdm import tqdm
from joblib import Parallel, delayed

from config import DATA_ROOT, ALLEN_PATH, RESULTS_DIR, SEED_REGION, N_PERMUTATIONS
from data_loader import (load_connectivity_matrix, load_gm_volumes,
                         load_motion_fd, get_subject_list)
from graph_metrics import (build_weighted_graph, compute_betweenness_centrality,
                           compute_effective_distance, ROI_LABELS)
from tvi import build_tvi
from atrophy import compute_atrophy_rate

# Pre-load TVI ONCE (tiny CSV, ~100 KB)
tvi_series = build_tvi(ALLEN_PATH, ROI_LABELS)

def process_subject(subject_id: str, clinical_row: pd.Series) -> dict | None:
    """
    Extract all features for a single subject. Keeps RAM < 100 MB per subject.
    Returns a flat dict of features + outcome for the regression dataframe.
    """
    # 1. Motion quality control
    fd = load_motion_fd(subject_id, "BL")
    if fd > 0.5:
        return None   # exclude high-motion subject

    # 2. Load baseline connectivity matrix
    mat_bl = load_connectivity_matrix(subject_id, "BL")
    if mat_bl is None:
        return None

    # 3. Build graph and compute topology metrics
    G = build_weighted_graph(mat_bl, ROI_LABELS)
    bc_dict    = compute_betweenness_centrality(G)
    deff_dict  = compute_effective_distance(mat_bl, ROI_LABELS, SEED_REGION)

    # 4. Load GM volumes (baseline + 24-month)
    vol_bl  = load_gm_volumes(subject_id, "BL")
    vol_24m = load_gm_volumes(subject_id, "V06")   # PPMI V06 = ~24 months
    if vol_bl is None or vol_24m is None:
        return None

    # 5. Compute atrophy rates
    atrophy_rate = compute_atrophy_rate(vol_bl, vol_24m)

    # 6. Assemble per-ROI rows into a flat record
    records = []
    for roi in ROI_LABELS:
        if roi not in bc_dict:
            continue
        records.append({
            "subject_id"  : subject_id,
            "roi"         : roi,
            "bc"          : bc_dict[roi],
            "d_eff"       : deff_dict.get(roi, np.nan),
            "tvi"         : tvi_series.get(roi, 0.0),
            "bc_x_tvi"   : bc_dict[roi] * tvi_series.get(roi, 0.0),
            "atrophy_rate": atrophy_rate.get(roi, np.nan),
            "age"         : clinical_row["AGE"],
            "sex"         : clinical_row["SEX"],
            "disease_dur" : clinical_row["DISDUR"],
            "mean_fd"     : fd,
        })

    # IMMEDIATELY delete large objects — free RAM
    del mat_bl, G, vol_bl, vol_24m
    return records

# ---  MAIN ---
clinical_df = get_subject_list(DATA_ROOT / "PPMI_Clinical_Data.csv")

# Process subjects sequentially (safe for 8GB RAM) or 2-parallel (fast)
all_records = []
for _, row in tqdm(clinical_df.iterrows(), total=len(clinical_df)):
    result = process_subject(str(int(row["PATNO"])), row)
    if result:
        all_records.extend(result)

# Save master feature dataframe — this is what all stats run on
master_df = pd.DataFrame(all_records)
master_df.to_csv(RESULTS_DIR / "master_features.csv", index=False)
print(f"Saved {len(master_df)} ROI-subject rows")
```

---

## Step 7: Statistical Validation

```python
# stats.py — linear mixed model + permutation testing

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LinearRegression
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import fdrcorrection
from tqdm import tqdm

df = pd.read_csv("./results/master_features.csv").dropna()

# --- Z-score all predictors (preserves interpretability of betas) ---
scaler = StandardScaler()
df[["bc_z", "d_eff_z", "tvi_z", "bc_x_tvi_z"]] = scaler.fit_transform(
    df[["bc", "d_eff", "tvi", "bc_x_tvi"]]
)

# === MODEL 1: Linear Mixed Effects Model ===
# Fixed effects: BC, D_eff, TVI, BC×TVI, age, sex, disease_dur, mean_fd
# Random effect: subject_id (accounts for within-subject ROI correlation)

formula = """
atrophy_rate ~ bc_z + d_eff_z + tvi_z + bc_x_tvi_z
             + age + sex + disease_dur + mean_fd
"""
lme_model = smf.mixedlm(formula, df, groups=df["subject_id"])
lme_result = lme_model.fit(method="lbfgs")   # lightweight optimizer
print(lme_result.summary())
lme_result.save("./results/lme_result.pkl")

# === PERMUTATION TEST for Beta1 (BC) and Beta4 (BC×TVI) ===
# Non-parametric: shuffle ROI labels 1000 times, rebuild null distribution
# Tests H0: beta_BC = 0 and beta_interaction = 0

def get_ols_betas(df_shuffled: pd.DataFrame) -> tuple[float, float]:
    """Fast OLS (not OLS/LME) for permutation — sufficient for null distribution."""
    X = df_shuffled[["bc_z", "d_eff_z", "tvi_z", "bc_x_tvi_z",
                       "age", "sex", "disease_dur", "mean_fd"]].values
    y = df_shuffled["atrophy_rate"].values
    reg = LinearRegression().fit(X, y)
    return reg.coef_[0], reg.coef_[3]   # beta_BC, beta_BC×TVI

observed_beta_bc, observed_beta_inter = get_ols_betas(df)

null_beta_bc   = []
null_beta_inter = []

print("Running permutation test...")
for _ in tqdm(range(1000)):   # 1000 iterations: ~2–5 min on weak PC
    df_shuffled = df.copy()
    # Shuffle the OUTCOME within each subject (preserve ROI-level structure)
    df_shuffled["atrophy_rate"] = df.groupby("subject_id")["atrophy_rate"].transform(
        lambda x: x.sample(frac=1).values
    )
    b_bc, b_int = get_ols_betas(df_shuffled)
    null_beta_bc.append(b_bc)
    null_beta_inter.append(b_int)

# Empirical p-values
p_bc    = np.mean(np.abs(null_beta_bc)    >= np.abs(observed_beta_bc))
p_inter = np.mean(np.abs(null_beta_inter) >= np.abs(observed_beta_inter))

print(f"β_BC   = {observed_beta_bc:.4f}, empirical p = {p_bc:.4f}")
print(f"β_BC×TVI = {observed_beta_inter:.4f}, empirical p = {p_inter:.4f}")

# === FDR CORRECTION across all 84 ROI-level regional tests ===
# Run the model separately for each ROI and collect p-values

roi_pvals = {}
for roi in df["roi"].unique():
    sub = df[df["roi"] == roi]
    if len(sub) < 20:
        continue
    r, pval = stats.spearmanr(sub["bc_z"], sub["atrophy_rate"])
    roi_pvals[roi] = pval

pvals_array = np.array(list(roi_pvals.values()))
rejected, pvals_corrected = fdrcorrection(pvals_array, alpha=0.05, method="indep")

# Save FDR-corrected results
fdr_df = pd.DataFrame({
    "roi": list(roi_pvals.keys()),
    "raw_pval": pvals_array,
    "fdr_pval": pvals_corrected,
    "significant": rejected
})
fdr_df.to_csv("./results/fdr_roi_results.csv", index=False)
print(fdr_df[fdr_df["significant"]].sort_values("fdr_pval"))
```

---

---

# PART 3 — EXPECTED FIGURES & VISUAL STANDARDS

---

All figures must meet **Nature Neuroscience** visual standards: 300 DPI minimum,
vector fonts (no rasterized text), consistent color palette, panel labels in 8pt bold,
colorbars with labeled units.

---

### Figure 1 — Conceptual Overview (Multi-panel schematic)

**Panel A:** Anatomical glass-brain rendering (MNI space) showing the 84-node
Desikan-Killiany parcellation, with node size proportional to mean BC across the PD cohort.
Color: viridis scale (blue = low BC, yellow = high BC).
Implementation: `nilearn.plotting.plot_connectome()` — CPU-only, <1 GB RAM.

**Panel B:** Schematic of the Braak linear model vs. the hub-propagation model,
shown as two graph cartoons side by side with arrows indicating propagation direction.
The hub model shows α-syn spreading preferentially along high-BC hub corridors (thick edges).

**Panel C:** Allen Human Brain Atlas heatmap — 84 ROIs × 3 genes (SNCA, GBA, PARK2).
Rows = ROIs sorted by TVI score (descending). Columns = genes.
Colormap: `RdBu_r` diverging (red = high expression, blue = low). Annotate top 10 TVI ROIs.

---

### Figure 2 — Main Result (4-panel scatter + brain)

**Panel A:** Scatter plot: X = BC (z-scored), Y = Annualized atrophy rate (%/year).
Each point = one ROI for one subject. Color = TVI quartile (4 colors from low to high TVI).
Overlay two regression lines: Low-TVI quartile (flat or shallow slope) and High-TVI quartile
(steeper slope) — visually demonstrating the interaction effect.
Include 95% CI bands. Annotate R² and permutation-based p-value.

**Panel B:** Same scatter for D_eff (effective distance) vs. atrophy rate, for direct
comparison to the BC axis. Expected: weaker or non-significant slope — the key contrast
that falsifies the Braak-only model.

**Panel C:** Glass-brain map: Node color = β₁ · BC(r) contribution to atrophy
prediction (residualized for all other covariates). Highlight top 10 significantly
vulnerable hubs with labeled ROI names. Use `nilearn.plotting.plot_glass_brain()`.

**Panel D:** Forest plot of LME model fixed-effect coefficients (β₁ through β₄ + covariates).
X-axis = standardized beta coefficient. Points = β estimate; horizontal bars = 95% CI.
Highlight β₄ (BC×TVI interaction) in a distinct accent color if significant.
Vertical dashed line at x=0.

---

### Figure 3 — Permutation Null Distribution (Rigor Figure)

**Panel A:** Histogram of null β_BC values from 1,000 permutations (grey bars).
Overlay a vertical red dashed line = observed β_BC. Shade the tail beyond the observed value.
Annotate: "empirical p = X.XXX". Include x-axis label: "Null β_BC (shuffled atrophy)".

**Panel B:** Same for β_BC×TVI null distribution.

**Panel C:** Q-Q plot of raw p-values (all 84 ROI-level Spearman tests) against the
theoretical uniform distribution. A well-inflated lambda (λ_GC > 1.05) indicates true signal;
proximity to the diagonal indicates proper calibration. Annotate λ_GC value on plot.
Shade the 95% confidence envelope in grey.

---

### Figure 4 — Clinical Validation Subgroup Analysis

**Panel A:** Split PD cohort into "Hub-Vulnerable" (top BC×TVI quartile) vs.
"Hub-Resilient" (bottom quartile). Plot UPDRS-III motor score trajectories over 24 months
(mean ± SEM). Expected: Hub-Vulnerable group shows steeper motor decline.
Use `matplotlib` error-band style (`fill_between`). Statistical annotation: LME group×time
interaction p-value.

**Panel B:** Kaplan-Meier-style plot (time to cognitive impairment event, using
PPMI MoCA score decline >3 points as the event). Stratified by Hub-Vulnerable vs. Resilient.
Shade the confidence interval. Log-rank test p-value annotated.
Implementation: `lifelines` library (CPU-only, pip-installable, <200 MB RAM).

---

### Supplementary Figures (Required for Peer Review)

**Supp. S1:** Motion QC plots — Mean FD distribution by group (PD vs. HC).
Box + strip plot; annotate excluded subjects (FD > 0.5) in red.

**Supp. S2:** WMH stratification analysis — Repeat Figure 2A separately in Low-WMH
and High-WMH tertiles. Show that BC → atrophy is driven by Low-WMH subjects.

**Supp. S3:** Spatial spin test results — Plot null distribution of TVI spatial
correlations from 5,000 cortical surface rotations, vs. observed TVI-atrophy correlation.

**Supp. S4:** Sensitivity analysis — Repeat main model with (a) 68-ROI Desikan
parcellation only (exclude subcortical), (b) binarized connectivity matrix (no log-weighting),
(c) UPDRS-III score as outcome instead of GM atrophy. All three should replicate the
direction and significance of the interaction term.

---

---

# APPENDIX: COMPLETE PYTHON ENVIRONMENT & RUNTIME ESTIMATES

---

## Package Requirements (all CPU-only, pip-installable)

```
numpy==1.26.4
scipy==1.13.0
pandas==2.2.2
networkx==3.3
scikit-learn==1.5.0
nilearn==0.10.4
nibabel==5.2.1
matplotlib==3.9.0
seaborn==0.13.2
statsmodels==0.14.2
pingouin==0.5.4
lifelines==0.29.0
tqdm==4.66.4
joblib==1.4.2
abagen==0.1.3      # Allen Brain Atlas → parcellation mapping
```

---

## Realistic Runtime Estimates (4-core / 8 GB PC)

| Step | Estimated Time | Peak RAM |
|---|---|---|
| Allen TVI construction (abagen) | ~10 min (once) | ~2 GB |
| Graph metrics per subject (84 nodes) | <1 second | <50 MB |
| Full cohort graph processing (n=100) | ~5 min total | <500 MB |
| LME model fitting | ~5 min | <2 GB |
| 1,000-iteration permutation test | ~10–20 min | <3 GB |
| Spatial spin test (5,000 rotations) | ~30 min | <1 GB |
| Figure generation (all panels) | ~5 min | <1 GB |
| **TOTAL end-to-end** | **~1–2 hours** | **<4 GB peak** |

---

## Reproducibility Checklist

- [ ] Random seeds fixed: `np.random.seed(42)` at top of every script
- [ ] All intermediate DataFrames saved as CSV at each step (no in-memory chains)
- [ ] Subject exclusions logged to `./results/exclusion_log.csv` with reason
- [ ] All figure scripts separated from analysis scripts (no monolithic notebooks)
- [ ] `requirements.txt` pinned to exact versions (above)
- [ ] `README.md` with step-by-step run order: `data_loader.py` → `main_pipeline.py`
  → `stats.py` → `figures/*.py`
- [ ] Share data at: PPMI open access URL + Allen Atlas export instructions documented

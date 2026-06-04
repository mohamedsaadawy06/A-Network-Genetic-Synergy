import numpy as np
import pandas as pd
from scipy import stats
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LinearRegression
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import fdrcorrection
from tqdm import tqdm
from pathlib import Path

# Paths
BASE_DIR = Path("C:/Users/MM/Desktop/nature isa")
RESULTS_DIR = BASE_DIR / "results"
master_csv = RESULTS_DIR / "master_features.csv"

# Load data
df = pd.read_csv(master_csv).dropna()
print(f"Loaded master features with {len(df)} rows.")

# Standardize predictors
scaler = StandardScaler()
df[["bc_z", "d_eff_z", "tvi_z", "bc_x_tvi_z"]] = scaler.fit_transform(
    df[["bc", "d_eff", "tvi", "bc_x_tvi"]]
)

# === MODEL 1: Linear Mixed Effects Model ===
print("\nFitting Linear Mixed Effects (LME) Model...")
formula = """
thickness ~ bc_z + d_eff_z + tvi_z + bc_x_tvi_z
             + age + sex + disease_dur + mean_fd
"""
lme_model = smf.mixedlm(formula, df, groups=df["subject_id"])
lme_result = lme_model.fit()
print(lme_result.summary())

# Save LME result
with open(RESULTS_DIR / "lme_summary.txt", "w") as f:
    f.write(str(lme_result.summary()))
lme_result.save(RESULTS_DIR / "lme_result.pkl")

# === MODEL 2: Fast OLS Permutation Test (1,000 iterations) ===
print("\nRunning OLS Permutation Test...")
def get_ols_betas(df_shuffled: pd.DataFrame) -> tuple[float, float]:
    """Fast OLS for permutation null distribution."""
    X = df_shuffled[["bc_z", "d_eff_z", "tvi_z", "bc_x_tvi_z",
                     "age", "sex", "disease_dur", "mean_fd"]].values
    y = df_shuffled["thickness"].values
    reg = LinearRegression().fit(X, y)
    return reg.coef_[0], reg.coef_[3]  # beta_BC, beta_BCxTVI

observed_beta_bc, observed_beta_inter = get_ols_betas(df)

null_beta_bc = []
null_beta_inter = []

for _ in tqdm(range(1000)):
    df_shuffled = df.copy()
    # Shuffle outcome within each subject to preserve regional autocorrelation
    df_shuffled["thickness"] = df.groupby("subject_id")["thickness"].transform(
        lambda x: x.sample(frac=1).values
    )
    b_bc, b_int = get_ols_betas(df_shuffled)
    null_beta_bc.append(b_bc)
    null_beta_inter.append(b_int)

# Empirical p-values
p_bc = np.mean(np.abs(null_beta_bc) >= np.abs(observed_beta_bc))
p_inter = np.mean(np.abs(null_beta_inter) >= np.abs(observed_beta_inter))

print(f"\nPermutation Test Results:")
print(f"Observed beta_BC = {observed_beta_bc:.4f}, empirical p-value = {p_bc:.4f}")
print(f"Observed beta_BCxTVI = {observed_beta_inter:.4f}, empirical p-value = {p_inter:.4f}")

# Save permutation results
np.savez(RESULTS_DIR / "permutation_results.npz",
         null_beta_bc=null_beta_bc,
         null_beta_inter=null_beta_inter,
         observed_beta_bc=observed_beta_bc,
         observed_beta_inter=observed_beta_inter,
         p_bc=p_bc,
         p_inter=p_inter)

# === MODEL 3: Group-level and Region-level Correlations ===
print("\nRunning group-level and region-level correlation analyses...")

# 1. Group-level correlation (average regional thickness vs BC across the 68 regions)
df_roi = df[['roi', 'bc', 'tvi', 'd_eff']].drop_duplicates().set_index('roi')
mean_thickness = df.groupby('roi')['thickness'].mean()
df_roi['mean_thickness'] = mean_thickness

r_group, p_group = stats.spearmanr(df_roi['bc'], df_roi['mean_thickness'])
print(f"\nGroup-level correlation across the 68 regions:")
print(f"Spearman's rho (BC vs Mean Thickness) = {r_group:.4f}, p-value = {p_group:.4e}")

# Save group-level results
df_roi.to_csv(RESULTS_DIR / "group_level_roi_metrics.csv")

# 2. Region-level correlation (regional thickness vs disease duration across subjects)
roi_pvals = {}
roi_rhos = {}
for roi in df["roi"].unique():
    sub = df[df["roi"] == roi]
    if len(sub) < 10:
        continue
    r, pval = stats.spearmanr(sub["disease_dur"], sub["thickness"])
    roi_rhos[roi] = r
    roi_pvals[roi] = pval

rois = list(roi_pvals.keys())
pvals_array = np.array([roi_pvals[r] for r in rois])
rhos_array = np.array([roi_rhos[r] for r in rois])

rejected, pvals_corrected = fdrcorrection(pvals_array, alpha=0.05, method="indep")

fdr_df = pd.DataFrame({
    "roi": rois,
    "spearman_rho": rhos_array,
    "raw_pval": pvals_array,
    "fdr_pval": pvals_corrected,
    "significant": rejected
})
fdr_df.to_csv(RESULTS_DIR / "fdr_roi_results.csv", index=False)
print(f"\nFDR Significant Regions (Thickness vs Disease Duration): {sum(rejected)}")
if sum(rejected) > 0:
    print(fdr_df[fdr_df["significant"]].sort_values("fdr_pval").to_string(index=False))

# === MODEL 4: Biophysical Generative Model Validation ===
print("\nFitting Biophysical Generative Model LME...")
df_sim = pd.read_csv(RESULTS_DIR / "biophysical_atrophy.csv").set_index("roi")
df = df.join(df_sim, on="roi", how="inner")

# Standardize simulated atrophy
df["sim_atrophy_z"] = scaler.fit_transform(df[["sim_atrophy"]])

formula_sim = """
thickness ~ sim_atrophy_z + age + sex + disease_dur + mean_fd
"""
lme_model_sim = smf.mixedlm(formula_sim, df, groups=df["subject_id"])
lme_result_sim = lme_model_sim.fit()
print(lme_result_sim.summary())

# Save summary
with open(RESULTS_DIR / "lme_biophysical_summary.txt", "w") as f:
    f.write(str(lme_result_sim.summary()))

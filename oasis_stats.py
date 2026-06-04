import numpy as np
import pandas as pd
from scipy import stats
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LinearRegression
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import fdrcorrection
from tqdm import tqdm
import os

# Paths
BASE_DIR = "C:/Users/MM/Desktop/nature isa"
RESULTS_DIR = os.path.join(BASE_DIR, "results")

def analyze_cohort(csv_filename: str, summary_txt_name: str, npz_name: str, cohort_label: str):
    csv_path = os.path.join(RESULTS_DIR, csv_filename)
    df = pd.read_csv(csv_path).dropna()
    print(f"\n=================== ANALYZING COHORT: {cohort_label} ===================")
    print(f"Loaded master features with {len(df)} rows ({len(df) // 68} subjects).")

    # Standardize predictors
    scaler = StandardScaler()
    df[["bc_z", "d_eff_z", "tvi_z", "bc_x_tvi_z"]] = scaler.fit_transform(
        df[["bc", "d_eff", "tvi", "bc_x_tvi"]]
    )

    # 1. Fit Linear Mixed Effects Model (LME)
    print(f"Fitting Linear Mixed Effects (LME) Model for {cohort_label}...")
    if cohort_label == "CN":
        formula = """
        atrophy_rate ~ bc_z + d_eff_z + tvi_z + bc_x_tvi_z
                     + age + sex + mean_fd
        """
    else:
        formula = """
        atrophy_rate ~ bc_z + d_eff_z + tvi_z + bc_x_tvi_z
                     + age + sex + disease_dur + mean_fd
        """
    lme_model = smf.mixedlm(formula, df, groups=df["subject_id"])
    lme_result = lme_model.fit()
    print(lme_result.summary())

    # Save LME result
    summary_txt_path = os.path.join(RESULTS_DIR, summary_txt_name)
    with open(summary_txt_path, "w") as f:
        f.write(str(lme_result.summary()))

    # 2. Run OLS Permutation Test
    print(f"Running OLS Permutation Test (1,000 iterations) for {cohort_label}...")
    def get_ols_betas(df_shuffled: pd.DataFrame) -> tuple[float, float]:
        if cohort_label == "CN":
            X = df_shuffled[["bc_z", "d_eff_z", "tvi_z", "bc_x_tvi_z",
                             "age", "sex", "mean_fd"]].values
        else:
            X = df_shuffled[["bc_z", "d_eff_z", "tvi_z", "bc_x_tvi_z",
                             "age", "sex", "disease_dur", "mean_fd"]].values
        y = df_shuffled["atrophy_rate"].values
        reg = LinearRegression().fit(X, y)
        return reg.coef_[0], reg.coef_[3]  # beta_BC, beta_BCxTVI

    observed_beta_bc, observed_beta_inter = get_ols_betas(df)

    null_beta_bc = []
    null_beta_inter = []

    for _ in tqdm(range(1000), desc=f"Permutations for {cohort_label}"):
        df_shuffled = df.copy()
        # Shuffle outcome within each subject to preserve regional autocorrelation
        df_shuffled["atrophy_rate"] = df.groupby("subject_id")["atrophy_rate"].transform(
            lambda x: x.sample(frac=1).values
        )
        b_bc, b_int = get_ols_betas(df_shuffled)
        null_beta_bc.append(b_bc)
        null_beta_inter.append(b_int)

    # Empirical p-values
    p_bc = np.mean(np.abs(null_beta_bc) >= np.abs(observed_beta_bc))
    p_inter = np.mean(np.abs(null_beta_inter) >= np.abs(observed_beta_inter))

    print(f"\nPermutation Test Results for {cohort_label}:")
    print(f"Observed beta_BC = {observed_beta_bc:.4f}, empirical p-value = {p_bc:.4f}")
    print(f"Observed beta_BCxTVI = {observed_beta_inter:.4f}, empirical p-value = {p_inter:.4f}")

    # Save permutation results
    npz_path = os.path.join(RESULTS_DIR, npz_name)
    np.savez(npz_path,
             null_beta_bc=null_beta_bc,
             null_beta_inter=null_beta_inter,
             observed_beta_bc=observed_beta_bc,
             observed_beta_inter=observed_beta_inter,
             p_bc=p_bc,
             p_inter=p_inter)

    # 3. Group-level correlation (average regional atrophy vs BC across the 68 regions)
    df_roi = df[['roi', 'bc', 'tvi', 'd_eff']].drop_duplicates().set_index('roi')
    mean_atrophy = df.groupby('roi')['atrophy_rate'].mean()
    df_roi['mean_atrophy'] = mean_atrophy
    
    r_group, p_group = stats.spearmanr(df_roi['bc'], df_roi['mean_atrophy'])
    print(f"\nGroup-level correlation across the 68 regions for {cohort_label}:")
    print(f"Spearman's rho (BC vs Mean Atrophy) = {r_group:.4f}, p-value = {p_group:.4e}")
    
    # Save group-level results
    roi_csv_path = os.path.join(RESULTS_DIR, f"oasis_{cohort_label.lower()}_group_level_roi_metrics.csv")
    df_roi.to_csv(roi_csv_path)

def main():
    # Analyze AD cohort
    analyze_cohort(
        csv_filename="oasis_master_features.csv",
        summary_txt_name="oasis_lme_summary.txt",
        npz_name="oasis_permutation_results.npz",
        cohort_label="AD"
    )

    # Analyze CN cohort
    analyze_cohort(
        csv_filename="oasis_cn_master_features.csv",
        summary_txt_name="oasis_cn_lme_summary.txt",
        npz_name="oasis_cn_permutation_results.npz",
        cohort_label="CN"
    )

if __name__ == "__main__":
    main()

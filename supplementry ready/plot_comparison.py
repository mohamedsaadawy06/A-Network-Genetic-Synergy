import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
import statsmodels.api as sm

# Set Nature plot style
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['axes.linewidth'] = 0.8
plt.rcParams['xtick.major.width'] = 0.8
plt.rcParams['ytick.major.width'] = 0.8
plt.rcParams['xtick.labelsize'] = 8
plt.rcParams['ytick.labelsize'] = 8
plt.rcParams['axes.labelsize'] = 9
plt.rcParams['axes.titlesize'] = 10
plt.rcParams['legend.fontsize'] = 7

# Paths
BASE_DIR = "C:/Users/MM/Desktop/nature isa"
RESULTS_DIR = os.path.join(BASE_DIR, "results")

def get_cohort_stats(csv_path: str) -> tuple[pd.Series, pd.Series]:
    """Fit OLS on standardized predictors and return coefficients and standard errors."""
    df = pd.read_csv(csv_path).dropna()
    scaler = StandardScaler()
    
    # Ensure columns exist and compute interaction if missing
    if "bc_x_tvi" not in df.columns:
        df["bc_x_tvi"] = df["bc"] * df["tvi"]
        
    df[["bc_z", "d_eff_z", "tvi_z", "bc_x_tvi_z"]] = scaler.fit_transform(
        df[["bc", "d_eff", "tvi", "bc_x_tvi"]]
    )
    
    covariates = ["age", "sex", "mean_fd"]
    if "disease_dur" in df.columns and df["disease_dur"].nunique() > 1:
        covariates.append("disease_dur")
        
    X = df[["bc_z", "d_eff_z", "tvi_z", "bc_x_tvi_z"] + covariates]
    X = sm.add_constant(X)
    if "thickness" in df.columns:
        y = df["thickness"]
    else:
        y = df["atrophy_rate"]
    
    model = sm.OLS(y, X)
    results = model.fit()
    
    predictors = ["bc_z", "d_eff_z", "tvi_z", "bc_x_tvi_z"]
    return results.params[predictors], results.bse[predictors]

def main():
    print("Loading cohort datasets...")
    pd_csv = os.path.join(RESULTS_DIR, "master_features.csv")
    oasis_ad_csv = os.path.join(RESULTS_DIR, "oasis_master_features.csv")
    oasis_cn_csv = os.path.join(RESULTS_DIR, "oasis_cn_master_features.csv")
    adni_ad_csv = os.path.join(RESULTS_DIR, "adni_master_features.csv")
    adni_cn_csv = os.path.join(RESULTS_DIR, "adni_cn_master_features.csv")
    
    # Check if files exist
    paths = {
        "PPMI PD": pd_csv,
        "OASIS AD": oasis_ad_csv,
        "OASIS CN": oasis_cn_csv,
        "ADNI AD": adni_ad_csv,
        "ADNI CN": adni_cn_csv
    }
    
    for name, path in paths.items():
        if not os.path.exists(path):
            print(f"Error: Required file {path} for {name} does not exist. Run pipelines first.")
            return
            
    pd_beta, pd_se = get_cohort_stats(pd_csv)
    oasis_ad_beta, oasis_ad_se = get_cohort_stats(oasis_ad_csv)
    oasis_cn_beta, oasis_cn_se = get_cohort_stats(oasis_cn_csv)
    adni_ad_beta, adni_ad_se = get_cohort_stats(adni_ad_csv)
    adni_cn_beta, adni_cn_se = get_cohort_stats(adni_cn_csv)
    
    print("\nPD (PPMI) Coefficients:")
    print(pd_beta)
    print("\nAD (ADNI) Coefficients:")
    print(adni_ad_beta)
    print("\nAD (OASIS 3) Coefficients:")
    print(oasis_ad_beta)
    print("\nCN (ADNI) Coefficients:")
    print(adni_cn_beta)
    print("\nCN (OASIS 3) Coefficients:")
    print(oasis_cn_beta)
    
    # ----------------------------------------------------
    # Plot Comparative Forest Plot
    # ----------------------------------------------------
    fig, ax = plt.subplots(figsize=(7, 5), dpi=300)
    
    predictors_labels = [
        "Betweenness Centrality (BC)",
        "Effective Distance (D_eff)",
        "Transcriptomic Vulnerability (TVI)",
        "BC x TVI (Interaction)"
    ]
    
    y_positions = np.arange(len(predictors_labels))
    width = 0.15  # Spacing between cohorts
    
    # Plot ADNI AD
    ax.errorbar(
        adni_ad_beta.values, y_positions + 2 * width, xerr=1.96 * adni_ad_se.values,
        fmt="o", color="#800000", capsize=3, elinewidth=1.2, markersize=5,
        label="AD Dementia (ADNI)"
    )
    
    # Plot OASIS AD
    ax.errorbar(
        oasis_ad_beta.values, y_positions + width, xerr=1.96 * oasis_ad_se.values,
        fmt="o", color="#d62728", capsize=3, elinewidth=1.2, markersize=5,
        label="AD Dementia (OASIS 3)"
    )
    
    # Plot PD (PPMI)
    ax.errorbar(
        pd_beta.values, y_positions, xerr=1.96 * pd_se.values,
        fmt="o", color="#1f77b4", capsize=3, elinewidth=1.2, markersize=5,
        label="Parkinson's Disease (PPMI)"
    )
    
    # Plot ADNI CN
    ax.errorbar(
        adni_cn_beta.values, y_positions - width, xerr=1.96 * adni_cn_se.values,
        fmt="o", color="#7f7f7f", capsize=3, elinewidth=1.2, markersize=5,
        label="Cognitively Normal (ADNI)"
    )
    
    # Plot OASIS Control
    ax.errorbar(
        oasis_cn_beta.values, y_positions - 2 * width, xerr=1.96 * oasis_cn_se.values,
        fmt="o", color="#2c3e50", capsize=3, elinewidth=1.2, markersize=5,
        label="Cognitively Normal (OASIS 3)"
    )
    
    # Reference line at beta=0
    ax.axvline(0, color="gray", linestyle="--", linewidth=0.8, zorder=0)
    
    ax.set_yticks(y_positions)
    ax.set_yticklabels(predictors_labels, fontsize=8)
    ax.set_xlabel("Standardized Effect Coefficient (Beta)", fontsize=8)
    ax.set_title("Cross-Disease Cohort Validation Comparison", weight="bold", pad=12)
    ax.legend(loc="upper right", frameon=True, facecolor="white", edgecolor="none")
    ax.grid(True, linestyle="--", alpha=0.3)
    
    # Invert y-axis to read top-to-bottom
    ax.invert_yaxis()
    
    plt.tight_layout()
    plot_path = os.path.join(RESULTS_DIR, "Figure_6_Comparative_Validation.png")
    fig.savefig(plot_path, bbox_inches="tight")
    plt.close(fig)
    print(f"\nSaved comparison forest plot to: {plot_path}")

if __name__ == "__main__":
    main()

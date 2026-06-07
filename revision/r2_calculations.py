import pandas as pd
import numpy as np
import statsmodels.formula.api as smf
from sklearn.preprocessing import StandardScaler
from pathlib import Path

def calculate_nakagawa_r2(model_result, df, formula):
    """
    Calculate Marginal and Conditional R^2 for LME models.
    Marginal R^2: Proportion of variance explained by fixed effects.
    Conditional R^2: Proportion of variance explained by fixed + random effects.
    """
    # 1. Variance of fixed effects predictions
    # We predict the response using only fixed effects
    fixed_effects = model_result.params.drop("Group Var")
    X = model_result.model.exog
    y_pred_fixed = np.dot(X, fixed_effects)
    var_fixed = np.var(y_pred_fixed)
    
    # 2. Variance of random effects
    var_random = model_result.cov_re.iloc[0,0]
    
    # 3. Residual variance
    var_resid = model_result.scale
    
    # 4. Compute R^2
    total_var = var_fixed + var_random + var_resid
    marginal_r2 = var_fixed / total_var
    conditional_r2 = (var_fixed + var_random) / total_var
    
    # Also calculate incremental variance of the interaction term specifically
    # Refit model without bc_x_tvi_z
    formula_reduced = formula.replace(" + bc_x_tvi_z", "").replace("bc_x_tvi_z + ", "")
    model_reduced = smf.mixedlm(formula_reduced, df, groups=df["subject_id"]).fit(method='cg', maxiter=100)
    
    fixed_effects_red = model_reduced.params.drop("Group Var")
    X_red = model_reduced.model.exog
    y_pred_fixed_red = np.dot(X_red, fixed_effects_red)
    var_fixed_red = np.var(y_pred_fixed_red)
    
    var_random_red = model_reduced.cov_re.iloc[0,0]
    var_resid_red = model_reduced.scale
    total_var_red = var_fixed_red + var_random_red + var_resid_red
    marginal_r2_red = var_fixed_red / total_var_red
    
    incremental_r2 = marginal_r2 - marginal_r2_red
    
    return marginal_r2, conditional_r2, incremental_r2

def process_cohort(name, csv_path, formula, y_col):
    df = pd.read_csv(csv_path).dropna()
    
    if name == "PPMI":
        cols_to_scale = ["bc", "d_eff", "tvi", "bc_x_tvi", "pvs_volume", "bc_x_pvs"]
        new_cols = ["bc_z", "d_eff_z", "tvi_z", "bc_x_tvi_z", "pvs_z", "bc_x_pvs_z"]
    elif "ADNI" in name:
        cols_to_scale = ["bc", "tvi", "bc_x_tvi"]
        new_cols = ["bc_z", "tvi_z", "bc_x_tvi_z"]
    else: # OASIS
        cols_to_scale = ["bc", "d_eff", "tvi", "bc_x_tvi"]
        new_cols = ["bc_z", "d_eff_z", "tvi_z", "bc_x_tvi_z"]
        
    scaler = StandardScaler()
    df[new_cols] = scaler.fit_transform(df[cols_to_scale])
    
    model = smf.mixedlm(formula, df, groups=df["subject_id"])
    result = model.fit(method='cg') # cg is more robust for these models
    
    m_r2, c_r2, inc_r2 = calculate_nakagawa_r2(result, df, formula)
    return {
        "Cohort": name,
        "Marginal R2": m_r2,
        "Conditional R2": c_r2,
        "Incremental R2 (Interaction)": inc_r2
    }

def main():
    base_dir = Path("C:/Users/MM/Desktop/nature isa/results")
    
    configs = [
        ("PPMI", base_dir / "master_features.csv", "thickness ~ bc_z + d_eff_z + tvi_z + bc_x_tvi_z + pvs_z + bc_x_pvs_z + age + sex + disease_dur + mean_fd", "thickness"),
        ("ADNI AD", base_dir / "adni_master_features.csv", "atrophy_rate ~ bc_z + tvi_z + bc_x_tvi_z + age + sex + disease_dur + mean_fd + icv", "atrophy_rate"),
        ("ADNI CN", base_dir / "adni_cn_master_features.csv", "atrophy_rate ~ bc_z + tvi_z + bc_x_tvi_z + age + sex + mean_fd + icv", "atrophy_rate"),
        ("OASIS-3 AD", base_dir / "oasis_master_features.csv", "atrophy_rate ~ bc_z + d_eff_z + tvi_z + bc_x_tvi_z + age + sex + disease_dur + mean_fd", "atrophy_rate"),
        ("OASIS-3 CN", base_dir / "oasis_cn_master_features.csv", "atrophy_rate ~ bc_z + d_eff_z + tvi_z + bc_x_tvi_z + age + sex + mean_fd", "atrophy_rate")
    ]
    
    results = []
    for name, path, form, y_col in configs:
        if not path.exists():
            print(f"Skipping {name}, file not found.")
            continue
        print(f"Processing {name}...")
        results.append(process_cohort(name, path, form, y_col))
        
    df_res = pd.DataFrame(results)
    print("\n--- Model R^2 Summary ---")
    print(df_res.to_string(index=False))
    df_res.to_csv(base_dir / "r2_summary.csv", index=False)

if __name__ == "__main__":
    main()

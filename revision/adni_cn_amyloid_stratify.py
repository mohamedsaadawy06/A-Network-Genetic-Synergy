import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from sklearn.preprocessing import StandardScaler
from pathlib import Path

def main():
    base_dir = Path("C:/Users/MM/Desktop/nature isa")
    
    # 1. Load ADNI CN features
    cn_csv = base_dir / "results" / "adni_cn_master_features.csv"
    df_cn = pd.read_csv(cn_csv)
    
    # 2. Load Amyloid PET PHC
    amy_csv = base_dir / "data" / "ADSP_PHC_PET_Amyloid_Simple_07Jun2026.csv"
    df_amy = pd.read_csv(amy_csv)
    
    # Format PTID to match subject_id
    # Wait, PTID in ADNI is exactly the same format: XXX_S_XXXX
    # Let's ensure types are strings
    df_amy['PTID'] = df_amy['PTID'].astype(str)
    df_cn['subject_id'] = df_cn['subject_id'].astype(str)
    
    # Keep the max amyloid status per subject (if they ever turned positive, consider them positive)
    # Alternatively, just merge on PTID directly, taking the first valid amyloid status
    # Let's group by PTID and take max of PHC_AMYLOID_STATUS
    amy_status = df_amy.dropna(subset=['PHC_AMYLOID_STATUS']).groupby('PTID')['PHC_AMYLOID_STATUS'].max().reset_index()
    
    # Merge
    df_merged = pd.merge(df_cn, amy_status, left_on='subject_id', right_on='PTID', how='inner')
    
    print(f"Total ADNI CN rows: {len(df_cn)}")
    print(f"Matched rows with Amyloid status: {len(df_merged)}")
    print(f"Matched subjects: {df_merged['subject_id'].nunique()} out of {df_cn['subject_id'].nunique()}")
    
    if len(df_merged) == 0:
        print("No matching subjects found!")
        return
        
    df_merged = df_merged.dropna(subset=['atrophy_rate', 'bc', 'tvi', 'age', 'sex', 'mean_fd', 'icv'])
    
    # Standardize predictors (fit on the whole matched dataset for consistency)
    scaler = StandardScaler()
    df_merged[["bc_z", "tvi_z", "bc_x_tvi_z"]] = scaler.fit_transform(
        df_merged[["bc", "tvi", "bc_x_tvi"]]
    )
    
    # Split
    df_pos = df_merged[df_merged['PHC_AMYLOID_STATUS'] == 1.0].copy()
    df_neg = df_merged[df_merged['PHC_AMYLOID_STATUS'] == 0.0].copy()
    
    print(f"Amyloid Positive Subjects: {df_pos['subject_id'].nunique()}")
    print(f"Amyloid Negative Subjects: {df_neg['subject_id'].nunique()}")
    
    formula = "atrophy_rate ~ bc_z + tvi_z + bc_x_tvi_z + age + sex + mean_fd + icv"
    
    print("\n--- Amyloid Positive Cohort ---")
    if len(df_pos) > 0:
        model_pos = smf.mixedlm(formula, df_pos, groups=df_pos["subject_id"])
        res_pos = model_pos.fit(method='cg') # 'cg' or 'lbfgs' sometimes needed for convergence, default is 'bfgs'
        print(res_pos.summary())
    else:
        print("Not enough subjects.")
        
    print("\n--- Amyloid Negative Cohort ---")
    if len(df_neg) > 0:
        model_neg = smf.mixedlm(formula, df_neg, groups=df_neg["subject_id"])
        res_neg = model_neg.fit(method='cg')
        print(res_neg.summary())
    else:
        print("Not enough subjects.")

if __name__ == "__main__":
    main()

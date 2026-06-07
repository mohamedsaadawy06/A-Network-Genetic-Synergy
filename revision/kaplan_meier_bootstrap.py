import numpy as np
import pandas as pd
from lifelines.statistics import logrank_test
from pathlib import Path
from tqdm import tqdm

def main():
    base_dir = Path("C:/Users/MM/Desktop/nature isa")
    df = pd.read_csv(base_dir / "results" / "master_features.csv").dropna()
    df_roi_metrics = pd.read_csv(base_dir / "results" / "group_level_roi_metrics.csv")
    
    # Subject stratification (from figures.py)
    tvi_vals = df_roi_metrics["tvi"].values
    bc_vals = df_roi_metrics["bc"].values
    bc_x_tvi = bc_vals * tvi_vals
    top_5_rois = df_roi_metrics.iloc[np.argsort(bc_x_tvi)[-5:]]["roi"].values
    
    df_top_5 = df[df["roi"].isin(top_5_rois)]
    subject_vulnerability = df_top_5.groupby("subject_id")["thickness"].mean()
    median_vulnerability = subject_vulnerability.median()
    vulnerable_subjects = subject_vulnerability[subject_vulnerability < median_vulnerability].index.astype(str)
    resilient_subjects = subject_vulnerability[subject_vulnerability >= median_vulnerability].index.astype(str)

    # 1. MoCA Event Data Construction (from figures.py)
    df_xl = pd.read_excel(base_dir / "data_ppmi/PPMI_Curated_Data_Cut_Public_20260319.xlsx", sheet_name="20260316")
    df_xl["PATNO_str"] = df_xl["PATNO"].astype(str)
    
    cognitive_events = []
    for pat in subject_vulnerability.index.astype(str):
        subj_rows = df_xl[(df_xl["PATNO_str"] == pat) & df_xl["EVENT_ID"].isin(["BL", "V04", "V06"])].sort_values("EVENT_ID")
        if len(subj_rows) < 2:
            continue
        bl_moca = subj_rows[subj_rows["EVENT_ID"] == "BL"]["moca"].values
        if len(bl_moca) == 0 or pd.isna(bl_moca[0]):
            continue
        bl_moca = bl_moca[0]
        
        event_observed = 0
        duration_months = 24.0
        for _, row in subj_rows.iterrows():
            visit = row["EVENT_ID"]
            if visit == "BL":
                continue
            m = row["moca"]
            if not pd.isna(m) and (bl_moca - m) >= 3:
                event_observed = 1
                duration_months = 12.0 if visit == "V04" else 24.0
                break
                
        cognitive_events.append({
            "subject_id": pat,
            "event": event_observed,
            "duration": duration_months,
            "group": "Vulnerable Hubs" if pat in vulnerable_subjects else "Resilient Hubs"
        })
        
    df_surv = pd.DataFrame(cognitive_events)
    
    g1 = df_surv[df_surv["group"] == "Vulnerable Hubs"]
    g2 = df_surv[df_surv["group"] == "Resilient Hubs"]
    
    # Original uncorrected Log-Rank
    lr_original = logrank_test(g1["duration"], g2["duration"], g1["event"], g2["event"])
    
    print("--- Kaplan-Meier MoCA Analysis ---")
    print(f"Original Log-rank stat: {lr_original.test_statistic:.4f}, p = {lr_original.p_value:.4f}")
    
    # Bootstrap Log-rank
    n_boot = 10000
    boot_stats = []
    
    np.random.seed(42)
    for _ in tqdm(range(n_boot), desc="Bootstrapping Log-Rank"):
        # Resample individuals with replacement from the entire cohort
        df_boot = df_surv.sample(frac=1.0, replace=True)
        g1_boot = df_boot[df_boot["group"] == "Vulnerable Hubs"]
        g2_boot = df_boot[df_boot["group"] == "Resilient Hubs"]
        
        # Only compute if we have variance in both
        if len(g1_boot) > 0 and len(g2_boot) > 0 and g1_boot["event"].sum() > 0 or g2_boot["event"].sum() > 0:
            lr_boot = logrank_test(g1_boot["duration"], g2_boot["duration"], g1_boot["event"], g2_boot["event"])
            boot_stats.append(lr_boot.test_statistic)
        else:
            boot_stats.append(0.0)
            
    boot_stats = np.array(boot_stats)
    ci_lower = np.percentile(boot_stats, 2.5)
    ci_upper = np.percentile(boot_stats, 97.5)
    
    print(f"Bootstrap 95% CI of log-rank test statistic: [{ci_lower:.4f}, {ci_upper:.4f}]")
    
    # Calculate bootstrap p-value
    p_boot = np.mean(boot_stats >= lr_original.test_statistic)
    print(f"Bootstrap empirical p-value: {p_boot:.4f}")
    
    # Bonferroni correction for 2 endpoints (MoCA and UPDRS-III)
    p_bonf = min(1.0, p_boot * 2)
    print(f"Bonferroni-corrected p-value (n=2): {p_bonf:.4f}")

if __name__ == "__main__":
    main()

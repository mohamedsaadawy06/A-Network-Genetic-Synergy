import pandas as pd
import numpy as np
import os
from pathlib import Path
from tqdm import tqdm

from config import DATA_ROOT, ALLEN_PATH, RESULTS_DIR, SEED_REGION, MOTION_FD_THRESH
from data_loader import (load_connectivity_matrix, load_real_cth,
                         load_motion_fd, get_subject_list, ROI_LABELS)
from graph_metrics import (build_weighted_graph, compute_betweenness_centrality,
                            compute_effective_distance)
from tvi import build_tvi

# Pre-load TVI and AQP4 once
tvi_series, aqp4_series = build_tvi(ALLEN_PATH, ROI_LABELS)

def process_subject(subject_id: str, clinical_row: pd.Series, cth_df: pd.DataFrame, exclusion_log: list) -> list | None:
    """
    Extract features and baseline cortical thickness for a single subject.
    """
    # 1. Motion quality control
    fd = load_motion_fd(subject_id, "BL")
    if fd > MOTION_FD_THRESH:
        exclusion_log.append({"PATNO": subject_id, "Reason": f"High motion: mean FD = {fd:.3f} > {MOTION_FD_THRESH}"})
        return None

    # 2. Load baseline connectivity matrix
    mat_bl = load_connectivity_matrix(subject_id, "BL")
    if mat_bl is None:
        exclusion_log.append({"PATNO": subject_id, "Reason": "Missing DTI connectivity matrix"})
        return None

    # 3. Build graph and compute topology metrics
    G = build_weighted_graph(mat_bl, ROI_LABELS)
    bc_dict = compute_betweenness_centrality(G)
    deff_dict = compute_effective_distance(mat_bl, ROI_LABELS, SEED_REGION)

    # 4. Load real cortical thickness (Baseline BL)
    cth_series = load_real_cth(subject_id, cth_df)
    if cth_series is None:
        exclusion_log.append({"PATNO": subject_id, "Reason": "Missing baseline cortical thickness"})
        return None

    pvs_dict = {}
    pvs_file = DATA_ROOT / f"sub-{subject_id}" / "ses-BL" / "anat" / f"sub-{subject_id}_PVS_volumes.csv"
    if os.path.exists(pvs_file):
        df_pvs = pd.read_csv(pvs_file)
        for _, r in df_pvs.iterrows():
            pvs_dict[r['roi']] = float(r['pvs_volume'])

    # 5. Assemble records
    records = []
    for roi in ROI_LABELS:
        if roi not in bc_dict:
            continue
        records.append({
            "subject_id"  : subject_id,
            "roi"         : roi,
            "bc"          : bc_dict.get(roi, 0.0),
            "d_eff"       : deff_dict.get(roi, np.nan),
            "tvi"         : tvi_series.get(roi, 0.0),
            "bc_x_tvi"    : bc_dict.get(roi, 0.0) * tvi_series.get(roi, 0.0),
            "aqp4"        : aqp4_series.get(roi, 0.0),
            "bc_x_aqp4"   : bc_dict.get(roi, 0.0) * aqp4_series.get(roi, 0.0),
            "pvs_volume"  : pvs_dict.get(roi, np.nan),
            "bc_x_pvs"    : bc_dict.get(roi, 0.0) * pvs_dict.get(roi, np.nan),
            "thickness"   : cth_series.get(roi, np.nan),
            "age"         : clinical_row["AGE"],
            "sex"         : clinical_row["SEX"],
            "disease_dur" : clinical_row["DISDUR"],
            "mean_fd"     : fd,
            "updrs3"      : clinical_row["UPDRS3"],
            "moca"        : clinical_row["MOCA"]
        })

    return records

def main():
    clinical_csv = DATA_ROOT / "PPMI_Clinical_Data.csv"
    clinical_df = get_subject_list(clinical_csv)
    
    # Load CTH CSV from workspace parent
    cth_csv = DATA_ROOT.parent / "FS7_APARC_CTH_04Jun2026.csv"
    if not cth_csv.exists():
        print(f"Error: {cth_csv} does not exist!")
        return
    cth_df = pd.read_csv(cth_csv)
    
    exclusion_log = []
    all_records = []
    
    print("Processing cohort subjects...")
    for _, row in tqdm(clinical_df.iterrows(), total=len(clinical_df)):
        subject_id = str(int(row["PATNO"]))
        result = process_subject(subject_id, row, cth_df, exclusion_log)
        if result:
            all_records.extend(result)
            
    # Save exclusion log
    pd.DataFrame(exclusion_log).to_csv(RESULTS_DIR / "exclusion_log.csv", index=False)
    print(f"Logged {len(exclusion_log)} subject exclusions to results/exclusion_log.csv")
    
    # Save master features
    master_df = pd.DataFrame(all_records)
    master_df.to_csv(RESULTS_DIR / "master_features.csv", index=False)
    print(f"Saved {len(master_df)} ROI-subject rows to results/master_features.csv")

if __name__ == "__main__":
    main()

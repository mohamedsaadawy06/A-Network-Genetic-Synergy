import pandas as pd
import numpy as np
from pathlib import Path
from config import DATA_ROOT, MOTION_FD_THRESH, COORDS_PATH

def load_roi_labels() -> list:
    """Load the 68 region labels in order from the coordinates file."""
    labels = []
    with open(COORDS_PATH, "r") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                labels.append(f"{parts[0]}_{parts[1]}")
    return labels

ROI_LABELS = load_roi_labels()

def load_connectivity_matrix(subject_id: str, session: str = "BL") -> np.ndarray | None:
    """
    Load the 68x68 parcellated connectivity matrix.
    Symmetrizes it and removes self-loops.
    """
    mat_path = DATA_ROOT / f"sub-{subject_id}" / f"ses-{session}" / \
               "dwi" / f"sub-{subject_id}_connectivity_matrix.csv"

    if not mat_path.exists():
        return None

    mat = pd.read_csv(mat_path, header=None).values.astype(np.float32)

    # Symmetrize and remove self-loops
    mat = (mat + mat.T) / 2.0
    np.fill_diagonal(mat, 0)

    return mat

def load_gm_volumes(subject_id: str, session: str) -> pd.Series | None:
    """Load regional GM volumes (68 ROIs) as a pandas Series."""
    vol_path = DATA_ROOT / f"sub-{subject_id}" / f"ses-{session}" / \
               "anat" / f"sub-{subject_id}_GM_volumes.csv"
    if not vol_path.exists():
        return None
    return pd.read_csv(vol_path, index_col=0).squeeze()

def load_real_cth(subject_id: str, cth_df: pd.DataFrame) -> pd.Series | None:
    """Extract baseline CTH for a subject and map to ROI labels."""
    sub_df = cth_df[(cth_df["PATNO"] == int(subject_id)) & (cth_df["EVENT_ID"] == "BL")]
    if sub_df.empty:
        return None
    row = sub_df.iloc[0]
    
    # Map to the 68 ROI labels
    vals = {}
    for roi in ROI_LABELS:
        hemi = "lh" if roi.startswith("L_") else "rh"
        reg_name = roi[2:]
        col = f"{hemi}_{reg_name}"
        if col in row.index:
            vals[roi] = float(row[col])
    return pd.Series(vals)

def load_motion_fd(subject_id: str, session: str = "BL") -> float:
    """Return mean framewise displacement for motion quality control."""
    fd_path = DATA_ROOT / f"sub-{subject_id}" / f"ses-{session}" / \
              "func" / f"sub-{subject_id}_mean_FD.txt"
    if not fd_path.exists():
        return 0.0
    with open(fd_path, "r") as f:
        val = f.read().strip()
    return float(val) if val else 0.0

def get_subject_list(clinical_csv: Path) -> pd.DataFrame:
    """Load clinical subjects from CSV."""
    df = pd.read_csv(clinical_csv)
    df = df.dropna(subset=["PATNO", "AGE", "SEX"])
    return df

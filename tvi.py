import pandas as pd
import numpy as np
from scipy.stats import zscore

def build_tvi(allen_csv_path: str, roi_labels: list):
    """
    Compute TVI(r) = z(SNCA_r) - z(GBA_r) - z(PARK2_r)
    for each ROI. Also return AQP4 series.
    """
    df = pd.read_csv(allen_csv_path, index_col=0)
    tvi = pd.Series(index=roi_labels, dtype=float)

    # Compute z-scores
    snca_z = zscore(df["SNCA"])
    gba_z = zscore(df["GBA"])
    park2_z = zscore(df["PARK2"])
    aqp4_z = zscore(df["AQP4"]) if "AQP4" in df.columns else pd.Series(0, index=df.index)

    df_z = pd.DataFrame({
        "SNCA_z": snca_z,
        "GBA_z": gba_z,
        "PARK2_z": park2_z,
        "AQP4_z": aqp4_z
    }, index=df.index)

    aqp4_series = pd.Series(index=roi_labels, dtype=float)
    for roi in roi_labels:
        if roi not in df_z.index:
            tvi[roi] = 0.0
            aqp4_series[roi] = 0.0
            continue
        row = df_z.loc[roi]
        tvi[roi] = row["SNCA_z"] - row["GBA_z"] - row["PARK2_z"]
        aqp4_series[roi] = row["AQP4_z"]

    return tvi, aqp4_series

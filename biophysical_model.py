import os
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from pathlib import Path
import networkx as nx

# Paths
BASE_DIR = Path("C:/Users/MM/Desktop/nature isa")
DATA_DIR = BASE_DIR / "data"
RESULTS_DIR = BASE_DIR / "results"
PPMI_DIR = DATA_DIR / "PPMI"

def main():
    print("Initializing Biophysical ODE Network Model...")
    
    # 1. Load Cammoun033 coordinates and label files to get labels in order
    roi_labels = []
    with open(DATA_DIR / "Cammoun033_coords.txt", "r") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                roi_labels.append(f"{parts[0]}_{parts[1]}")
                
    subdirs = [d for d in os.listdir(PPMI_DIR) if d.startswith("sub-")]
    sample_sub = subdirs[0]
    mat_path = PPMI_DIR / sample_sub / "ses-BL" / "dwi" / f"{sample_sub}_connectivity_matrix.csv"
    C = pd.read_csv(mat_path, header=None).values.astype(np.float32)
    
    # Symmetrize and remove self-loops
    C = (C + C.T) / 2.0
    np.fill_diagonal(C, 0)
    
    # Log-transform and normalize adjacency
    W = np.log1p(C)
    # Scale W to prevent extreme diffusion rates
    W = W / np.max(W) if np.max(W) > 0 else W
    
    # 3. Load gene expression data
    allen_df = pd.read_csv(DATA_DIR / "allen_expression.csv", index_col=0)
    
    # Align gene expression profiles to ROI labels
    snca_z = allen_df.loc[roi_labels, "SNCA"].values
    gba_z = allen_df.loc[roi_labels, "GBA"].values
    
    # 4. Define ODE parameters
    N = len(roi_labels)
    rho = 0.05    # Global diffusion rate
    gamma_0 = 0.1 # Base seeding rate
    d_0 = 0.03    # Base clearance rate
    
    # Regional synthesis, clearance, and normal protein levels
    gamma = gamma_0 * (1.0 + 0.2 * snca_z)
    d = d_0 * (1.0 + 0.2 * gba_z)
    y = 1.0 * (1.0 + 0.2 * snca_z)
    
    # 5. Define ODE system
    # x_i is the toxic protein concentration in region i
    def ode_system(t, x):
        dxdt = np.zeros(N)
        for i in range(N):
            # Local growth: misfolding - clearance
            local_growth = gamma[i] * x[i] * y[i] - d[i] * x[i]
            
            # Graph Laplacian diffusion: rho * sum_j W_ij * (x_j - x_i)
            diffusion = rho * np.sum(W[i, :] * (x - x[i]))
            
            dxdt[i] = local_growth + diffusion
        return dxdt
        
    # 6. Set initial conditions (Seed pathology in L_entorhinal)
    x0 = np.zeros(N)
    seed_idx = roi_labels.index("L_entorhinal")
    x0[seed_idx] = 0.1 # Seed concentration
    
    # 7. Solve the system over simulated time (T = 100)
    t_span = (0, 100)
    t_eval = np.linspace(0, 100, 200)
    print("Solving network ODE system...")
    sol = solve_ivp(ode_system, t_span, x0, t_eval=t_eval, method="RK45")
    
    # Extract concentration timeline
    # sol.y has shape (N, len(t_eval))
    
    # 8. Calculate simulated atrophy (cumulative toxic burden)
    # Integral of x_i(t) over time
    import scipy.integrate as integrate
    sim_atrophy = integrate.trapezoid(sol.y, sol.t, axis=1)
    
    # Normalize simulated atrophy to percent-like scale for correlation checks
    sim_atrophy = (sim_atrophy - sim_atrophy.min()) / (sim_atrophy.max() - sim_atrophy.min()) * -2.0
    
    # Save simulation results
    df_sim = pd.DataFrame({
        "roi": roi_labels,
        "sim_atrophy": sim_atrophy
    })
    df_sim.to_csv(RESULTS_DIR / "biophysical_atrophy.csv", index=False)
    print("Saved results to results/biophysical_atrophy.csv")
    
    # 9. Verify correlation with observed patient thickness
    df_master = pd.read_csv(RESULTS_DIR / "master_features.csv")
    mean_thickness = df_master.groupby("roi")["thickness"].mean()
    
    df_compare = df_sim.set_index("roi").join(mean_thickness, how="inner")
    
    from scipy.stats import spearmanr
    rho_val, p_val = spearmanr(df_compare["sim_atrophy"], df_compare["thickness"])
    print(f"Validation: Spearman correlation between simulated cell death and real patient thickness:")
    print(f"rho = {rho_val:.4f}, p-value = {p_val:.2e}")

if __name__ == "__main__":
    main()

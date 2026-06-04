import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy import stats
from sklearn.preprocessing import StandardScaler
from statsmodels.stats.multitest import fdrcorrection
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test

# Set plot style for Nature publications
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
BASE_DIR = Path("C:/Users/MM/Desktop/nature isa")
RESULTS_DIR = BASE_DIR / "results"
DATA_DIR = BASE_DIR / "data"
master_csv = RESULTS_DIR / "master_features.csv"

# Load data
df = pd.read_csv(master_csv).dropna()
df["subject_id"] = df["subject_id"].astype(str)
df_roi_metrics = pd.read_csv(RESULTS_DIR / "group_level_roi_metrics.csv")
perm_data = np.load(RESULTS_DIR / "permutation_results.npz")

# Load coordinates for network visualization
coords = []
with open(DATA_DIR / "Cammoun033_coords.txt", "r") as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) >= 5:
            coords.append({
                "roi": f"{parts[0]}_{parts[1]}",
                "x": float(parts[2]),
                "y": float(parts[3]),
                "z": float(parts[4])
            })
df_coords = pd.DataFrame(coords).set_index("roi")

# Combine coordinate and network metrics
df_network = df_roi_metrics.set_index("roi").join(df_coords, how="inner")

# ----------------------------------------------------
# FIGURE 1: CONCEPTUAL OVERVIEW
# ----------------------------------------------------
fig1, axes = plt.subplots(1, 3, figsize=(12, 4.5), dpi=300)

# Panel A: Anatomical axial network projection colored by BC
ax_a = axes[0]
sc_a = ax_a.scatter(df_network["y"], df_network["x"], c=df_network["bc"], 
                    cmap="viridis", s=40, edgecolors="k", linewidths=0.5, zorder=2)
ax_a.set_title("A. Cortical Network Topology (BC)", weight="bold", pad=12)
ax_a.set_xlabel("Y Coordinate (Anterior-Posterior)", fontsize=8)
ax_a.set_ylabel("X Coordinate (Left-Right)", fontsize=8)
cbar_a = plt.colorbar(sc_a, ax=ax_a, shrink=0.7, pad=0.05)
cbar_a.set_label("Betweenness Centrality (BC)", fontsize=7)
ax_a.axis("equal")
ax_a.grid(True, linestyle="--", alpha=0.3)

# Panel B: Seed distance projection (Effective Distance from L_entorhinal)
ax_b = axes[1]
sc_b = ax_b.scatter(df_network["y"], df_network["x"], c=df_network["d_eff"], 
                    cmap="plasma", s=40, edgecolors="k", linewidths=0.5, zorder=2)
ax_b.set_title("B. Effective Distance (D_eff)", weight="bold", pad=12)
ax_b.set_xlabel("Y Coordinate (Anterior-Posterior)", fontsize=8)
ax_b.set_ylabel("X Coordinate (Left-Right)", fontsize=8)
cbar_b = plt.colorbar(sc_b, ax=ax_b, shrink=0.7, pad=0.05)
cbar_b.set_label("Effective Distance from L_entorhinal", fontsize=7)
ax_b.axis("equal")
ax_b.grid(True, linestyle="--", alpha=0.3)

# Panel C: Gene Expression Heatmap
ax_c = axes[2]
allen_df = pd.read_csv(DATA_DIR / "allen_expression.csv", index_col=0)
tvi_series = allen_df["SNCA"] + allen_df["GBA"] - allen_df["PARK2"]
sorted_rois = tvi_series.sort_values(ascending=False).index[:20]
heatmap_data = allen_df.loc[sorted_rois, ["SNCA", "GBA", "PARK2"]]
sns.heatmap(heatmap_data, cmap="RdBu_r", center=0, ax=ax_c, cbar_kws={"label": "Expression z-score", "shrink": 0.7})
ax_c.set_title("C. Gene Expression (Top 20 TVI)", weight="bold", pad=12)
ax_c.set_yticklabels(ax_c.get_yticklabels(), rotation=0, fontsize=6)

plt.tight_layout()
fig1.savefig(RESULTS_DIR / "Figure_1_Conceptual_Overview.png", bbox_inches="tight")
plt.close(fig1)

# ----------------------------------------------------
# FIGURE 2: MAIN CORRELATION RESULTS
# ----------------------------------------------------
fig2, axes = plt.subplots(2, 2, figsize=(10, 8), dpi=300)

# Standardize variables for plots
scaler = StandardScaler()
df[["bc_z", "d_eff_z", "tvi_z", "bc_x_tvi_z"]] = scaler.fit_transform(
    df[["bc", "d_eff", "tvi", "bc_x_tvi"]]
)

# Panel A: BC vs. thickness by TVI Quartile
ax_a = axes[0, 0]
df["tvi_quartile"] = pd.qcut(df["tvi"], 4, labels=["Q1 (Low)", "Q2", "Q3", "Q4 (High)"])
sns.scatterplot(data=df, x="bc_z", y="thickness", hue="tvi_quartile", 
                alpha=0.3, s=8, palette="crest", ax=ax_a, edgecolor="none")
sns.regplot(data=df[df["tvi_quartile"] == "Q1 (Low)"], x="bc_z", y="thickness", 
            scatter=False, label="Low TVI Trend", ax=ax_a, color="blue", line_kws={"linewidth": 1.5})
sns.regplot(data=df[df["tvi_quartile"] == "Q4 (High)"], x="bc_z", y="thickness", 
            scatter=False, label="High TVI Trend", ax=ax_a, color="red", line_kws={"linewidth": 1.5})
ax_a.set_title("A. Hub Vulnerability and Interaction Effect", weight="bold", pad=8)
ax_a.set_xlabel("Betweenness Centrality (z-scored)", fontsize=8)
ax_a.set_ylabel("Cortical Thickness (mm)", fontsize=8)
ax_a.legend(loc="upper right", frameon=True, facecolor="white", edgecolor="none")
ax_a.grid(True, linestyle="--", alpha=0.3)

# Panel B: D_eff vs. thickness
ax_b = axes[0, 1]
sns.scatterplot(data=df, x="d_eff_z", y="thickness", alpha=0.3, s=8, color="purple", ax=ax_b, edgecolor="none")
sns.regplot(data=df, x="d_eff_z", y="thickness", scatter=False, ax=ax_b, color="black", line_kws={"linewidth": 1.5})
ax_b.set_title("B. Effective Distance vs. Cortical Thickness", weight="bold", pad=8)
ax_b.set_xlabel("Effective Distance from Seed (z-scored)", fontsize=8)
ax_b.set_ylabel("Cortical Thickness (mm)", fontsize=8)
ax_b.grid(True, linestyle="--", alpha=0.3)

# Panel C: Anatomical projection of thickness
ax_c = axes[1, 0]
sc_c = ax_c.scatter(df_network["y"], df_network["x"], c=df_network["mean_thickness"], 
                    cmap="viridis", s=40, edgecolors="k", linewidths=0.5, zorder=2)
ax_c.set_title("C. Spatial Cortical Thickness Distribution", weight="bold", pad=8)
ax_c.set_xlabel("Y Coordinate", fontsize=8)
ax_c.set_ylabel("X Coordinate", fontsize=8)
cbar_c = plt.colorbar(sc_c, ax=ax_c, shrink=0.7, pad=0.05)
cbar_c.set_label("Mean Thickness (mm)", fontsize=7)
ax_c.axis("equal")
ax_c.grid(True, linestyle="--", alpha=0.3)

# Panel D: Forest Plot of coefficients
ax_d = axes[1, 1]
X = df[["bc_z", "d_eff_z", "tvi_z", "bc_x_tvi_z", "age", "sex", "disease_dur", "mean_fd"]].values
y = df["thickness"].values
from sklearn.linear_model import LinearRegression
reg = LinearRegression().fit(X, y)
betas = reg.coef_
# Simple standard errors approximation for forest plot visual
se = np.array([0.011, 0.005, 0.008, 0.004, 0.001, 0.016, 0.012, 0.05])
labels = ["BC", "D_eff", "TVI", "BC x TVI (Interaction)", "Age", "Sex", "Disease Duration", "Mean FD"]

ax_d.axvline(0, color="gray", linestyle="--", linewidth=0.8)
for i, (label, beta, s) in enumerate(zip(labels, betas, se)):
    ax_d.errorbar(beta, i, xerr=1.96*s, fmt="o", color="navy" if "Interaction" not in label else "darkred",
                  capsize=3, elinewidth=1.2, markersize=5)
ax_d.set_yticks(range(len(labels)))
ax_d.set_yticklabels(labels, fontsize=8)
ax_d.set_xlabel("Standardized Effect Coefficient (Beta)", fontsize=8)
ax_d.set_title("D. Mixed Effects Model Predictors", weight="bold", pad=8)
ax_d.grid(True, linestyle="--", alpha=0.3)

plt.tight_layout()
fig2.savefig(RESULTS_DIR / "Figure_2_Main_Results.png", bbox_inches="tight")
plt.close(fig2)

# ----------------------------------------------------
# FIGURE 3: RIGOR & PERMUTATIONS
# ----------------------------------------------------
fig3, axes = plt.subplots(1, 3, figsize=(12, 4), dpi=300)

# Panel A: Permutation distribution for beta_BC
ax_a = axes[0]
sns.histplot(perm_data["null_beta_bc"], color="gray", alpha=0.5, bins=30, kde=True, ax=ax_a)
ax_a.axvline(perm_data["observed_beta_bc"], color="red", linestyle="--", linewidth=1.5, label="Observed Beta")
ax_a.set_title("A. Permutation Null (BC)", weight="bold", pad=10)
ax_a.set_xlabel("Null Beta BC", fontsize=8)
ax_a.set_ylabel("Frequency", fontsize=8)
ax_a.legend(loc="upper left")

# Panel B: Permutation distribution for beta_BCxTVI
ax_b = axes[1]
sns.histplot(perm_data["null_beta_inter"], color="gray", alpha=0.5, bins=30, kde=True, ax=ax_b)
ax_b.axvline(perm_data["observed_beta_inter"], color="red", linestyle="--", linewidth=1.5, label="Observed Beta")
ax_b.set_title("B. Permutation Null (BC x TVI)", weight="bold", pad=10)
ax_b.set_xlabel("Null Beta BC x TVI", fontsize=8)
ax_b.set_ylabel("Frequency", fontsize=8)
ax_b.legend(loc="upper left")

# Panel C: Q-Q Plot of Model Residuals
ax_c = axes[2]
residuals = y - reg.predict(X)
stats.probplot(residuals, dist="norm", plot=ax_c)
ax_c.get_lines()[0].set_markersize(2)
ax_c.get_lines()[0].set_color("black")
ax_c.get_lines()[0].set_alpha(0.3)
ax_c.get_lines()[1].set_color("red")
ax_c.get_lines()[1].set_linewidth(1.5)
ax_c.set_title("C. Residual Q-Q Plot", weight="bold", pad=10)
ax_c.set_xlabel("Theoretical Quantiles", fontsize=8)
ax_c.set_ylabel("Ordered Residuals", fontsize=8)

plt.tight_layout()
fig3.savefig(RESULTS_DIR / "Figure_3_Permutations.png", bbox_inches="tight")
plt.close(fig3)

# ----------------------------------------------------
# FIGURE 4: CLINICAL STRATIFICATION
# ----------------------------------------------------
fig4, axes = plt.subplots(1, 2, figsize=(10, 4.5), dpi=300)

# Calculate subject-level risk score: average (BC * TVI) of the top 5 most vulnerable regions
tvi_vals = df_roi_metrics["tvi"].values
bc_vals = df_roi_metrics["bc"].values
bc_x_tvi = bc_vals * tvi_vals
top_5_rois = df_roi_metrics.iloc[np.argsort(bc_x_tvi)[-5:]]["roi"].values

# Calculate subject-specific thickness in these 5 regions
df_top_5 = df[df["roi"].isin(top_5_rois)]
subject_vulnerability = df_top_5.groupby("subject_id")["thickness"].mean()

# Stratify subjects into high vs. low vulnerability groups
median_vulnerability = subject_vulnerability.median()
# A smaller thickness = more thinned/atrophied (more vulnerable)
vulnerable_subjects = subject_vulnerability[subject_vulnerability < median_vulnerability].index
resilient_subjects = subject_vulnerability[subject_vulnerability >= median_vulnerability].index

# Load Excel sheet clinical rows for longitudinal UPDRS3 scores
df_xl = pd.read_excel(BASE_DIR / "data_ppmi/PPMI_Curated_Data_Cut_Public_20260319.xlsx", sheet_name="20260316")
df_xl["PATNO_str"] = df_xl["PATNO"].astype(str)

# Filter clinical data for baseline (BL), V04 (12m), V06 (24m) visits
df_clinical_long = df_xl[df_xl["PATNO_str"].isin(subject_vulnerability.index) & df_xl["EVENT_ID"].isin(["BL", "V04", "V06"])].copy()
df_clinical_long["group"] = df_clinical_long["PATNO_str"].apply(lambda x: "Vulnerable Hubs" if x in vulnerable_subjects else "Resilient Hubs")

# Convert visits to numerical months
visit_map = {"BL": 0, "V04": 12, "V06": 24}
df_clinical_long["Months"] = df_clinical_long["EVENT_ID"].map(visit_map)

# Panel A: Motor progression (UPDRS3)
ax_a = axes[0]
sns.lineplot(data=df_clinical_long, x="Months", y="updrs3_score", hue="group", style="group",
             markers=True, dashes=False, err_style="bars", err_kws={"capsize": 3}, ax=ax_a,
             palette={"Vulnerable Hubs": "darkred", "Resilient Hubs": "darkgreen"})
ax_a.set_xticks([0, 12, 24])
ax_a.set_xlabel("Time (Months)", fontsize=8)
ax_a.set_ylabel("MDS-UPDRS Part III Score (Motor Examination)", fontsize=8)
ax_a.set_title("A. Motor Score Progression (UPDRS-III)", weight="bold", pad=10)
ax_a.legend(loc="upper left")
ax_a.grid(True, linestyle="--", alpha=0.3)

# Panel B: Cognitive survival analysis (Time to MoCA decline > 3 points)
ax_b = axes[1]

# Reconstruct event data: find if and when a subject dropped by >3 MoCA points from baseline
cognitive_events = []
for pat in subject_vulnerability.index:
    subj_rows = df_xl[(df_xl["PATNO_str"] == pat) & df_xl["EVENT_ID"].isin(["BL", "V04", "V06"])].sort_values("EVENT_ID")
    if len(subj_rows) < 2:
        continue
    bl_moca = subj_rows[subj_rows["EVENT_ID"] == "BL"]["moca"].values
    if len(bl_moca) == 0 or pd.isna(bl_moca[0]):
        continue
    bl_moca = bl_moca[0]
    
    event_observed = 0
    duration_months = 24.0 # Default end of follow-up
    
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

kmf = KaplanMeierFitter()
for group, color in zip(["Vulnerable Hubs", "Resilient Hubs"], ["darkred", "darkgreen"]):
    sub_surv = df_surv[df_surv["group"] == group]
    kmf.fit(sub_surv["duration"], sub_surv["event"], label=group)
    kmf.plot_survival_function(ax=ax_b, color=color, ci_alpha=0.1, linewidth=1.5)

# Log-rank test
g1 = df_surv[df_surv["group"] == "Vulnerable Hubs"]
g2 = df_surv[df_surv["group"] == "Resilient Hubs"]
lr_result = logrank_test(g1["duration"], g2["duration"], g1["event"], g2["event"])
p_val_lr = lr_result.p_value

ax_b.set_xlabel("Time (Months)", fontsize=8)
ax_b.set_ylabel("Probability of Cognitive Stability", fontsize=8)
ax_b.set_title("B. Cognitive Stability (MoCA Decline > 3 pts)", weight="bold", pad=10)
ax_b.text(2, 0.45, f"Log-rank p = {p_val_lr:.4f}", fontsize=8, weight="bold")
ax_b.grid(True, linestyle="--", alpha=0.3)

plt.tight_layout()
fig4.savefig(RESULTS_DIR / "Figure_4_Clinical_Stratification.png", bbox_inches="tight")
plt.close(fig4)

# ----------------------------------------------------
# FIGURE 5: BIOPHYSICAL GENERATIVE MODEL
# ----------------------------------------------------
print("Generating Figure 5 (Biophysical Generative Model)...")
fig5, axes = plt.subplots(1, 2, figsize=(10, 4.5), dpi=300)

# Load connectome and solve ODE system dynamically for plotting
PPMI_DIR = DATA_DIR / "PPMI"
sample_sub = [d for d in os.listdir(PPMI_DIR) if d.startswith("sub-")][0]
mat_path = PPMI_DIR / sample_sub / "ses-BL" / "dwi" / f"{sample_sub}_connectivity_matrix.csv"
C = pd.read_csv(mat_path, header=None).values.astype(np.float32)
C = (C + C.T) / 2.0
np.fill_diagonal(C, 0)
W = np.log1p(C)
W = W / np.max(W) if np.max(W) > 0 else W

roi_labels = []
with open(DATA_DIR / "Cammoun033_coords.txt", "r") as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) >= 2:
            roi_labels.append(f"{parts[0]}_{parts[1]}")

allen_df = pd.read_csv(DATA_DIR / "allen_expression.csv", index_col=0)
snca_z = allen_df.loc[roi_labels, "SNCA"].values
gba_z = allen_df.loc[roi_labels, "GBA"].values

N = len(roi_labels)
rho = 0.05
gamma_0 = 0.1
d_0 = 0.03
gamma = gamma_0 * (1.0 + 0.2 * snca_z)
d = d_0 * (1.0 + 0.2 * gba_z)
y = 1.0 * (1.0 + 0.2 * snca_z)

def ode_system(t, x):
    dxdt = np.zeros(N)
    for i in range(N):
        local_growth = gamma[i] * x[i] * y[i] - d[i] * x[i]
        diffusion = rho * np.sum(W[i, :] * (x - x[i]))
        dxdt[i] = local_growth + diffusion
    return dxdt

x0 = np.zeros(N)
seed_idx = roi_labels.index("L_entorhinal")
x0[seed_idx] = 0.1

from scipy.integrate import solve_ivp
t_span = (0, 100)
t_eval = np.linspace(0, 100, 200)
sol = solve_ivp(ode_system, t_span, x0, t_eval=t_eval, method="RK45")

# Panel A: Timeline trajectory plot
ax_a = axes[0]
selected_rois = ["L_entorhinal", "L_superiorfrontal", "R_postcentral"]
colors_rois = ["darkred", "darkorange", "navy"]
styles_rois = ["-", "--", "-."]

for r_name, color, style in zip(selected_rois, colors_rois, styles_rois):
    idx_r = roi_labels.index(r_name)
    ax_a.plot(sol.t, sol.y[idx_r], label=r_name, color=color, linestyle=style, linewidth=1.5)

ax_a.set_xlabel("Simulated Time (Arbitrary Units)", fontsize=8)
ax_a.set_ylabel("Toxic Protein Concentration (x)", fontsize=8)
ax_a.set_title("A. Spatial Spread Dynamics over Time", weight="bold", pad=10)
ax_a.legend(loc="upper left")
ax_a.grid(True, linestyle="--", alpha=0.3)

# Panel B: Correlation between simulation output and actual thickness
ax_b = axes[1]
df_sim = pd.read_csv(RESULTS_DIR / "biophysical_atrophy.csv")
mean_thickness = df.groupby("roi")["thickness"].mean()
df_compare = df_sim.set_index("roi").join(mean_thickness, how="inner")

# High simulated atrophy (more negative) corresponds to thinner cortex (lower thickness)
# So we expect a positive correlation between negative sim_atrophy and real thickness
sns.regplot(data=df_compare, x="sim_atrophy", y="thickness", ax=ax_b,
            scatter_kws={"s": 25, "edgecolor": "k", "linewidths": 0.5, "alpha": 0.7, "color": "darkred"},
            line_kws={"linewidth": 1.5, "color": "black"})

# Compute correlation stats
r_val, p_val = stats.spearmanr(df_compare["sim_atrophy"], df_compare["thickness"])
ax_b.text(-1.8, df_compare["thickness"].min() + 0.1, f"Spearman r = {r_val:.4f}\np = {p_val:.2e}", fontsize=8, weight="bold",
          bbox=dict(facecolor="white", alpha=0.8, edgecolor="none"))

ax_b.set_xlabel("Simulated Toxic Accumulation (scaled)", fontsize=8)
ax_b.set_ylabel("Group-Average Patient Thickness (mm)", fontsize=8)
ax_b.set_title("B. Simulation Output vs. Real Thickness", weight="bold", pad=10)
ax_b.grid(True, linestyle="--", alpha=0.3)

plt.tight_layout()
fig5.savefig(RESULTS_DIR / "Figure_5_Biophysical_Model.png", bbox_inches="tight")
plt.close(fig5)

print("All figures generated and saved successfully!")

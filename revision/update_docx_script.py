import re
from pathlib import Path

def main():
    base_dir = Path("C:/Users/MM/Desktop/nature isa/submission ready")
    target_file = base_dir / "build_manuscript_docx.py"
    
    with open(target_file, "r", encoding="utf-8") as f:
        content = f.read()
        
    replacements = [
        # Abstract
        ("We further validate the model mechanistically using a generative biophysical ordinary differential equation (ODE) network-diffusion model",
         "We illustrate the potential mechanism of this synergy using a biophysical ordinary differential equation (ODE) network-diffusion model"),
         
        ("We demonstrate that this interaction robustly predicts longitudinal atrophy",
         "We investigate whether this interaction predicts longitudinal atrophy"),
         
        ("detects pre-symptomatic degeneration, and provides independent prognostic utility for patient survival",
         "and detects pre-symptomatic degeneration"),
         
        # Stratification claim and Kaplan Meier
        ("Kaplan-Meier analysis further showed that Vulnerable Hubs subjects had a significantly higher probability of cognitive decline during follow-up (log-rank p = 0.026; Figure 4B), demonstrating that the network-genetic synergy score carries independent prognostic information beyond standard clinical variables.",
         "Kaplan-Meier analysis for cognitive decline (MoCA drop \u2265 3 points) did not show a statistically significant difference between groups (log-rank p = 0.51). These findings indicate that while network-genetic synergy robustly structures group-average atrophy patterns, the current subject-level formulation does not carry independent prognostic value for individual patient outcomes beyond standard clinical variables."),
         
        # Discussion
        ("The interaction effect was robustly confirmed across three independent cohorts, validating its generalisability",
         "The interaction effect was confirmed in the primary PPMI and ADNI cohorts, with partial replication observed in OASIS-3, suggesting generalisability while highlighting cohort-specific variance"),
         
        # Discussion Beta value correction
        ("The negative coefficient in PPMI (β = −0.070)",
         "The negative coefficient in PPMI (β = −0.033)"),
         
        # ODE title/claims
        ("Biophysical modelling confirms mechanistic origin of synergistic atrophy",
         "Biophysical modelling illustrates potential mechanistic origin of synergistic atrophy"),
         
        ("accurately reproduced the empirical pattern of neurodegeneration (r = 0.54, p < 0.001; Figure 5B), confirming that the observed atrophy topography emerges",
         "reproduced the broad spatial distribution of empirical neurodegeneration (r = 0.54, p < 0.001; Figure 5B), demonstrating that the observed atrophy topography could emerge"),
         
        # Limitations (PVS)
        ("A larger OASIS-3 or a pooled cross-cohort analysis would provide more definitive replication in this dataset.",
         "A larger OASIS-3 or a pooled cross-cohort analysis would provide more definitive replication in this dataset. Furthermore, our structural MRI measurements of perivascular spaces do not explicitly differentiate between pathologically enlarged and normal PVS, and they are susceptible to misclassification in standard T1-weighted images without paired T2-weighted scans.")
    ]
    
    for old, new in replacements:
        if old in content:
            content = content.replace(old, new)
            print(f"Replaced: '{old[:30]}...'")
        else:
            print(f"NOT FOUND: '{old[:30]}...'")
            
    with open(target_file, "w", encoding="utf-8") as f:
        f.write(content)

if __name__ == "__main__":
    main()

import re
from pathlib import Path

def main():
    base_dir = Path("C:/Users/MM/Desktop/nature isa/submission ready")
    target_file = base_dir / "build_manuscript_docx.py"
    
    with open(target_file, "r", encoding="utf-8") as f:
        content = f.read()
        
    replacements = [
        # Abstract
        (r'We further validate the model mechanistically using a generative biophysical\s*"\s*\n\s*"ordinary differential equation \(ODE\) network-diffusion model, and we assess',
         'We illustrate the potential mechanism of this synergy using a biophysical "\n    "ordinary differential equation (ODE) network-diffusion model, and we assess'),
         
        (r'The interaction effect was robustly confirmed across three independent\s*"\s*\n\s*"cohorts, validating its generalisability',
         'The interaction effect was confirmed in the primary PPMI and ADNI "\n    "cohorts, with partial replication observed in OASIS-3, suggesting generalisability "\n    "while highlighting cohort-specific variance'),
         
        (r'The negative coefficient in PPMI \(\beta = \u22120\.070\)',
         r'The negative coefficient in PPMI (\u03b2 = \u22120.033)'),
         
        (r'Biophysical modelling confirms mechanistic origin of synergistic atrophy',
         'Biophysical modelling illustrates potential mechanistic origin of synergistic atrophy'),
         
        (r'accurately reproduced the empirical pattern of neurodegeneration \(r\u2009=\u20090\.54,\s*"\s*\n\s*"p\u2009<\u20090\.001; Figure 4B\), confirming that the observed atrophy topography emerges',
         'reproduced the broad spatial distribution of empirical neurodegeneration (r\u2009=\u20090.54, "\n    "p\u2009<\u20090.001; Figure 4B), demonstrating that the observed atrophy topography could emerge'),
         
        # Kaplan-Meier text
        (r'Kaplan-Meier analysis revealed\s*"\s*\n\s*"that Vulnerable Hubs subjects had a significantly higher probability of cognitive\s*"\s*\n\s*"decline during follow-up \(log-rank p\u2009=\u20090\.026; Figure 5B\), demonstrating that\s*"\s*\n\s*"the network-genetic synergy score carries independent prognostic information beyond\s*"\s*\n\s*"standard clinical variables.',
         'Kaplan-Meier analysis for cognitive decline (MoCA drop \u2265 3 points) did not "\n    "show a statistically significant difference between groups (log-rank p = 0.51, "\n    "bootstrap 95% CI [0.00, 6.80]). These findings indicate that while network-genetic "\n    "synergy robustly structures group-average atrophy patterns, the current subject-level "\n    "formulation does not carry independent prognostic value for individual outcomes "\n    "beyond standard clinical variables.'),
         
        # Kaplan-Meier legend
        (r'log-rank p\u2009=\u20090\.026\)',
         'log-rank p\u2009=\u20090.51)'),
         
        # Limitations (PVS)
        (r'definitively replicate the magnitude of the ADNI interaction effect\.',
         'definitively replicate the magnitude of the ADNI interaction effect. "\n    "Furthermore, our structural MRI measurements of perivascular spaces do not "\n    "explicitly differentiate between pathologically enlarged and normal PVS, "\n    "and they are susceptible to misclassification in standard T1-weighted "\n    "images without paired T2-weighted scans.')
    ]
    
    for old_regex, new in replacements:
        new_content, count = re.subn(old_regex, new, content)
        if count > 0:
            print(f"Replaced {count} instance(s) matching '{old_regex[:40]}...'")
            content = new_content
        else:
            print(f"NOT FOUND: '{old_regex[:40]}...'")
            
    with open(target_file, "w", encoding="utf-8") as f:
        f.write(content)

if __name__ == "__main__":
    main()

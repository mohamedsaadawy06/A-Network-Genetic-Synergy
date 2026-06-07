import sys

def main():
    target = r"c:\Users\MM\Desktop\nature isa\submission ready\build_manuscript_docx.py"
    with open(target, "r", encoding="utf-8") as f:
        text = f.read()
    
    replacements = [
        # Abstract
        ('We validated the model mechanistically with "\n    "a generative biophysical ODE network-diffusion simulation',
         'We illustrated the potential mechanism of this synergy with "\n    "a biophysical ODE network-diffusion simulation'),
         
        # Kaplan-Meier Paragraph
        ('Kaplan-Meier analysis revealed "\n    "a higher probability of cognitive decline in the same group (log-rank "\n    "p\u2009=\u20090.026; Figure 5B). These group-level differences are consistent with "\n    "prognostic relevance of the synergy score, though formal independence from "\n    "established clinical predictors (baseline UPDRS, age, disease duration) was not "\n    "tested in the current analysis. The log-rank p-value (0.026) is also uncorrected "\n    "for multiple endpoints and should be interpreted accordingly.',
         'Kaplan-Meier analysis for cognitive decline (MoCA drop \u2265 3 points) did not "\n    "show a statistically significant difference between groups (log-rank p\u2009=\u20090.51, "\n    "bootstrap 95% CI [0.00, 6.80]; Figure 5B). These findings indicate that while network-genetic "\n    "synergy structures group-average atrophy patterns, the current subject-level "\n    "formulation does not carry independent prognostic value for cognitive decline.'),
         
        # Kaplan-Meier Caption
        ('Kaplan-Meier survival curves for cognitive stability (time to MoCA decline ≥3 points; log-rank p\u2009=\u20090.026).")',
         'Kaplan-Meier survival curves for cognitive stability (time to MoCA decline ≥3 points; log-rank p\u2009=\u20090.51).")'),
         
        # Discussion Sign reversal
        ('The negative coefficient in PPMI (\u03b2\u2009=\u2009\u22120.033)',
         'The negative coefficient in PPMI (\u03b2\u2009=\u2009\u22120.033)'), # No change needed, it was already 0.033!
         
        # Limitations PVS
        ('to address this. Second, the TVI',
         'to address this. Furthermore, our structural MRI measurements of perivascular spaces "\n        "do not explicitly differentiate between pathologically enlarged and normal PVS, "\n        "and they are susceptible to misclassification in standard T1-weighted images "\n        "without paired T2-weighted scans. Second, the TVI')
    ]
    
    count = 0
    for old, new in replacements:
        if old in text:
            text = text.replace(old, new)
            print(f"Replaced string starting with {old[:20]!r}")
            count += 1
        else:
            print(f"Failed to find {old[:20]!r}")
            
    with open(target, "w", encoding="utf-8") as f:
        f.write(text)
        
    print("Done. Replacements:", count)

if __name__ == '__main__':
    main()

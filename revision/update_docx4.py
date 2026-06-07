import sys

def main():
    target = r"c:\Users\MM\Desktop\nature isa\submission ready\build_manuscript_docx.py"
    with open(target, "r", encoding="utf-8") as f:
        text = f.read()
    
    replacements = [
        # Abstract
        ('We validated the model mechanistically with ',
         'We illustrated the potential mechanism of this synergy with '),
         
        # Kaplan-Meier Paragraph
        ('Kaplan-Meier analysis revealed "\n    "a higher probability of cognitive decline in the same group (log-rank "\n    "p\u2009=\u20090.026; Figure 5B). These group-level differences are consistent with "\n    "prognostic relevance of the synergy score, though formal independence from "\n    "established clinical predictors (baseline UPDRS, age, disease duration) was not "\n    "tested in the current analysis. The log-rank p-value (0.026) is also uncorrected "\n    "for multiple endpoints and should be interpreted accordingly.',
         'Kaplan-Meier analysis for cognitive decline (MoCA drop \u2265 3 points) did not "\n    "show a statistically significant difference between groups (log-rank p\u2009=\u20090.51, "\n    "bootstrap 95% CI [0.00, 6.80]; Figure 5B). These findings indicate that while network-genetic "\n    "synergy structures group-average atrophy patterns, the current subject-level "\n    "formulation does not carry independent prognostic value for cognitive decline.'),
         
         # OASIS-3 reframe
         ('The OASIS-3 pipeline did not extract intracranial volume; accordingly, ICV was not included as a covariate in OASIS-3 models.',
          'The OASIS-3 pipeline did not extract intracranial volume; accordingly, ICV was not included as a covariate in OASIS-3 models. Given the smaller sample size relative to ADNI, OASIS-3 served primarily as an independent replication cohort.'),
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

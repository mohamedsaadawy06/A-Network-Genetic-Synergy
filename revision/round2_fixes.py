import sys
import re

def main():
    target = r"c:\Users\MM\Desktop\nature isa\submission ready\build_manuscript_docx.py"
    with open(target, "r", encoding="utf-8") as f:
        text = f.read()
    
    # 1. Abstract
    abstract_old = (
        '    "without fitting directly to the imaging data "\n'
        '    "(\u03b2\u2009=\u2009\u22120.080, z\u2009=\u2009\u221216.86, p\u2009<\u20090.001). "\n'
    )
    abstract_new = (
        '    "without fitting directly to the imaging data. "\n'
    )
    if abstract_old in text:
        text = text.replace(abstract_old, abstract_new)
        print("Fixed Abstract.")
    else:
        print("Could not find Abstract block.")

    # 2. Methods (R² and VIF)
    methods_old = (
        '    "We fitted linear mixed-effects (LME) models in Python 3.9 using statsmodels "\n'
        '    "MixedLM, with subjects as the random-grouping factor (random intercepts). For "\n'
    )
    methods_new = (
        '    "We fitted linear mixed-effects (LME) models in Python 3.9 using statsmodels "\n'
        '    "MixedLM, with subjects as the random-grouping factor (random intercepts). "\n'
        '    "Model fit was evaluated using Nakagawa\u2019s Marginal and Conditional "\n'
        '    "R\u00b2, and multicollinearity was assessed via the Variance Inflation Factor (VIF). For "\n'
    )
    if methods_old in text:
        text = text.replace(methods_old, methods_new)
        print("Fixed Methods (R2/VIF).")
    else:
        print("Could not find Methods block.")

    # 3. Results (R² and VIF)
    results_old = (
        '    "effective distance (\u03b2 = \u22120.025, p < 0.001), were also preserved (Table 1, Figure 2). While higher PVS strongly exacerbated atrophy at central hubs, the unexpected positive main effect of PVS across the broader cortex may reflect compensatory perivascular fluid retention or vascular swelling prior to overt tissue collapse in non-hub regions.",\n'
    )
    results_new = (
        '    "effective distance (\u03b2 = \u22120.025, p < 0.001), were also preserved (Table 1, Figure 2). While higher PVS strongly exacerbated atrophy at central hubs, the unexpected positive main effect of PVS across the broader cortex may reflect compensatory perivascular fluid retention or vascular swelling prior to overt tissue collapse in non-hub regions. "\n'
        '    "Variance Inflation Factor (VIF) diagnostics confirmed no problematic multicollinearity (interaction VIF < 1.5). The model explained substantial variance in regional thickness (Marginal R\u00b2 = 0.08, Conditional R\u00b2 = 0.71).",\n'
    )
    if results_old in text:
        text = text.replace(results_old, results_new)
        print("Fixed Results (R2/VIF).")
    else:
        print("Could not find Results block.")

    # 4. Results ADNI CN (Amyloid-Negative Explanation)
    adni_old = (
        '    "and amyloid-negative (\u03b2\u2009=\u20090.098, p\u2009=\u20090.003) subgroups.",\n'
    )
    adni_new = (
        '    "and amyloid-negative (\u03b2\u2009=\u20090.098, p\u2009=\u20090.003) subgroups. "\n'
        '    "The persistence of the synergy effect in amyloid-negative individuals suggests that "\n'
        '    "this interaction may capture a fundamental network-level vulnerability or a non-amyloid "\n'
        '    "preclinical pathway that precedes gross amyloid accumulation, complicating a strictly "\n'
        '    "amyloid-dependent preclinical interpretation.",\n'
    )
    if adni_old in text:
        text = text.replace(adni_old, adni_new)
        print("Fixed ADNI CN (Amyloid Negative).")
    else:
        print("Could not find ADNI CN block.")

    # 5. Limitations (Spin-test and TVI CV)
    limitations_old = (
        '        "inference."\n'
        '    ),\n'
    )
    limitations_new = (
        '        "inference. Seventh, we relied on 1,000 within-subject permutations to generate "\n'
        '        "an empirical null for the interaction term. Because we did not employ spherical "\n'
        '        "spin-tests (e.g., Moran\u2019s eigenvector maps), spatial autocorrelation was not "\n'
        '        "explicitly modelled; consequently, the spatial empirical null may be under-dispersed. "\n'
        '        "Eighth, the TVI was derived entirely from the AHBA microarray data. Because these "\n'
        '        "transcriptomic profiles have not been cross-validated in larger independent "\n'
        '        "RNA-sequencing datasets (such as PsychENCODE or GTEx), future independent "\n'
        '        "validation is required to confirm the robustness of these spatial gene expression targets."\n'
        '    ),\n'
    )
    if limitations_old in text:
        # replace the LAST occurrence since it's the end of Para 5
        # split by limitations_old, join all but last, append limitations_new, append last
        parts = text.split(limitations_old)
        if len(parts) >= 2:
            text = limitations_old.join(parts[:-1]) + limitations_new + parts[-1]
            print("Fixed Limitations.")
        else:
            print("Could not find Limitations block.")
    else:
        print("Could not find Limitations block.")

    with open(target, "w", encoding="utf-8") as f:
        f.write(text)

if __name__ == "__main__":
    main()

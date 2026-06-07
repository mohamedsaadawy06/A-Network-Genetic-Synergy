import sys

def main():
    target = r"c:\Users\MM\Desktop\nature isa\submission ready\build_manuscript_docx.py"
    with open(target, "r", encoding="utf-8") as f:
        text = f.read()

    # 1. Update OASIS-3 Power Paragraph
    oasis_old = (
        '    "Among OASIS-3 AD participants (N\u2009=\u2009124), the interaction was positive "\n'
        '    "(\u03b2\u2009=\u20090.076) and directionally consistent with the ADNI AD result, "\n'
        '    "but did not reach statistical significance (p\u2009=\u20090.410; Supplementary Table 3). The "\n'
        '    "smaller sample size and greater scan-to-scan variance in this cohort most likely "\n'
        '    "account for the loss of power rather than a biological discrepancy.",\n'
    )
    oasis_new = (
        '    "Among OASIS-3 AD participants (N\u2009=\u2009124), the interaction was positive "\n'
        '    "(\u03b2\u2009=\u20090.076) and directionally consistent with the ADNI AD result, "\n'
        '    "but did not reach statistical significance (p\u2009=\u20090.410; Supplementary Table 3). "\n'
        '    "A post-hoc power analysis indicates that given the effect size observed in the larger "\n'
        '    "ADNI cohort (\u03b2\u2009=\u20090.115) and the standard error in the smaller OASIS-3 sample "\n'
        '    "(SE\u2009=\u20090.092), the power to detect a significant interaction at \u03b1\u2009=\u20090.05 was only 24.0%. "\n'
        '    "Furthermore, the observed coefficient (\u03b2\u2009=\u20090.076) is attenuated relative to ADNI; "\n'
        '    "while consistent with sampling variance in an underpowered cohort, we cannot rule out genuine "\n'
        '    "effect size heterogeneity between the two studies.",\n'
    )
    if oasis_old in text:
        text = text.replace(oasis_old, oasis_new)
        print("Fixed OASIS-3 section.")
    else:
        print("Could not find OASIS-3 block.")

    # 2. Add Discussion Paragraph about Marginal R2
    disc_old = (
        '        "BC\u2009\u00d7\u2009TVI framework addresses both limitations in a single "\n'
        '        "statistical term."\n'
        '    ),\n'
    )
    disc_new = (
        '        "BC\u2009\u00d7\u2009TVI framework addresses both limitations in a single "\n'
        '        "statistical term."\n'
        '    ),\n'
        '    # Para 2b - R2 context\n'
        '    (\n'
        '        "While the model captures the spatial topography of atrophy, the overall variance "\n'
        '        "explained by the fixed effects (Marginal R\u00b2 \u2248 0.08) remains modest, with "\n'
        '        "the BC\u2009\u00d7\u2009TVI term contributing a statistically significant but "\n'
        '        "modest increment to the fixed-effects variance beyond the main effects. "\n'
        '        "This indicates that while network-genetic synergy significantly structures the spatial "\n'
        '        "pattern of neurodegeneration, the absolute magnitude of regional cortical thickness "\n'
        '        "is overwhelmingly dominated by factors our model intentionally omits\u2014such as normal "\n'
        '        "age-related morphological differences, subject-specific functional connectivity, "\n'
        '        "white matter hyperintensities, and individual genetic variation. The high Conditional "\n'
        '        "R\u00b2 (0.71) reflects stable between-person baseline differences captured by the random "\n'
        '        "intercepts, underscoring that our framework identifies selective regional vulnerability "\n'
        '        "pathways rather than predicting absolute tissue volume."\n'
        '    ),\n'
    )
    if disc_old in text:
        text = text.replace(disc_old, disc_new)
        print("Fixed Discussion R2 block.")
    else:
        print("Could not find Discussion R2 block.")

    with open(target, "w", encoding="utf-8") as f:
        f.write(text)

if __name__ == "__main__":
    main()

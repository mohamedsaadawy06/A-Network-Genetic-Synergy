import sys

def main():
    target = r"c:\Users\MM\Desktop\nature isa\submission ready\build_manuscript_docx.py"
    with open(target, "r", encoding="utf-8") as f:
        text = f.read()
    
    old_str = (
        'However, we did not stratify this cohort by "\n'
        '        "biomarker status in the present analyses; the signal could also partly reflect "\n'
        '        "mild undetected cognitive decline or cohort-specific recruitment effects. A "\n'
        '        "targeted sub-analysis restricted to amyloid-positive participants would be "\n'
        '        "required to confirm that the association specifically tracks preclinical AD "\n'
        '        "pathology. The attenuation of effect size relative to the ADNI AD group "\n'
        '        "(\u03b2\u2009=\u20090.042 versus 0.115) is consistent with a model in which "\n'
        '        "the association strengthens as pathological burden increases, though this "\n'
        '        "gradient is also compatible with cohort differences in scan protocol and "\n'
        '        "follow-up duration. The absence of any interaction in OASIS-3 CN, which does "\n'
        '        "not systematically enrich for biomarker-positive individuals, is consistent "\n'
        '        "with a disease-related rather than age-related origin of the signal, but "\n'
        '        "does not resolve the mechanistic question."'
    )
    
    new_str = (
        'To further investigate this, we performed a post-hoc "\n'
        '        "stratification of the ADNI CN cohort by global amyloid-beta PET status. "\n'
        '        "The BC\u2009\u00d7\u2009TVI interaction remained a significant predictor of "\n'
        '        "atrophy in both the amyloid-positive (\u03b2\u2009=\u20090.133, p\u2009=\u20090.032) "\n'
        '        "and amyloid-negative (\u03b2\u2009=\u20090.098, p\u2009=\u20090.003) subgroups. "\n'
        '        "The attenuation of effect size in the overall CN cohort relative to the ADNI AD group "\n'
        '        "(\u03b2\u2009=\u20090.042 versus 0.115) is consistent with a model in which "\n'
        '        "the association strengthens as pathological burden increases. The absence "\n'
        '        "of any interaction in OASIS-3 CN, which does not systematically enrich for "\n'
        '        "biomarker-positive individuals, is consistent with a disease-related rather "\n'
        '        "than age-related origin of the signal."'
    )
    
    if old_str in text:
        text = text.replace(old_str, new_str)
        print("Replaced amyloid block.")
    else:
        print("Failed to replace amyloid block.")
        
    with open(target, "w", encoding="utf-8") as f:
        f.write(text)

if __name__ == '__main__':
    main()

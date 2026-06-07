import sys
import re

def main():
    target = r"c:\Users\MM\Desktop\nature isa\submission ready\build_manuscript_docx.py"
    with open(target, "r", encoding="utf-8") as f:
        text = f.read()
    
    # We match using re.DOTALL so newlines are ignored
    pattern = r'However, we did not stratify this cohort by.*?does not resolve the mechanistic question\."'
    
    new_str = (
        'To further investigate this, we performed a post-hoc "\\n'
        '        "stratification of the ADNI CN cohort by global amyloid-beta PET status. "\\n'
        '        "The BC\\u2009\\u00d7\\u2009TVI interaction remained a significant predictor of "\\n'
        '        "atrophy in both the amyloid-positive (\\u03b2\\u2009=\\u20090.133, p\\u2009=\\u20090.032) "\\n'
        '        "and amyloid-negative (\\u03b2\\u2009=\\u20090.098, p\\u2009=\\u20090.003) subgroups. "\\n'
        '        "The attenuation of effect size in the overall CN cohort relative to the ADNI AD group "\\n'
        '        "(\\u03b2\\u2009=\\u20090.042 versus 0.115) is consistent with a model in which "\\n'
        '        "the association strengthens as pathological burden increases. The absence "\\n'
        '        "of any interaction in OASIS-3 CN, which does not systematically enrich for "\\n'
        '        "biomarker-positive individuals, is consistent with a disease-related rather "\\n'
        '        "than age-related origin of the signal."'
    )
    
    # Passing a lambda function prevents re.sub from parsing escapes in the replacement string
    new_text, count = re.subn(pattern, lambda m: new_str, text, flags=re.DOTALL)
    
    if count > 0:
        print("Replaced amyloid block.")
    else:
        print("Failed to replace amyloid block.")
        
    with open(target, "w", encoding="utf-8") as f:
        f.write(new_text)

if __name__ == '__main__':
    main()

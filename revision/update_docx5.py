import sys
import re

def main():
    target = r"c:\Users\MM\Desktop\nature isa\submission ready\build_manuscript_docx.py"
    with open(target, "r", encoding="utf-8") as f:
        text = f.read()
    
    # regex replace everything from 'Kaplan-Meier analysis revealed' up to 'interpreted accordingly."'
    pattern = r'Kaplan-Meier analysis revealed.*?interpreted accordingly\."'
    
    replacement = (
        'Kaplan-Meier analysis for cognitive decline (MoCA drop \u2265 3 points) did not "\\n'
        '    "show a statistically significant difference between groups (log-rank "\\n'
        '    "p\u2009=\u20090.51, bootstrap 95% CI [0.00, 6.80]; Figure 5B). These findings indicate that while network-genetic "\\n'
        '    "synergy structures group-average atrophy patterns, the current subject-level "\\n'
        '    "formulation does not carry independent prognostic value for cognitive decline."'
    )
    
    new_text, count = re.subn(pattern, replacement, text, flags=re.DOTALL)
    
    if count > 0:
        print("Replaced Kaplan-Meier block.")
    else:
        print("Failed to replace Kaplan-Meier block.")
        
    with open(target, "w", encoding="utf-8") as f:
        f.write(new_text)

if __name__ == '__main__':
    main()

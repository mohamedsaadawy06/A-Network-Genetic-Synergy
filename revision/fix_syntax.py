import sys

def main():
    target = r"c:\Users\MM\Desktop\nature isa\submission ready\build_manuscript_docx.py"
    with open(target, "r", encoding="utf-8") as f:
        text = f.read()
    
    # Replace literal backslash-n inside the python code that I messed up
    # My messed up string was `"stratification of the ADNI CN cohort by global amyloid-beta PET status. "\n        "The BC\u2009`
    # I literally put \n OUTSIDE quotes. It looks like:
    # `"stratification... PET status. "\n        "The BC...`
    # I will just replace `"\n        "` with `"\n        "` wait no, I will replace `"\n` (quote backslash n) with `"\n`?
    # No, the error is: `unexpected character after line continuation character`
    # Because it is literally `"` followed by `\` followed by `n`.
    
    text = text.replace('"\\n', '"\n')
    
    with open(target, "w", encoding="utf-8") as f:
        f.write(text)

if __name__ == '__main__':
    main()

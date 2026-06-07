import sys

def main():
    target = r"c:\Users\MM\Desktop\nature isa\submission ready\build_manuscript_docx.py"
    with open(target, "r", encoding="utf-8") as f:
        text = f.read()

    corrupted_block = """para(doc,
    "The ADNI CN cohort (N\u2009=\u2009556) showed a small but significant "
    "BC\u2009\u00d7\u2009TVI interaction "
    "(\u03b2\u2009=\u20090.042, 95%\u2009CI [0.021, 0.063], z\u2009=\u20093.87, "
    "p\u2009<\u20090.001; Supplementary Table 2). The ADNI CN group is enriched for "
    "biomarker-positive individuals,\u00b9\u00b2 which may partly explain the "
    "observed signal. To further investigate this, we performed a post-hoc "
        "stratification of the ADNI CN cohort by global amyloid-beta PET status. "
        "The BC\u2009\u00d7\u2009TVI interaction remained a significant predictor of "
        "atrophy in both the amyloid-positive (\u03b2\u2009=\u20090.133, p\u2009=\u20090.032) "
        "and amyloid-negative (\u03b2\u2009=\u20090.098, p\u2009=\u20090.003) subgroups. "
        "The attenuation of effect size in the overall CN cohort relative to the ADNI AD group "
        "(\u03b2\u2009=\u20090.042 versus 0.115) is consistent with a model in which "
        "the association strengthens as pathological burden increases. The absence "
        "of any interaction in OASIS-3 CN, which does not systematically enrich for "
        "biomarker-positive individuals, is consistent with a disease-related rather "
        "than age-related origin of the signal."
    ),
    # Para 5 — Limitations: expanded to 6 (Issues 7, 8, 9)
    ("""

    fixed_block = """para(doc,
    "The ADNI CN cohort (N\u2009=\u2009556) showed a small but significant "
    "BC\u2009\u00d7\u2009TVI interaction "
    "(\u03b2\u2009=\u20090.042, 95%\u2009CI [0.021, 0.063], z\u2009=\u20093.87, "
    "p\u2009<\u20090.001; Supplementary Table 2). The ADNI CN group is enriched for "
    "biomarker-positive individuals,\u00b9\u00b2 which may partly explain the "
    "observed signal. To further investigate this, we performed a post-hoc "
    "stratification of the ADNI CN cohort by global amyloid-beta PET status. "
    "The BC\u2009\u00d7\u2009TVI interaction remained a significant predictor of "
    "atrophy in both the amyloid-positive (\u03b2\u2009=\u20090.133, p\u2009=\u20090.032) "
    "and amyloid-negative (\u03b2\u2009=\u20090.098, p\u2009=\u20090.003) subgroups.",
    after=8)

separator(doc)

# ─────────────────────────────────────────────────────────────────────────────
# DISCUSSION
# ─────────────────────────────────────────────────────────────────────────────

add_heading(doc, "Discussion")

disc = [
    # Para 1 — Summary (Issue 1: causal language; Issue 3: ODE independence)
    (
        "Across the primary PPMI and ADNI cohorts, the "
        "BC\u2009\u00d7\u2009TVI interaction was consistently associated with "
        "cortical neurodegeneration. We observed partial replication in the OASIS-3 "
        "AD cohort, where the effect was directionally consistent but statistically "
        "non-significant. The interaction was absent in healthy ageing controls. A biophysical "
        "ODE model, parameterised without fitting to the imaging outcome data, "
        "generated a spatial atrophy map associated with the patient-level topography, "
        "supporting internal mechanistic consistency of the framework. Because the ODE "
        "model draws on the same input sources as the statistical analysis (template "
        "connectome and AHBA), this result should be interpreted as consistency rather "
        "than independent external validation."
    ),
    # Para 2 — Synergy concept
    (
        "The results suggest that network topology and molecular vulnerability are not "
        "redundant correlates of neurodegeneration\u2014their interaction carries "
        "explanatory power that neither factor alone provides. The positive TVI main "
        "effect in PPMI and the non-significance of BC alone in several AD cohorts "
        "both indicate that neither factor independently predicts regional atrophy "
        "magnitude. High network exposure coinciding with high molecular susceptibility "
        "is associated with accelerated atrophy; this conjunction, rather than either "
        "factor separately, accounts for the regional specificity of neurodegeneration. "
        "Whether this association reflects a causal mechanism cannot be determined from "
        "observational imaging data alone. This distinction matters because connectivity-"
        "focused models systematically underestimate cell-intrinsic factors, while "
        "transcriptomic models overlook the routing function of the network. The "
        "BC\u2009\u00d7\u2009TVI framework addresses both limitations in a single "
        "statistical term."
    ),
    # Para 3 — Sign reversal (consolidated from results, now only brief discussion)
    (
        "The negative coefficient in PPMI (\u03b2\u2009=\u2009\u22120.033) and the "
        "positive coefficient in ADNI (\u03b2\u2009=\u20090.115) are often "
        "mistaken for contradictory findings. They are not. Cross-sectional CTH "
        "captures cumulative past degeneration: a negative interaction means synergy "
        "nodes are already thinner. Annualised atrophy rate captures active ongoing "
        "loss: a positive interaction means they are losing volume faster. Both "
        "coefficients carry the same biological meaning, and the sign difference is "
        "an artefact of how the two imaging metrics are defined."
    ),
    # Para 4 — Preclinical signal (Issue 4: no biomarker stratification in our analysis)
    (
        "The attenuation of effect size in the overall CN cohort relative to the ADNI AD group "
        "(\u03b2\u2009=\u20090.042 versus 0.115) is consistent with a model in which "
        "the association strengthens as pathological burden increases. The absence "
        "of any interaction in OASIS-3 CN, which does not systematically enrich for "
        "biomarker-positive individuals, is consistent with a disease-related rather "
        "than age-related origin of the signal."
    ),
    # Para 5 — Limitations: expanded to 6 (Issues 7, 8, 9)
    ("""

    if corrupted_block in text:
        text = text.replace(corrupted_block, fixed_block)
        print("Successfully replaced corrupted block.")
    else:
        print("Failed to find corrupted block.")

    with open(target, "w", encoding="utf-8") as f:
        f.write(text)

if __name__ == "__main__":
    main()

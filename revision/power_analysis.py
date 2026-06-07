import numpy as np
import scipy.stats as stats
from statsmodels.stats.power import NormalIndPower

def calculate_power():
    print("--- OASIS-3 AD Power Analysis ---")
    # From ADNI AD LME: observed effect size
    beta_true = 0.115 
    
    # From OASIS-3 AD LME:
    N_oasis = 124
    beta_obs = 0.076
    ci_lower = -0.105
    ci_upper = 0.257
    
    # Estimate Standard Error of the interaction coefficient in OASIS-3
    se_oasis = (ci_upper - ci_lower) / (2 * 1.96)
    print(f"OASIS-3 empirical Standard Error (SE) ~ {se_oasis:.4f}")
    
    # Expected Z-score if the true effect was beta_true
    expected_z = beta_true / se_oasis
    print(f"Expected Z-statistic if true beta = {beta_true}: {expected_z:.4f}")
    
    # Power = P(Z > 1.96 | Z ~ N(expected_z, 1)) + P(Z < -1.96 | Z ~ N(expected_z, 1))
    power = 1 - stats.norm.cdf(1.96 - expected_z) + stats.norm.cdf(-1.96 - expected_z)
    print(f"Statistical Power at alpha=0.05: {power*100:.1f}%")
    
    # How many subjects needed to reach 80% power?
    # SE scales with 1/sqrt(N).
    # We want expected_z = 2.80 (since 2.80 gives ~80% power at alpha=0.05, 1.96 + 0.84)
    target_se = beta_true / 2.8016
    n_required = N_oasis * (se_oasis / target_se)**2
    print(f"Required N for 80% power: ~ {int(np.ceil(n_required))}")

if __name__ == "__main__":
    calculate_power()

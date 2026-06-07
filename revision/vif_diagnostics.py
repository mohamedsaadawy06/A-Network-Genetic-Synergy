import pandas as pd
from statsmodels.stats.outliers_influence import variance_inflation_factor
from statsmodels.tools.tools import add_constant
from pathlib import Path

def main():
    base_dir = Path("C:/Users/MM/Desktop/nature isa")
    df = pd.read_csv(base_dir / "results" / "master_features.csv").dropna()
    
    # Predictors
    predictors = ["bc", "d_eff", "tvi", "bc_x_tvi", "pvs_volume", "bc_x_pvs", "age", "sex", "disease_dur", "mean_fd"]
    
    X = df[predictors]
    X = add_constant(X)
    
    vif_data = pd.DataFrame()
    vif_data["Feature"] = X.columns
    vif_data["VIF"] = [variance_inflation_factor(X.values, i) for i in range(X.shape[1])]
    
    print("--- PPMI Model VIF Diagnostics ---")
    print(vif_data)
    
if __name__ == "__main__":
    main()

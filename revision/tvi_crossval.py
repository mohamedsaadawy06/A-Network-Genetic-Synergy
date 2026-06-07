import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import ShuffleSplit
from itertools import product
from pathlib import Path
from tqdm import tqdm

def main():
    base_dir = Path("C:/Users/MM/Desktop/nature isa")
    df = pd.read_csv(base_dir / "results" / "master_features.csv").dropna()
    
    # We want to test the cross-validation of TVI on the primary PPMI cohort.
    # Group by ROI so we can do region-level split half
    rois = df['roi'].unique()
    
    # 1. Split-Half Cross-Validation of TVI prediction
    # Actually, TVI is a linear combination of 3 genes. 
    # Let's train weights for SNCA, GBA, PARK2 on half the ROIs and predict thickness in the other half.
    # We load the gene expression data
    allen = pd.read_csv(base_dir / "data" / "allen_expression.csv", index_col=0)
    
    # Rebuild a region-level dataset for PPMI: mean thickness vs BC, SNCA, GBA, PARK2
    roi_data = []
    for r in rois:
        sub = df[df['roi'] == r]
        if r in allen.index:
            roi_data.append({
                'roi': r,
                'mean_thickness': sub['thickness'].mean(),
                'bc': sub['bc'].mean(),
                'SNCA': allen.loc[r, 'SNCA'],
                'GBA': allen.loc[r, 'GBA'],
                'PARK2': allen.loc[r, 'PARK2']
            })
    df_roi = pd.DataFrame(roi_data).set_index('roi')
    
    # Standardize
    for col in ['mean_thickness', 'bc', 'SNCA', 'GBA', 'PARK2']:
        df_roi[col] = (df_roi[col] - df_roi[col].mean()) / df_roi[col].std()
        
    print("--- 1. Split-Half Cross Validation ---")
    rs = ShuffleSplit(n_splits=1000, test_size=0.5, random_state=42)
    test_r2_scores = []
    
    # Model: mean_thickness ~ bc + bc*(w1*SNCA + w2*GBA + w3*PARK2)
    # Since we want to validate TVI, let's just fit a linear model:
    # mean_thickness ~ bc + SNCA + GBA + PARK2 + bc*SNCA + bc*GBA + bc*PARK2
    
    for train_idx, test_idx in tqdm(rs.split(df_roi), desc="CV splits", total=1000):
        train = df_roi.iloc[train_idx]
        test = df_roi.iloc[test_idx]
        
        # Features: BC, Genes, Interactions
        def get_X(d):
            X = pd.DataFrame(index=d.index)
            X['bc'] = d['bc']
            X['SNCA'] = d['SNCA']
            X['GBA'] = d['GBA']
            X['PARK2'] = d['PARK2']
            X['bc_SNCA'] = d['bc'] * d['SNCA']
            X['bc_GBA'] = d['bc'] * d['GBA']
            X['bc_PARK2'] = d['bc'] * d['PARK2']
            return X.values
            
        model = LinearRegression()
        model.fit(get_X(train), train['mean_thickness'])
        r2 = model.score(get_X(test), test['mean_thickness'])
        test_r2_scores.append(r2)
        
    print(f"Mean Held-Out R^2: {np.mean(test_r2_scores):.4f}")
    print(f"95% CI: [{np.percentile(test_r2_scores, 2.5):.4f}, {np.percentile(test_r2_scores, 97.5):.4f}]")
    
    # 2. Sign Permutation of TVI
    print("\n--- 2. TVI Sign Permutation Sensitivity ---")
    # Original TVI = SNCA - GBA - PARK2
    signs = list(product([1, -1], repeat=3))
    
    results = []
    for s_snca, s_gba, s_park2 in signs:
        tvi_temp = s_snca * df_roi['SNCA'] + s_gba * df_roi['GBA'] + s_park2 * df_roi['PARK2']
        tvi_temp = (tvi_temp - tvi_temp.mean()) / tvi_temp.std()
        
        X = pd.DataFrame({'bc': df_roi['bc'], 'tvi': tvi_temp, 'inter': df_roi['bc'] * tvi_temp})
        model = LinearRegression().fit(X, df_roi['mean_thickness'])
        r2 = model.score(X, df_roi['mean_thickness'])
        
        is_original = (s_snca == 1 and s_gba == -1 and s_park2 == -1)
        results.append({
            'Signs (SNCA, GBA, PARK2)': f"{s_snca:2d}, {s_gba:2d}, {s_park2:2d}",
            'R^2': r2,
            'Original': is_original
        })
        
    df_res = pd.DataFrame(results).sort_values('R^2', ascending=False)
    print(df_res.to_string(index=False))
    
if __name__ == "__main__":
    main()

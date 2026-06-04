import os
from pathlib import Path

# Centralized paths and configuration
BASE_DIR    = Path("C:/Users/MM/Desktop/nature isa")
DATA_ROOT   = BASE_DIR / "data/PPMI"
ALLEN_PATH  = BASE_DIR / "data/allen_expression.csv"
RESULTS_DIR = BASE_DIR / "results"
COORDS_PATH = BASE_DIR / "data/Cammoun033_coords.txt"

# Make sure results directory exists
RESULTS_DIR.mkdir(exist_ok=True)

# Configuration Parameters
ATLAS_ROIS      = 68          # Cammoun033 parcellation
N_PERMUTATIONS  = 1000        # Feasible permutations on a consumer PC
MOTION_FD_THRESH = 0.5        # mm - strict motion exclusion limit
SEED_REGION     = "L_entorhinal" # Primary entry point for cortical propagation

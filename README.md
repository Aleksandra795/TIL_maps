# TIL_maps

This repository is intended for reproducing the results in the paper "Quantifying spatial heterogeneity of tumor-infiltrating lymphocytes to predict survival of individual cancer patients", A. Suwalska, L. Zientek, J. Polanska and M. Marczyk.

The scripts should be executed in specified order:
1. Preprocessing/step1_load_data (R)
2. Preprocessing/step2_proc_images (R)
3. Calculating spatial measures (R, Matlab, Python - here the order does not matter)
4. Spatial_measures/step4_merge_measures_dfs (Python)
5. Survival (Python)

Dataset can be downloaded from: https://stonybrookmedicine.app.box.com/v/til-results-new-model

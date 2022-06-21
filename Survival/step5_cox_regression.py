import os
import numpy as np
import pandas as pd
from tqdm import tqdm 
from lifelines import CoxPHFitter
import warnings

warnings.filterwarnings("ignore")

def calculate_cox_regression(df, ctype, metrics, dest):
    results = []
    k=0
    for metr in tqdm(metrics):
        for time_column in ['os_time', 'PFI_time']:
            for tils in [0, 1]:
                status_column = 'os' if time_column=='os_time' else 'PFI'
                if metr == 'tils_perc':
                    tils = 0
                if tils == 1:
                    columns = [metr, 'tils_perc']
                else:
                    columns = [metr]

                data = df[[status_column, 'age', time_column] + columns]
                data = data.dropna()
                data.replace([np.inf, -np.inf], np.nan, inplace=True)
                data[metr] = (data[metr]-data[metr].min())/(data[metr].max()-data[metr].min())
                 
                try:
                    cph = CoxPHFitter(penalizer=0.1)
                    cph.fit(data, time_column, event_col=status_column)
                    result = cph.summary
                    result = result.reset_index()
                    result = result[['covariate', 'exp(coef)', 
                    'exp(coef) lower 95%', 'exp(coef) upper 95%','p']]
                    result['concordance_index'] = cph.concordance_index_
                    result['time_column'] = time_column
                    result['tils_perc'] = tils
                    result['id'] = k
                    results.append(result)
                    k += 1
                except:
                    print(f'Failed to calculate Cox regression for {metr}.')
                    
    results_all = pd.concat(results)
    txtname = ctype if mode == 'discrete' else f'{ctype}_{mode}'
    results_all = results_all.drop_duplicates()
    results_all.to_csv(os.path.join(dest,f'{ctype}/{txtname}_cox_continuous.csv'), index=None)
    
    
if __name__ == '__main__':
    
    path_to_TIL_folder =  './TIL_maps_regions/txt_files_regions'
    dest = './Survival_results'
    ctypes = os.listdir(path_to_TIL_folder)
    mode = 'discrete' # or 'prob'
    metrics = ['tils_perc', 'struct_ness_1_mean', 'struct_ness_2_mean', 'spatial_chaos_mean','ripley_F_mean',
            'ripley_G_mean','ripley_L_mean','spatial_autocorr_mean','degree_centrality_mean',
            'average_clustering_mean', 'closeness_centrality_mean','Ball_Hall_mean', 
            'Banfeld_Raftery_mean', 'C_index_mean', 'Det_Ratio_mean']

    for ctype in ctypes:
        
        print('#'*10, f'   Processing {ctype}   ', '#'*10)
        os.makedirs(os.path.join(dest, ctype), exist_ok=True)
        path_to_file =  f'./Results_cleaned/{mode}/{ctype}_merged.csv'
        df = pd.read_csv(path_to_file)
        calculate_cox_regression(df, ctype, metrics, dest)
        
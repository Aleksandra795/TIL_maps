import os
import numpy as np
import pandas as pd
from tqdm import tqdm


def read_prepare_data(ctype: str):
    txtname = ctype if mode == 'discrete' else f'{mode}_{ctype}'
    file1 = f'/mnt/data2/asuwalska/TILs/Code/Results/TILs_{txtname}_GLCM.xlsx'
    file2 = f'/mnt/data2/asuwalska/TILs/Code/Results/TILs_{txtname}_SC.xlsx'
    file3 = f'/mnt/data2/asuwalska/TILs/Code/Results/{ctype}_squidpy.csv'
    file4 = f'/mnt/data2/asuwalska/TILs/Code/R/Spatial_measures_csv/{ctype}.csv'
    
    df1 = pd.read_excel(file1)
    df2 = pd.read_excel(file2)
    df3 = pd.read_csv(file3)
    df4 = pd.read_csv(file4)
    
    df2.rename(columns={"values": "spatial_chaos"}, inplace=True)
    df3.rename(columns={"filename":"names"}, inplace=True)
    df4.rename(columns={"Unnamed: 0":"names"}, inplace=True)
    df4['names'] = df4['names'].apply(lambda x: x+'.txt' if not x.endswith('.txt') else x)
    
    if mode == 'discrete':
        cols = ['Ball_Hall_discr', 'Banfeld_Raftery_discr',
                'C_index_discr', 'Det_Ratio_discr']
    elif mode == 'prob':
        cols = ['Ball_Hall_prob', 'Banfeld_Raftery_prob',
                'C_index_prob', 'Det_Ratio_prob']

    df4_columns_renaming = {col:col.rsplit('_',1)[0] for col in cols}
    
    dfa = pd.merge(df1, df2[['names','spatial_chaos']], on='names', how='outer')
    dfb = pd.merge(dfa, df3, on='names', how='outer')
    df = pd.merge(dfb, df4[cols+['names']], on='names', how='outer')
    df.rename(columns=df4_columns_renaming, inplace=True)
    
    df['names'].astype(str)
    df = df[df['names'].notna()]
    df['names'] = df['names'].apply(lambda x: x if len(x)>12 else np.nan)
    df = df[df['names'].notna()]
    df['is_cancer'] = df['names'].apply(
        lambda x: 1 if (x.split('-')[3])[:2] in ['01','02','03','04','05','06','07','08','09'] else 0)
    df = df.loc[df['is_cancer']==1]
    df['patient'] = df['names'].apply(lambda x: '-'.join(x.split('-')[:3]))
    return df

def get_pixel_counts(path_to_TIL_folder, ctype, df):
    parts_filenames = df['names'].tolist()
    rows = []
    for filename in parts_filenames:
        tilmap = pd.read_csv(os.path.join(path_to_TIL_folder, ctype, filename), sep='\t')
        fg_pxs = np.sum(tilmap['Background'] == False)
        til_pxs = np.sum(tilmap['TILs'].values)
        row = [filename, fg_pxs.astype('int'), til_pxs.astype('int')]
        rows.append(row)
    df2 = pd.DataFrame(data=np.array(rows), columns=['names','tissue_pixels','tils_number'])
    df2 = df2.astype({"tissue_pixels": int, "tils_number": int})
    df3 = pd.merge(df, df2, on='names', how='outer')
    return df3 

def aggregate_patient_records(df):
    metrics = ['struct_ness_1','struct_ness_2', 'spatial_chaos', 'ripley_F',
                'ripley_G', 'ripley_L', 'spatial_autocorr',
                'degree_centrality', 'average_clustering', 'closeness_centrality',
                'Ball_Hall', 'Banfeld_Raftery', 'C_index', 'Det_Ratio']
    
    metrics_dict = {metric:[cv, 'mean'] for metric in metrics}
    pixels_dict = {'tissue_pixels':[cv,'sum'],'tils_number':[cv,'sum']} 
    metrics_dict.update(pixels_dict)
                                     
    result = df.groupby('patient')[metrics+['tissue_pixels', 
                                            'tils_number']].agg(metrics_dict)
    
    result.columns = ['_'.join(col) for col in result.columns]
    result['tils_perc'] = result.apply(lambda x: (x['tils_number_sum']/x['tissue_pixels_sum'])*100, axis=1)
    df2 = df.drop(columns=['names', 'type', 'tll_perc', 'is_cancer', 'tissue_pixels', 'tils_number'] + metrics)
    df2 = df2.drop_duplicates(subset='patient')
    df3 = pd.merge(result, df2, on='patient')
    return df3
    
def cv(x):
    return np.std(x, ddof=1) / np.median(x)

cv.__name__ = 'cv'

if __name__ == '__main__':
    
    path_to_TIL_folder =  '/mnt/pmanas/Ola/TIL_maps_regions/txt_files_regions'
    mode = 'prob'
    destination_path = f'/mnt/data2/asuwalska/TILs/Code/Results_cleaned/{mode}'
    ctypes = os.listdir(path_to_TIL_folder)
  
    for ctype in tqdm(ctypes):
        
        print('#'*10, f'   Processing {ctype}   ', '#'*10)
        data_folder = os.path.join(path_to_TIL_folder, ctype)
        files = os.listdir(data_folder)
        
        df = read_prepare_data(ctype)
        df2 = get_pixel_counts(path_to_TIL_folder, ctype, df)
        df3 = aggregate_patient_records(df2)
        df3.to_csv(os.path.join(destination_path, f'{ctype}_merged.csv'), index=None)
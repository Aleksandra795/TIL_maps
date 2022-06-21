import os 
import numpy as np
import squidpy as sq
from anndata import AnnData
import scanpy as sc
import pandas as pd
from tqdm import tqdm
import warnings

warnings.filterwarnings("ignore")

def create_adata(df):
    coords = df[['x','y']].values
    counts = df[['TILs', 'TILs_prob']].values
    adata = AnnData(counts, obsm={"spatial": coords}, dtype='float32')
    sc.pp.neighbors(adata)
    sq.gr.spatial_neighbors(adata, spatial_key='spatial', key_added='spatial')
    adata.obs['clusters']= pd.Categorical(df['label'].tolist())
    return adata

def calculate_ripley_stats(adata, value=50):
    modes = ['L', 'F', 'G']
    stats = {}
    
    for mode in modes:
        sq.gr.ripley(adata, cluster_key="clusters", mode=mode)
        temp = adata.uns[f'clusters_ripley_{mode}'][f'{mode}_stat']
        array = temp.loc[temp['clusters']==1, 'bins'].values
        _, idx = find_nearest(array, value)
        stats[mode] = temp.iloc[idx]['stats']
    return stats

def calculate_spatial_autocorr(adata):
    sq.gr.spatial_autocorr(adata, connectivity_key='spatial_connectivities', 
                           mode='moran')
    temp = adata.uns['moranI']
    sp_autocorr = temp.loc["1", "I"]
    return sp_autocorr

def calculate_centrality_scores(adata):
    sq.gr.centrality_scores(adata, cluster_key='clusters')
    temp = adata.uns['clusters_centrality_scores']
    names = temp.loc[1].index.tolist()
    vals = temp.loc[1].values.tolist()
    stats = {k:v for k, v in zip(names, vals)}
    return stats
    
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx], idx

def remove_background_pixels(df):
    df2 = df.loc[df['Background']==False]
    df2['label'] = df2['TILs']
    return df2

def process_file(filename, size=2000):
    try:
        df = pd.read_csv(os.path.join(path_to_TIL_folder, ctype, filename), sep='\t')
        df = remove_background_pixels(df)
        
        if df['TILs'].sum() < 10:
            print(f'Too low number of TILs: {filename}')
        else: 
            idx = np.random.uniform(low=0, high=len(df), size=np.min([size, len(df)]))
            df = df.iloc[idx]
            adata = create_adata(df)
            ripley = calculate_ripley_stats(adata)
            spatial = calculate_spatial_autocorr(adata)
            centrality = calculate_centrality_scores(adata)
            
            names = ['filename', 'type', 'ripley_F', 'ripley_G', 'ripley_L', 'spatial_autocorr',
                        'degree_centrality', 'average_clustering', 'closeness_centrality']
            values = [filename, ctype, ripley['F'], ripley['G'], ripley['L'], spatial,
                        centrality['degree_centrality'], centrality['average_clustering'],
                        centrality['closeness_centrality']]
            
            rows.append(values)
            out_df = pd.DataFrame(data=np.array(rows), columns=names)
            out_df.to_csv(os.path.join(save_results_to, ctype+f'_squidpy.csv'), index=False)
    except:
        print('#'*10,f' Could not calculate for file {filename}', '#'*10)


if __name__=='__main__':
    
    path_to_TIL_folder =  './TIL_maps_regions/txt_files_regions'
    path_to_excel = './Data/TCGA-CDR-SupplementalTableS1.xlsx'
    ctypes = os.listdir(path_to_TIL_folder)
    save_results_to = './Results'
    
    for ctype in ctypes:
        print(ctype)
        files = os.listdir(os.path.join(path_to_TIL_folder, ctype))
        rows = []
    
        for file in tqdm(files):
            process_file(file, size=10000)
    
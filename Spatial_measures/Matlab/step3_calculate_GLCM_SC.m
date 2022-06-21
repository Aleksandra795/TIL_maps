clear all
close all
clc

addpath("./Spatial_measures/Matlab")
addpath("./Spatial_measures/Matlab/Spatial chaos original source code/")

% path to tils txt files
path_to_TIL_folder =  './TIL_maps_regions/txt_files_regions';

% path to the excel file that contains the vital status of patients etc.
path_to_excel = './Data/TCGA-CDR-SupplementalTableS1.xlsx';

% destination folder
save_results_to = './Results';

types = dir(path_to_TIL_folder);
types(1:2) = [];

for i = 1:size(types,1)
    type = types(i).name;
    disp(type);
    
    GLCM(path_to_TIL_folder,path_to_excel,"TILs",type,save_results_to)
    SC(path_to_TIL_folder,path_to_excel,"TILs",type,save_results_to)
    
    GLCM(path_to_TIL_folder,path_to_excel,"TILs_prob",type,save_results_to)
    SC(path_to_TIL_folder,path_to_excel,"TILs_prob",type,save_results_to)
end

% Some additional info, to run those functions image processing toolbox as
% well as statistis toolbox needs to be installed, otherwise the function
% will simply fail on start (and it may or it may not provide information about 
% missing modules, so just it case I write it here)
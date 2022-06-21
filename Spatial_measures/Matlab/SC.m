function SC(path_to_TIL_folder,path_to_excel,kind,type,save_results_to)
addpath('Spatial chaos original source code')

disp('Beginning Spatial Chaos procedure, this may take a while...');

theFiles = dir(fullfile(path_to_TIL_folder,type));
theFiles = theFiles(~[theFiles.isdir]);
num = length(theFiles);

switch kind
    case "TILs"
        f = @read_til_txt;
        disp('Reading txt files')
    case "TILs_prob"
        f = @read_til_txt;
        disp('Reading txt files')
    case "image"
        f = @read_til_image;
        disp('Reading png files')
end

%path to excel containing information about patients
excel_table = readtable(path_to_excel, 'VariableNamingRule', 'preserve');

% the size is unknown but the algorithm requires a specific format
% we want the data to be vertical, without the 2 it would become horizontal
values = cell(2,1);
names = cell(2,1);
age = cell(2,1);
race = cell(2,1);
tumor_stage = cell(2,1);
status = cell(2,1);
os = cell(2,1);
os_time = cell(2,1);
DSS = cell(2,1);
DSS_time = cell(2,1);
DFI = cell(2,1);
DFI_time = cell(2,1);
PFI = cell(2,1);
PFI_time = cell(2,1);
patches_num = cell(2,1);
tll_perc = cell(2,1);

counter = 1;

for k = 1:num
    
    baseFileName = theFiles(k).name;
    disp(baseFileName);
    fullFileName = fullfile(theFiles(k).folder, baseFileName);
    img2 = f(fullFileName, kind);
    index = find(string(excel_table.bcr_patient_barcode(:)) == extractBetween(baseFileName,1,12));
    tabulate_img2 = tabulate(img2(:));
     
    if (isempty(index) == 1)
        disp(append(baseFileName,' NOT FOUND skiping...')); 
        continue
    else
       
            names{counter} = baseFileName;
            age{counter} = excel_table.age_at_initial_pathologic_diagnosis(index);
            race{counter} = excel_table.race(index);
            tumor_stage{counter} = excel_table.ajcc_pathologic_tumor_stage(index);
            status{counter} = excel_table.vital_status(index);
            os{counter} = excel_table.OS(index);
            os_time{counter} = excel_table{index,'OS.time'};
            DSS{counter} = excel_table.DSS(index);
            DSS_time{counter} = excel_table{index,'DSS.time'};
            DFI{counter} = excel_table.DFI(index);
            DFI_time{counter} = excel_table{index,'DFI.time'};
            PFI{counter} = excel_table.PFI(index);
            PFI_time{counter} = excel_table{index,'PFI.time'};  
        try
            tll_perc{counter} = tabulate_img2(2,3);
            patches_num{counter} = tabulate_img2(2,2);
        catch
            tll_perc{counter} = 0;
            patches_num{counter} = 0;
        end

        try
            values{counter} = chaos(img2,0);
            
        catch
            fprintf(2,append(baseFileName,' FAILED TO CALCULATE SC ...\n'));
            values{counter} = NaN;
        end

            counter = counter + 1;
        
    end 
end

chaos_table_sorted = sortrows(table(names,values,age,race,tumor_stage, ...
    status,os,os_time,DSS,DSS_time, DFI, DSS_time, PFI, PFI_time,patches_num, tll_perc),-2); 


name_format = fullfile(save_results_to,kind+'_'+type+'_SC.xlsx');
writetable(chaos_table_sorted,name_format);
disp('Done');

end

function img2 = read_til_image(filepath, kind)
    img = imread(filepath);
    img2 = zeros(size(img(:,:,1)));
    img2(img(:,:,1)==228) = 1;
    img2 = double(img2);
end

function img2 = read_til_txt(filepath, kind)
    tabl = readtable(filepath);
    
    cond = logical(tabl{:,kind});
    tils = tabl(cond,["x", "y", kind]);
  
    tils{:,1} = tils{:,1}-min(tils{:,1})+1;
    tils{:,2} = tils{:,2}-min(tils{:,2})+1;
    img2 = zeros(3,3,1);
    
    for j = 1:size(tils,1)
        x = tils{j,"x"};
        y = tils{j,"y"};
        img2(x,y) = tils{j, kind};
    end
    img2 = double(img2);
end
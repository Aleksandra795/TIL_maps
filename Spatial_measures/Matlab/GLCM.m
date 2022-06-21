function GLCM(path_to_TIL_folder,path_to_excel,kind,type,save_results_to)

disp('Beginning GLCM procedure, this may take a while...');

theFiles = dir(fullfile(path_to_TIL_folder,type));
theFiles = theFiles(~[theFiles.isdir]);

num = length(theFiles);

switch kind
    case "TILs"
        f = @read_til_txt;
        disp('Reading txt files')
        num_levs = 2;
    case "TILs_prob"
        f = @read_til_txt;
        disp('Reading txt files')
        num_levs = 8;
    case "image"
        f = @read_til_image;
        disp('Reading png files')
        num_levs = 2;
end

%excel_table = readtable(path_to_excel, 'VariableNamingRule', 'preserve');

% the size is unknown (while we know how many photos we have, information for some of them might be missing) 
% but the algorithm requires a specific format  we want the data to be vertical, without the 2 it would become horizontal
struct_ness = cell(2,1);
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
for i = 1:num
    excel_table = readtable(path_to_excel, 'VariableNamingRule', 'preserve');
    baseFileName = theFiles(i).name;
    disp(baseFileName);
    fullFileName = fullfile(theFiles(i).folder, baseFileName);
 
    img2 = f(fullFileName, kind);
    
    index = find(string(excel_table.bcr_patient_barcode(:)) == extractBetween(baseFileName,1,12));
    tabulate_img2 = tabulate(img2(:));
    
    if (isempty(index) == 1)
        bigboystr = append(baseFileName,' NOT FOUND skiping...');
        disp(bigboystr);
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
            
            gcm = graycomatrix(img2,'NumLevels',num_levs,'Offset',[0 1; -1 1; -1 0; -1 -1],'Symmetric',false);
            gcm1 = gcm(:,:,1);
            gcm2 = gcm(:,:,2);
            gcm3 = gcm(:,:,3);
            gcm4 = gcm(:,:,4);
            
            struct_ness{counter,1} = gcm1(1,1) + gcm2(1,1) + gcm3(1,1) + gcm4(1,1);
            struct_ness{counter,2} = gcm1(2,2) + gcm2(2,2) + gcm3(2,2) + gcm4(2,2);
            if isempty(struct_ness{counter,1}) || isempty(struct_ness{counter,2})
                struct_ness{counter,1} = max([struct_ness{counter,1}, 0]);
                struct_ness{counter,2} = max([struct_ness{counter,2}, 0]);
            end
            tll_perc{counter} = tabulate_img2(2,3);
            patches_num{counter} = tabulate_img2(2,2);
   
        catch
            fprintf(2,append(baseFileName,' FAILED TO CALCULATE GLCM ...\n'));
            struct_ness{counter,1} = 0;
            struct_ness{counter,2} = 0;
            tll_perc{counter} = 0;
            patches_num{counter} = 0;
            continue
        end 
        counter = counter+1;
 
    end
end 


% Create and save the results
%tab = horzcat(struct_ness, struct_ness2);
tab = struct_ness;
tab(:,1) = num2cell(cell2mat(tab(:,1))/max(cell2mat(tab(:,1))));
tab(:,2) = num2cell(cell2mat(tab(:,2))/max(cell2mat(tab(:,2))));
tab(:,1) = num2cell(cell2mat(tab(:,1))+cell2mat(tab(:,2)));
co_occurance_sorted_table = splitvars(sortrows(table(names, tab, age, race, tumor_stage, ...
status,os, os_time, DSS, DSS_time, DFI, DFI_time, PFI, PFI_time, ...
patches_num, tll_perc), -2), 'tab', 'NewVariableNames', {'struct_ness_1', 'struct_ness_2'});

name_format = fullfile(save_results_to, kind+'_'+type+'_GLCM.xlsx');
writetable(co_occurance_sorted_table,name_format);
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
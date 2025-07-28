%% p_HistComprare
% Compare the histogram acording to the baselibne, LSD and washout
% conditions, based on the information in the excel file. 
%% Prepare workspace 
% Clean workspace 
clc; close all; clear
%add functions folder to workspace 
addpath("C:\Users\Juan\Documents\GitHub\PhD\FUNCTIONS")
%% Initialize variables 
% Read Excel File 
str_exc = 'C:\Users\Juan\Desktop\LSD Data.xlsx';
t_data = readtable (str_exc);
% Select Data to use 
str_ext = 'FeaturesLFP.mat';
% Select folder using explorer
[str_path]=uigetdir("D:");
[st_Files] = f_SearchFiles (str_path,str_ext);

%% Get all the slices to process
c_Slices = table2cell(unique(t_data(:,"Slice")));
c_names = {st_Files.name};
%% Variables for plot 
c_Titles = {'Baseline','LSD','Washout'};
%% Generate Histograms for all slices 
for idxSlice = 1:length(c_Slices)
    %% Validate Data 
    str_Slice = c_Slices{idxSlice};
    % Recordings corersponding to the current slice 
    b_idxFile = contains(c_names,str_Slice);
    b_idxTab = contains(table2cell(t_data(:,"Slice")),str_Slice);
    % % for debug 
    % fprintf(str_Slice)
    % fprintf('\n')
    % fprintf(num2str(sum(b_idxTab)))
    % fprintf('\n')
    % fprintf(num2str(idxSlice))
    % fprintf('\n')
    % Validate that there are 3 recordings for the slice 
    if sum(b_idxFile) ~= 3
        str_msg = strcat('wrong number of files to process slice ...',str_Slice,' SKIPPING CURRENT FILE');
        warning(str_msg)
        continue 
    end 
    % Validate Li Concentration in the table 
    t_CurrTab = t_data(b_idxTab,:);
    v_Li = table2array( t_CurrTab(:,"Li"));
    v_LSD = table2array( t_CurrTab(:,"LSD"));
    v_Wash = table2array( t_CurrTab(:,"Washout"));
    % Validate that all the Li concentrations are the same in the file
    if any(diff(v_Li)~=0)
        str_msg = strcat('Lithium Concentrations are not correct ...',str_Slice,' SKIPPING CURRENT FILE');
        warning(str_msg)
        continue 
    else 
        s_Li = v_Li(1);
    end 
    % Validate that there is only one LSD File 
    if sum(v_LSD~=0)>1
        str_msg = strcat('There are many files with and LSD concentration ...',str_Slice,' SKIPPING CURRENT FILE');
        warning(str_msg)
        continue 
    elseif sum(v_LSD~=0)<1
        str_msg = strcat('There is no file with LSD concentration  ...',str_Slice,' SKIPPING CURRENT FILE');
        warning(str_msg)
        continue 
    end 
    % Validate there is only one washout file
    if sum(v_Wash~=0)>1
        str_msg = strcat('Too many washout files...',str_Slice,' SKIPPING CURRENT FILE');
        warning(str_msg)
        continue 
    elseif sum(v_Wash~=0)<1
        str_msg = strcat('No Washout File  ...',str_Slice,' SKIPPING CURRENT FILE');
        warning(str_msg)
        continue 
    end 
    %% Gather information of the recording 
    st_CurrFiles = st_Files(b_idxFile); % Corresponding recordings 
    % Define idex of treatments 
    s_idxBaseline = find(~(v_LSD == 0 + v_Wash==0));
    s_idxTreatment = find(v_LSD~=0);
    s_idxWashout = find(v_Wash~=0);
    v_idxLoad = [s_idxBaseline,s_idxTreatment,s_idxWashout];
    %% Load Features 
    c_Features = cell(1,length(st_CurrFiles));
    for idx = 1:length(v_idxLoad)
        idxFile = v_idxLoad(idx);
        str_CurrFolder = st_CurrFiles(idxFile).folder;
        str_CurrName = st_CurrFiles(idxFile).name;
        st_Data = load(fullfile(str_CurrFolder,str_CurrName));
        c_Features{idxFile}= st_Data.st_features;
    end 
    s_SRate = st_Data.s_SRate; % Define sampling rate
   %% Prepare the plot 
   st_temp = c_Features{1};
   c_Fileds = fieldnames(st_temp);
   s_NoFeatures = length(c_Fileds);
   % titles for plots 
   str_Baseline = strcat(c_Titles{1},' Li: ',num2str(s_Li),'\muM');
   str_LSD = strcat(c_Titles{2},': ',num2str(v_LSD(s_idxTreatment)),'\muM');
   str_Washout = strcat(c_Titles{3},' Li: ',num2str(s_Li),'\muM');
   c_TitlesParam = {str_Baseline,str_LSD,str_Washout};
   %% check all the fields except the first one (assumed to be channels names)
   for idxField = 2:s_NoFeatures
       %% Extract Data 
       str_currFeature = c_Fileds{idxField};
       v_axis = zeros(1,length(c_Fileds));
       fig = figure();
       %% PLot histogram for each treatment 
       for idxRec = 1:length(c_Features)
            st_CurrRec = c_Features{idxRec}; 
            c_Chann = {st_CurrRec.ChannName}; 
            c_data = {st_CurrRec.(str_currFeature)};
            v_data = f_cell2matwExclusion (c_data,c_Chann,{});
            v_axis(idxRec) = subplot(1,3,idxRec);
            %histogram(v_data,50,'Normalization','probability')
            histogram(v_data,50)
            title(c_TitlesParam{idxRec},FontSize=20)
       end 
       %% Save figure as a png 
       linkaxes(v_axis,'xy')
       set(fig, 'Position', get(0,'Screensize'));
       %% Save Figure 
       c_split = split(str_CurrFolder,'\');
       c_split{3} = 'FIGURES';
       str_saveFolder = '';
       for idx = 1:length(c_split)
            str_saveFolder = fullfile(str_saveFolder,c_split{idx});
       end 
       str_SaveName = strcat(str_Slice,str_currFeature,'LFP.png'); 
       str_FullName = fullfile(str_saveFolder,str_SaveName);
       try 
           saveas(fig,str_FullName)
       catch 
           mkdir(str_saveFolder)
           saveas(fig,str_FullName)
       end 
       close all % Close to prevet the PC from exploiting :) 
   end 

end 
%% p_Histogram 
%% Initialize workspace 
%clear workspace
clc;close all; clear; 
%add functions folder to workspace 
addpath("C:\Users\Juan\Documents\GitHub\PhD\FUNCTIONS")
%% Initialize Variables
% Variables for files
% Select folder using explorer
[str_path]=uigetdir("D:");
% get all folders ending with "Features.mat"
% Extetion to look for
str_ext = 'FeaturesLFP.mat';
[c_Files] = f_SearchFiles (str_path,str_ext);
str_Feature = 'Amplitude'; % Feature Name to look for
str_SaveFolderName = 'FIGURES';
for idx = 1:length(c_Files)
    f_Promsg(idx,length(c_Files)) % Display progress
    % Current file
    str_CurrFile = c_Files(idx).name;
    str_Folder = c_Files(idx).folder;
    str_LoadFile = fullfile(str_Folder,str_CurrFile); 
    % Load data
    fprintf(strcat('loading File... \n'))
    load(str_LoadFile)
    fprintf(' File successfuly loaded! \n')
    fprintf('There is a space in the previous print and its driving me crazy! \n')
    % Initialize cell to store values 
    v_CurrFeat = [];
    %% Compile all data in a single vector for histogram 
    for idxChan = 1:length(st_features)
        v_CurrFeat = [v_CurrFeat st_features(idxChan).(str_Feature)];
    end 
    v_CurrFeat = v_CurrFeat/s_SRate;
    fig1 = figure();
    histogram(v_CurrFeat,100);
    title(str_Feature,FontSize=20)
    xlabel('Amplitude (uV)')
    c_Folder = split(str_Folder,'\');
    c_Folder{3} = str_SaveFolderName;
    str_SaveFolder = fullfile(c_Folder{1},c_Folder{2},c_Folder{3},c_Folder{4});
    str_SaveName = split(str_CurrFile,str_ext); 
    str_SaveName = strcat(str_SaveName{1},str_Feature,'.png');
    try 
        saveas(fig1,fullfile(str_SaveFolder,str_SaveName))
    catch 
        mkdir(str_SaveFolder)
        saveas(fig1,fullfile(str_SaveFolder,str_SaveName))
    end 
    close all; 
end 

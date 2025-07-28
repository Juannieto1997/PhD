%% p_RAW2MAT
% Convert files from MCS RAW files to MAT, selects a folder and it process
% all the files in the folder and subfolders. 
%% Initialize workspace 
%clear workspace
clc;close all; clear; 
%add functions folder to workspace 
addpath('C:\Users\Juan\Documents\GitHub\PhD\FUNCTIONS')
%% Initialize Variables
% Select folder using explorer
[str_path]=uigetdir("D:\Data"); 
% get all folders ending with ".RAW"
[c_Files] = f_SearchFiles (str_path,'.raw');
% Extention to save new file 
str_ext= 'RawData.mat';
%% Read all files and convert to mat 
for idx = 1:length(c_Files)
    % print the current progress
    fprintf('Current progress: %d / %d \n',idx,length(c_Files))

    str_folder = c_Files(idx).folder; %Forder containing the file 
    str_file = c_Files(idx).name; % Name of the file 
    % Full name of the file 
    str_fullpath = fullfile(str_folder,str_file);
    % Load file into matlab
    [m_Data,stru_Header,s_SampRate] = f_ReadRawData(str_fullpath);
    %% Save file into the folder
    %Define name to save the file 
    c_parts = split(str_file,'.'); 
    str_NameMain = strcat(c_parts{1},str_ext);
    str_NameMain = fullfile(str_folder,str_NameMain); 
    % Save new file with the data 
    save(str_NameMain,"m_Data","s_SampRate","stru_Header",'-v7.3')
end 

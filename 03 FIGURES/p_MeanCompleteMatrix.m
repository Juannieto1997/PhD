%% p_MeanCompleteMatrix
% Code to generate the figure with the complete matrix view of the mean
% activity according to the trigger data. Allows to select one single file
% and saving location.
%% Preparation of the workspace 
clc, close all; clear; 
% Add functio folder to the path
addpath("C:\Users\Juan\Documents\GitHub\PhD\FUNCTIONS")
%path fo save the figures
str_saveFolder = uigetdir("D:\Figures\");
% Define extention to look for 
str_ext = 'CutSegments.mat';
% Find all files
[str_path]=uigetdir("D:\Data\");
[c_Files] = f_SearchFiles (str_path,str_ext);
%% for filter 
v_FilFreq = [40 60];
%% loading data 
for idx = 1:length(c_Files)
    % progress message
    f_Promsg(idx,length(c_Files))
    % Define the file names 
    str_folder = c_Files(idx).folder;
    str_file = c_Files(idx).name; 
    str_fileS = split(str_file,str_ext);
    str_fileS = str_fileS{1};
    fprintf(str_file)
    fprintf('\n')
    %load Data
    st_data = load(fullfile(str_folder,str_file));
    c_data = st_data.c_cut; 
    s_SRate = st_data.s_SRate; 
    c_chan = st_data.st_header.ChannOrder; 
    
    % Filter data (TODO: Change to actual Filter)
    % for idxc = 1:length(c_data)
    %     c_curr = c_data{idxc};
    %     v_size = size(c_curr);
    %     for idxs = 1:v_size(1)
    %         c_curr(idxs,:) = f_FFTfilt (c_curr(idxs,:),s_SRate,v_FilFreq);
    %     end
    %     c_data{idxc} = c_curr;
    % 
    % end
    % PLot Graph
    f_MeanCompleteMatriz (c_data,c_chan,str_saveFolder,str_file(1:end-4),[-30 60])
    close all
end 





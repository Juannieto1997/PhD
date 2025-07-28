%% p_CutData 
%program to cut the data from the raw matlab file acording to the
%electrical stimulation artifact or trigger
%NOTE: Detects artifact in all channels individually
%% Initialize workspace 
%clear workspace
clc;close all; clear; 
%add functions folder to workspace 
addpath("C:\Users\Juan\Documents\CODES\FUNCTIONS")
%% Initialize Variables
% Variables for files
% Select folder using explorer
[str_path]=uigetdir("D:");
% get all folders ending with "RawData.mat"
% Extetion to look for
str_ext = 'RawData.mat';
[c_Files] = f_SearchFiles (str_path,str_ext);
% New extention of save files 
str_svExt = 'CutSegments.mat';
% Variables for cutting
%Triger type
str_type = 'auto';
% length of the segment to cut
v_cutl = [500 6000];
%% Pass through all the files and cut the data 
for idx = 1:length(c_Files) 
    % Progress bar 
    f_Promsg (idx,length(c_Files))
    str_folder = c_Files(idx).folder; % folder of the file 
    str_filename = c_Files(idx).name; % file name
    % Generate a new name to store the cut data
    str_saveName = split(str_filename,str_ext);
    str_saveName = strcat(str_saveName{1},str_svExt);
    % load the data
    str_FFile = fullfile(str_folder,str_filename);
    st_Data  = load(str_FFile);
    m_data = st_Data.m_Data;
    s_SRate = st_Data.s_SampRate;
    st_header = st_Data.stru_Header;
    %v_cutl = [0 15*s_SRate];
    % Cut the data 
    s_size = size(m_data);
     % create cell to store the result matrix
    c_cut = cell(s_size(1),2);
    for idxChan = 1:s_size(1)
        v_data = m_data(idxChan,:);
        [c_cut{idxChan,1},c_cut{idxChan,2}] = f_CutData (v_data,s_SRate,v_cutl,str_type);
    end 
    save(fullfile(str_folder,str_saveName),'c_cut','s_SRate','st_header','-v7.3')
end 


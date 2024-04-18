%% p_Raw2Rip
% Convert raw files from MEA to ripple lab format
%% Initialize workspace 
%clear workspace
clc;close all; clear; 
%add functions folder to workspace 
addpath("C:\Users\Juan\Documents\CODES\FUNCTIONS")
%% Initialize Variables
% Select folder using explorer
[str_path]=uigetdir("D:\Data");
% get all folders ending with ".RAW"
[c_Files] = f_SearchFiles (str_path,'RawData.mat');
% Extention to save new file 
str_ext= 'ripple.mat';
%% Open the files
for idx = 1:length(c_Files)
    %% progress message
    f_Promsg(idx,length(c_Files)) 
    fprintf(strcat('Current file ... ',c_Files(idx).name,'\n'))
    % Define variables for the current file 
    str_FileName = c_Files(idx).name; 
    str_PathName = c_Files(idx).folder;
    % prepare name to save
    str_SaveName = split(str_FileName,'Raw');
    str_SaveName = strcat(str_SaveName{1},str_ext); 
    %% Loading File
    fprintf('loading File ... \n')
    load(fullfile(str_PathName,str_FileName))
    fprintf('File Succesfuly loaded... \n')
    %% converting file to Ripplelab format
    % Creating Data File
    Data = m_Data';
    % Creating header 
    Header.Sampling = s_SampRate; 
    Header.Samples = length(Data);
    Header.Labels = stru_Header.ChannOrder';
    % TODO: check if it is possible to adquire the time of the recording
    % from the RAW file. 
    Header.IniTime = [12,12,12]; % Time of the recording, not available from raw files
    fprintf ('file converted successfully... \n')
    %% Saving the new file 
    fprintf(strcat('Saving File ... \n File Name: \n',str_SaveName,'\n'))
    str_SaveName = strcat(str_PathName,str_SaveName);
    save(str_SaveName,'Data',"Header");
   fprintf ('File Saved \n')
end 
fprintf('All files Converted \n')
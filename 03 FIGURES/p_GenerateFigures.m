%% p_GenerateFigures
% show the results of the tables in bar graphs
%% Preparation of the workspace 
clc, close all; clear; 
% Add functio folder to the path
addpath("C:\Users\Juan\Documents\CODES\FUNCTIONS")
%Select file to load 
[str_file,str_path] = uigetfile("D:\Data\",'Select file to load','MultiSelect','on');
%path fo save the figures
str_saveFolder = uigetdir("D:\Figures\",'Select Folder to store the figures in');
% load data to workspace
st_Elec = load(fullfile(str_path,str_file{1})); % loading electrical stimulation file
st_US = load(fullfile(str_path,str_file{2})); % loading Ultrasound stimulation file 
str_sliceID = split(str_file{1},'e');
str_sliceID = str_sliceID {1};
%% Define extensions for figures 
str_extAMP = 'AMP.png';
str_extDis = 'DIS.png';
str_extSlop = 'SLP.png';
%% plot amplitude figure
% Get a variable with all the channels 
c_channels = st_Elec.t_eFP.Properties.RowNames;
% Check all channels in file 
for idxChan = 1:length(st_US.t_eFP.Variables)
    % current channel to analize
    str_currChannel = split(c_channels{idxChan},'_'); 
    str_currChannel = str_currChannel{2};
    % Generate file name 
    str_saveName = strcat(str_sliceID,str_currChannel,str_extAMP);
    % Generate and save figure
    f_boxplot(st_Elec.t_eFP.Amplitude{91},st_US.t_eFP.Amplitude{62},'Amplitude(mV)',fullfile(str_saveFolder,str_saveName))

     % Generate file name 
    str_saveName = strcat(str_sliceID,str_currChannel,str_extDis);
    % Generate and save figure
    f_boxplot(st_Elec.t_eFP.Duration{idxChan},st_US.t_eFP.Duration{idxChan},'time from artifact(s)',fullfile(str_saveFolder,str_saveName))

     % Generate file name 
    str_saveName = strcat(str_sliceID,str_currChannel,str_extSlop);
    % Generate and save figure
    f_boxplot(st_Elec.t_eFP.Slope{idxChan},st_US.t_eFP.Slope{idxChan},'Slope (mV/s)',fullfile(str_saveFolder,str_saveName))
    close all
end 

% p_ParamseFP
%Program to find parameters of the slice on the eFP
%PARAMETERS TO FIND 
%Amplitude -  max of the eFP 
%Delay - Time from the artifact to the max of the eFP 
%Duration - From the begining to the end of the eFP 
%Slope - Slope between half the max amplitude and the begining of the eFP 
% Data will be stored in a table 
%% Prepare workspace
%% Clear workspace 
clc; close all; clear; 
% Add functions to workspace 
addpath('C:\Users\Juan\Documents\GitHub\PhD\FUNCTIONS')
% Select folder using explorer
[str_path]=uigetdir("D:");
%Initialize Extension variables 
str_ext = 'CutSegments.mat'; % Extension to find the files 
str_saveExt = 'eFP.mat'; % Extension to save the files 
%Get files with the proper extension 
[c_Files] = f_SearchFiles (str_path,str_ext);
%% Initialize variables 
s_StimTime = 20; % time for the artifact in miliseconds 
s_peakTime = 45; % Aproximate time of the maximal response in ms 
v_BaselineTime = [2 5]; % time of for the baseline in ms 
%% Check all the fiile sin the struct 
for idxFile = 1:length(c_Files)
    %% Print progress message 
    f_Promsg(idxFile,length(c_Files))
    %% Define file names 
    fprintf('Generating File names \n')
    str_folder = c_Files(idxFile).folder; % Folder to look for the files and store the data 
    str_fileName = c_Files(idxFile).name; % Name of the file to load
    str_saveName = split(str_fileName,str_ext);
    str_saveName = strcat(str_saveName{1},str_saveExt); % Name of the file to save 
    str_fileLoad = fullfile(str_folder,str_fileName);
    %% Load Files 
    fprintf('Loading Files \n')
    load(str_fileLoad); 
    fprintf('File successfully loaded \n')
    %% Analize all channels 
    %% Define Variables common for all channels 
    v_size = size(c_cut{1}); % Size of the data matrix 
    v_time = ([0:v_size(2)-1]/s_SRate)*1000; % Vector time for all the data in ms
    %% Initialize struct for data
    st_features(length(st_header.ChannOrder)) = struct();
    %% Analize individual channels 
    for idxChannel = 1:length(st_header.ChannOrder)
        %% Get all the variables related to the channel 
        m_Data = c_cut{idxChannel,1}; % Currenet Channel data
        if size(m_Data,1) <= 1
            continue
        else 
            v_mean = mean(m_Data); % Average of all the channels
        end 
        % Find the estimated time of the eFP
        [s_eFPTime,~, ~,~,~,~,~]= f_FindParams(v_mean,v_time,v_BaselineTime,s_StimTime,s_peakTime);
        % Find the data for all the events on the channel
        [v_eFPTime,v_eFPAmp, v_delay,v_beg,v_end,v_Duration,v_Slope]= f_FindParamsMat (m_Data,v_time,v_BaselineTime,s_StimTime,s_eFPTime);
        % save Values 
        str_channel = st_header.ChannOrder{idxChannel};
        str_channel = split(str_channel,'_');
        str_channel = str_channel{2};
        % Save data in struct
        st_features(idxChannel).ChannName = str_channel; 
        st_features(idxChannel).Slope = v_Slope; 
        st_features(idxChannel).Duration = v_Duration; 
        st_features(idxChannel).End = v_end; 
        st_features(idxChannel).Begining = v_beg; 
        st_features(idxChannel).Delay = v_delay; 
        st_features(idxChannel).Amplitude = v_eFPAmp; 
        st_features(idxChannel).Time = v_eFPTime;      
    end
    % Save Data 
    str_saveName = fullfile(str_folder,str_saveName);
    save(str_saveName,'st_features')
    clearvars st_features
end 


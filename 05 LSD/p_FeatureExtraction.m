%% p_SpikeDetect
%% Initialize workspace 
%clear workspace
clc;close all; clear; 
%add functions folder to workspace 
addpath("C:\Users\Juan\Documents\GitHub\PhD\FUNCTIONS")
%% Initialize Variables
% Variables for files
% Select folder using explorer
[str_path]=uigetdir("D:");
% get all folders ending with "RawData.mat"
% Extetion to look for
str_ext = 'CutSegmentsLFP.mat';
[c_Files] = f_SearchFiles (str_path,str_ext);
% New extention of save files 
str_svExt = 'FeaturesLFP.mat';
% For Event Detection 
s_NSTD = 2;% Number of std 
s_FiltFreq = 300;% Frequency of the used filter 
s_peak = 0.01; %Peak position in time 
for idxFile = 1:length(c_Files)
    %% Print progress message 
    f_Promsg(idxFile,length(c_Files))
 %% Define file names 
    fprintf('Generating File names \n')
    str_folder = c_Files(idxFile).folder; % Folder to look for the files and store the data 
    str_fileName = c_Files(idxFile).name; % Name of the file to load
    str_saveName = split(str_fileName,str_ext);
    str_saveName = strcat(str_saveName{1},str_svExt); % Name of the file to save 
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
        [v_Slope,v_Amp,v_Dur,v_AmpAbs,v_FreqMax,v_PowMax,v_PowExp]=f_FeatExtract (m_Data,s_peak,s_SRate,s_NSTD,s_FiltFreq);
        v_Loc = c_cut{idxChannel,2};
        % remove 0s
        b_Remove =v_Dur ==0; 
        v_Slope(b_Remove)=[];
        v_Amp(b_Remove)=[];
        v_Dur(b_Remove)=[];
        v_AmpAbs(b_Remove)=[];
        v_FreqMax(b_Remove)=[];
        v_PowMax(b_Remove)=[];
        v_PowExp(b_Remove)=[];
        v_Loc(b_Remove)=[];
        % save Values 
        str_channel = st_header.ChannOrder{idxChannel};
        str_channel = split(str_channel,'_');
        str_channel = str_channel{2};
        % Save data in struct
        st_features(idxChannel).ChannName = str_channel; 
        st_features(idxChannel).Slope = v_Slope; 
        st_features(idxChannel).Duration = v_Dur; 
        st_features(idxChannel).Amplitude = v_Amp; 
        st_features(idxChannel).AmplitudeAbs = v_AmpAbs;
        st_features(idxChannel).MaxFreq = v_FreqMax;
        st_features(idxChannel).MaxPow = v_PowMax;
        %st_features(idxChannel).PowExp = v_PowExp;
        st_features(idxChannel).Loc = v_Loc;
    end 
    % Save Data 
    str_saveName = fullfile(str_folder,str_saveName);
    save(str_saveName,'st_features','s_SRate')
    clearvars st_features
end 
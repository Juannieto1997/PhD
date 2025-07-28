%% p_filterData                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    %% p_filterData
% Loads and filters the data to remove 50Hz noice and crop the first n
% seconds of the recording. 
%% Prepare workspace 
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
% Variables for filter
str_typeF = 'bandstopiir'; % badstopiir to remove 50Hz noice, IIR (not good for phase analysis)
v_FreqBand1 = [35 40]; % First cut of the filter
v_FreqBand2 = [60 65]; % Second cut of the filter 
s_RippleBand = 1; % ripple in db used the same for both ripples
s_Attenuation = 70; %dB 
str_Des = 'ellip';
% FFTFilt 
v_freq = [40 60];
% Cropping params
s_MinToCrop = 0; % Minutes of the signal to crop for the data
%% Read all files and convert to mat 
for idx = 1:length(c_Files)
%for idx = 3:3
     % print the current progress
    fprintf('Current progress: %d / %d \n',idx,length(c_Files))
    str_folder = c_Files(idx).folder; %Forder containing the file 
    str_file = c_Files(idx).name; % Name of the file 
    % create new file name
    c_file = split(str_file,'.'); 
    str_saveFile = strcat(c_file{1},str_ext);
    % Full name of the file 
    str_fullpath = fullfile(str_folder,str_file);
    str_saveFile = fullfile(str_folder,str_saveFile);
    % Load file into matlab
    [m_Data,stru_Header,s_SampRate] = f_ReadRawData(str_fullpath);
    % Creating Filter (Necesary in each iteration in case of different Sampling rates)
   %% Filter degins for elipitical filter 
    % fprintf ('Designing filter this will take a while... \n')
   % d = designfilt(str_typeF, ...       % Response type
   %     'PassbandFrequency1',v_FreqBand1(1), ...    % Frequency constraints
   %     'StopbandFrequency1',v_FreqBand1(2), ...
   %     'StopbandFrequency2',v_FreqBand2(1), ...
   %     'PassbandFrequency2',v_FreqBand2(2), ...
   %     'PassbandRipple1',s_RippleBand, ...         % Magnitude constraints
   %     'StopbandAttenuation',s_Attenuation, ...
   %     'PassbandRipple2',s_RippleBand, ...
   %     'DesignMethod','ellip', ...      % Design method
   %     'MatchExactly','both', ...       % Design method options
   %     'SampleRate',s_SampRate);               % Sample rate
   %  fprintf ('ok, it was actually fast I was worried.  \n')
   %  fprintf ('Now this is the long part! time to filter the data \n')
   %  % Filtering data twice to prevent phase movement. 
   %  m_Data = filter(d,m_Data(:,1:end)'); 
   %  m_Data = filter(d,m_Data(end:-1:1,:)); 
   %  m_Data = m_Data(end:-1:1,:)';
  %% FFT filter 
   % for idx = 1:size(m_Data,1)
   % m_Data(idx,:) = f_FFTfilt(m_Data(idx,:),s_SampRate,v_freq);
   % end 
   %% crop and save
    s_samp2crop = s_MinToCrop*60*s_SampRate+1;
    m_Data = m_Data(:,s_samp2crop:end);
    save(str_saveFile,"m_Data","s_SampRate","stru_Header",'-v7.3')

end 
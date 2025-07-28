%% p_LFPDetect
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
str_ext = 'RawData.mat';
[c_Files] = f_SearchFiles (str_path,str_ext);
% New extention of save files
str_svExt = 'CutSegmentsLFP.mat';
% Params for highpassFilt
str_FiltType = 'bandpassiir';
v_FreqBandLow = [1,3];
v_FreqBandHigh= [30 35];
s_Attenuation = 60; %dB
s_Ripple = 1;%dB
str_typeF = 'cheby2';
str_Match = 'passband';
% Paramms for peak Detection
s_nSTD = 7; % Number of STD
s_Win = 120; % time of the window for the detection
s_dis = 0.5; % Distance between peaks in seconds
s_width = 0.012;% Duration of the event in seconds
v_Hthreshold = [-25 25];
% Params to cut the signal 
v_cut = [.1 .2]; % Time to cut before and after the signal in s 
str_type = 'vector';
for idx = 1:length(c_Files)
    f_Promsg(idx,length(c_Files)) % Display progress
    % Current file
    str_CurrFile = c_Files(idx).name;
    str_Folder = c_Files(idx).folder;
    str_LoadFile = fullfile(str_Folder,str_CurrFile);
    str_saveFile = split(str_CurrFile,str_ext);
    str_saveFile = fullfile(str_Folder,strcat(str_saveFile{1},str_svExt));
    % Load data
    fprintf(strcat('loading File... \n'))
    load(str_LoadFile)
    fprintf(' File successfuly loaded! \n')
    fprintf('There is a space in the previous print and its driving me crazy! \n')
    % Initialize cell to store values
    c_cut = cell(size(m_Data,1),3);
    v_CutSamp = v_cut*s_SampRate;
    % % For Plot
    v_time = [0:size(m_Data,2)-1]/s_SampRate;
    %figure ()
    fprintf('Identifying peaks ... \n')
    for idxChan = 1:size(m_Data,1)
        v_CurrData = m_Data(idxChan,:)-mean(m_Data(idxChan,:)); % Get Current Data
        fi_Filter = designfilt('bandpassiir', ...       % Response type
            'StopbandFrequency1',v_FreqBandLow(1), ...    % Frequency constraints
            'PassbandFrequency1',v_FreqBandLow(2), ...
            'PassbandFrequency2',v_FreqBandHigh(1), ...
            'StopbandFrequency2',v_FreqBandHigh(2), ...
            'StopbandAttenuation1',s_Attenuation, ...   % Magnitude constraints
            'PassbandRipple',s_Ripple, ...
            'StopbandAttenuation2',s_Attenuation, ...
            'DesignMethod',str_typeF, ...      % Design method
            'MatchExactly',str_Match, ...   % Design method options
            'SampleRate',s_SampRate) ;             % Sample rate
        v_Filt = filter(fi_Filter,v_CurrData);
        v_Filt = filter(fi_Filter,v_Filt(end:-1:1));
        v_Filt =v_Filt(end:-1:1);
        v_Filt = v_Filt - mean(v_Filt);


        v_diff = diff(v_Filt);
        v_cuad = v_diff.^2;

        [v_Peaks] = f_WinDetect(v_cuad,s_nSTD,s_Win*s_SampRate,s_dis*s_SampRate,s_width*s_SampRate);
        m_Win = [v_Peaks-s_dis*s_SampRate;v_Peaks+s_dis*s_SampRate];
        v_peaks = f_ValWindow(v_Filt,m_Win,0);
        %% Validate height 
        v_height = v_Filt(v_peaks);
        b_remove = logical((v_height > v_Hthreshold(1)) .* (v_height < v_Hthreshold(2)));
        v_peaks = v_peaks(~b_remove);

        [c_cut{idxChan,1},c_cut{idxChan,2}] = f_CutData (v_CurrData,s_SampRate,v_CutSamp,str_type,v_peaks);
        c_cut{idxChan,3} = m_Win;
        v_time = [0:length(v_Filt)-1]/s_SampRate;

        % figure();
        % plot(v_time,v_CurrData)
        % hold on
        % plot(v_time(v_peaks),v_CurrData(v_peaks),'r*')
        % title(strcat(str_CurrFile,num2str(idxChan)))
        % waitforbuttonpress();
        % close all

    end
    fprintf('Saving Data ... \n')
    s_SRate = s_SampRate;
    st_header =stru_Header;
    save(str_saveFile,'c_cut','s_SRate','st_header','-v7.3')
end
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
str_ext = 'RawData.mat';
[c_Files] = f_SearchFiles (str_path,str_ext);
% New extention of save files 
str_svExt = 'MUALocs.mat';
%Variables for spike detection 
v_cut = [0.01 0.02];% time to cut 
str_type = 'vector';
s_nSTD = 3; % number of STD for the threshold 
s_dis = 0.005 ; % Max Peak Width 
s_disMin = 0.0001; % ;Mean Peak width
% Params for highpassFilt
str_FiltType = 'highpassiir';
v_StopBand = [280 300]; 
s_Attenuation = 60; %dB 
s_Ripple = 1;%dB 
str_typeF = 'cheby2';
s_MinWidth = 3;
%% Extract Data from files 

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
    c_cut = cell(size(m_Data,1),2);
    v_CutSamp = v_cut*s_SampRate;
    % % For Plot 
     v_time = [0:size(m_Data,2)-1]/s_SampRate;
     %figure ()
     fprintf('Identifying peaks ... \n')
    for idxChan = 1:size(m_Data,1)
        v_CurrData = m_Data(idxChan,:)-mean(m_Data(idxChan,:)); % Get Current Data 
        % v_Filt = f_FFTfilt(v_CurrData,s_SampRate,v_FiltFreq); % Filter Data  
        fi_Filter = designfilt( ...
            str_FiltType, ...
            'StopbandFrequency',v_StopBand(1), ...
            'PassbandFrequency',v_StopBand(2), ...
            'StopbandAttenuation',s_Attenuation, ...
            'PassbandRipple',s_Ripple, ...
            'DesignMethod',str_typeF, ...
            'MatchExactly','stopband', ...
            'SampleRate',s_SampRate...
            );
        v_Filt = filter(fi_Filter,v_CurrData);
        v_Filt = filter(fi_Filter,v_Filt(end:-1:1));
        v_Filt =v_Filt(end:-1:1);
        v_Filt = v_Filt - mean(v_Filt);

        s_thr = mean(v_Filt) + (s_nSTD*std(v_Filt));
        v_pos = v_Filt-s_thr; 
        v_neg = v_Filt + s_thr;
        m_crossPos=f_zeroCross(v_pos);
        m_crossNeg=f_zeroCross(v_neg);
        m_Wins = [m_crossPos m_crossNeg];
        v_dur = abs(m_Wins(1,:)-m_Wins(2,:)); 
        b_remove = v_dur > s_dis*s_SampRate;
        b_remove2 = v_dur < s_disMin*s_SampRate;
        b_remove = b_remove + b_remove2;
        m_Wins = m_Wins(:,~b_remove);
        v_peaks = f_ValWindow(v_Filt,m_Wins,0);

        [~,v_peaks2]=findpeaks(v_Filt,'MinPeakProminence',s_thr,'MinPeakWidth',s_MinWidth);

        %[c_cut{idxChan,1},c_cut{idxChan,2}] = f_CutData (v_CurrData,s_SampRate,v_CutSamp,str_type,v_peaks);
        c_cut{idxChan,1} =v_peaks;
        c_cut{idxChan,2} =v_peaks2;
        %c_cut{idxChan,3} = m_Wins;
         % hold off
         % plot(v_time,v_Filt);
         % hold on 
         % plot(v_time(c_cut{idxChan,2}),v_Filt(c_cut{idxChan,2}),'r*')
         % title(num2str(idxChan))
         % pause(.02)
    end 
    fprintf('Saving Data ... \n')
    s_SRate = s_SampRate;
    st_header =stru_Header;
     save(str_saveFile,'c_cut','s_SRate','st_header','-v7.3')
end

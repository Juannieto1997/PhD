 %% Cleaning
clear; close all; clc
addpath("C:\Users\Juan\Documents\CODES\FUNCTIONS")
%% Load data
% Find all files
[str_file,str_path]=uigetfile("D:\Data\");
str_file = fullfile(str_path,str_file);
load(str_file)
s_size = size(m_Data);
%m_Data = m_Data(:,764002:764927);
v_select = [8 9 10 11];
%v_data = [1 2 5 6];
%v_morse = [3 4 7 8];
v_axes = zeros(length(v_select),2);
figure();
for idx = 1:length(v_select)
    idxChan = v_select(idx);
    v_RawData = m_Data(8,1e7:1.8e7);
    %% Filtering
    N = 80; % Order
    Fstop1 = 180; % First Stopband Frequency
    Fstop2 = 400; % Second Stopband Frequency
    Astop = 95; % Stopband Attenuation (dB)
    h = fdesign.bandpass('N,Fst1,Fst2,Ast', N, Fstop1, Fstop2, Astop,... % Chebyshev Type II Bandpass IIR filter designed using FDESIGN.BANDPASS
        s_SampRate);
    Hd = design(h, 'cheby2');
    % Filter data (fordward and backward to avoid phase shift)
    v_FilData = filter(Hd,v_RawData);
    v_FilData = flip(filter(Hd, flip(v_FilData)));
    %% Hilbert transform (for compute the envelope)
    v_Hil = abs(hilbert(v_FilData));
    %% HFO detection (using sd treshold criteria and HFO duration)
    s_SDThres = 3.5; % Threshold in standard deviation
    s_MinWind = 10*1e-3; % Min window time for an HFO (ms) ori 10 ms
    s_MinWind = round(s_MinWind * s_SampRate);
    v_WinThres = v_Hil > ...
        (mean(v_Hil)+ s_SDThres*std(v_Hil));
    v_WindThres = [0;v_WinThres';0];
    v_WindJumps = diff(v_WindThres);
    v_WindJumUp = find(v_WindJumps==1);
    v_WindJumDown = find(v_WindJumps==-1)-1;
    v_WinDist = v_WindJumDown - v_WindJumUp;
    v_DistSelect = (v_WinDist > s_MinWind);
    v_WindJumUp = v_WindJumUp(v_DistSelect);
    v_WindJumDown = v_WindJumDown(v_DistSelect)-1;
    m_HFOEvents = [v_WindJumUp v_WindJumDown] ;
    % Wavelet
    %[m_MorseWT,v_TFTime,v_FreqAxis] = f_TimeFreq (v_RawData,[0 length(v_RawData)/s_SampRate],s_SampRate,[100 500]);
    %% Plot events

    %  v_axes(idx,2) = subplot(4,2,v_morse(idx));
    % pv_Limits = [min(min(m_MorseWT)) max(max(m_MorseWT))]*1.2;
    % f_ImageMatrix(m_MorseWT,v_TFTime,v_FreqAxis,pv_Limits)

    %v_axes(idx,1) = subplot(4,2,v_data(idx));
    v_axes(idx,1) = subplot(4,1,idx);
    plot(v_RawData)
    hold on
    plot(v_FilData,'r')
    plot(v_Hil,'g')
    try
        xline(m_HFOEvents(:,1),'m')
        xline(m_HFOEvents(:,2),'k')
    catch
        continue
    end
     
end
linkaxes(v_axes,'x')
% xlim([7638170 7704110])
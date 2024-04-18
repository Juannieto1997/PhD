%% p_TimeFrequency
% Generate the time frequency análisis of the selected file. 
%% Initialize workspace 
%clear workspace
clc;close all; clear; 
%add functions folder to workspace 
addpath("C:\Users\Juan\Documents\CODES\FUNCTIONS")
%% Loading File to process
%look for file using GUI
[str_File,str_Path] = uigetfile('D:\');
% Select Folder to save Figures 
str_FigPath = uigetdir ('D:\Figures\') ;
% load file into matlab
fprintf(strcat('Loading File... ',str_File,'\n'))
load(fullfile(str_Path,str_File))
fprintf('File Loaded \n')
%% Preprocessing data 
% Original data Size
s_Size = size(m_Data);
% New Sampling Rate
s_nSRate = 1000;
% Data Filtering 
% Notch Filter (50Hz StopBand) - [45 55]
NotchFilt = designfilt('bandstopiir', ...       % Response type
       'PassbandFrequency1',45, ...    % Frequency constraints
       'StopbandFrequency1',47, ...
       'StopbandFrequency2',52, ...
       'PassbandFrequency2',55, ...
       'PassbandRipple1',1, ...         % Magnitude constraints
       'StopbandAttenuation',55, ...
       'PassbandRipple2',1, ...
       'DesignMethod','cheby2', ...      % Design method
       'MatchExactly','stopband', ...       % Design method options
       'SampleRate',s_SampRate)  ;            % Sample rate
% Beta Filtering (10 30 bandPass) - [5 40]

BetaFilt = designfilt('bandpassiir', ...
    'FilterOrder',20,'HalfPowerFrequency1',5, ...
    'HalfPowerFrequency2',40,'SampleRate',s_nSRate);

% Gama Filtering (30 100 bandPass) - [20 110]

GamaFilt = designfilt('bandpassiir', ...
    'FilterOrder',20,'HalfPowerFrequency1',20, ...
    'HalfPowerFrequency2',110,'SampleRate',s_nSRate);
% Matriz to store the new data 
m_PreData = zeros(s_Size(1),s_Size(2)/(s_SampRate/s_nSRate));
% Generate time vector 
v_time = [0:(s_Size(2)/(s_SampRate/s_nSRate)-1)]/s_nSRate;
% time limit for the Spectrogram
v_TimeLims = [0 v_time(end)];
% Frequency limits for the sectrogram 
v_Lims = [10 100];
% Define a tuekey window to remove border artifact 
v_tuk = tukeywin(length(v_time),0.1);
v_tuk = v_tuk';
% Baseline duration in minutes 
s_Baseline = 15;
v_LimEnd = [(v_time(end-1)-((2*s_Baseline)*60)) (v_time(end)-((s_Baseline)*60))];
b_End = logical((v_time>v_LimEnd(1) & (v_time<v_LimEnd(2))));
% Frequency Vector for FFT 
v_freq = linspace(0,s_nSRate,(s_Baseline*60*s_nSRate));
b_cut = logical((v_freq>v_Lims(1)) .* (v_freq<v_Lims(2)));
v_freqpl = v_freq(b_cut);
%% Preprocessing
for idx = 1:s_Size(1)
    %% Get Data
    % Channel Name 
    str_channel = stru_Header.ChannOrder{idx};
    c_Channel = split(str_channel,'_');
    %% Filtering Data
    % Print Messssage 
    f_Promsg(idx,s_Size(1))
    fprintf(strcat('Filtering Channel ... ',c_Channel{2},'\n'))
    % forward Filter to remove 50Hz noise 
    v_fitlData = filtfilt(NotchFilt, m_Data(idx,:));
    % Backwards filter to remove phase shift
    v_fitlData = filtfilt(NotchFilt,flip(v_fitlData));
    %% Validate Fitler using FFT 
    % v_FFT = fft(m_Data(idx,:));
    % v_FFT = v_FFT(1:floor(s_Size(2)/2));
    % v_freq = linspace(0,s_SampRate/2,floor(s_Size(2)/2));
    % plot(v_freq,abs(v_FFT))
    %% Subsample data
    m_PreData(idx,:) = resample(v_fitlData,s_nSRate,s_SampRate);
    %% Filtering Data for Display 
    % Beta Data 
    fprintf(strcat('Filtering data for Beta \n Filtering Channel ... ',c_Channel{2},'\n'))
    v_Beta = filtfilt(BetaFilt,m_PreData(idx,:));
    % gama Data 
    fprintf(strcat('Filtering data for Gama \n Filtering Channel ... ',c_Channel{2},'\n'))
    v_Gama = filtfilt(GamaFilt,m_PreData(idx,:));
    %% Wavelet 
    % Calculating Time Frequency analysis
    fprintf(strcat('Calculating Wavelet! Channel ... ',c_Channel{2},'\n'))
    fprintf('this is going to take quite some time ... you should go get a coffe \n')
    fprintf('I forgot i dont drink coffe ... \n')
    [m_MorseWT,v_TFTime,v_FreqAxis] = f_TimeFreq (m_PreData(idx,:),v_TimeLims,s_nSRate,v_Lims);
    fprintf('ok its done now \n')
    % Calculate limits 
    pv_Limits = [min(min(m_MorseWT)) max(max(m_MorseWT))]*0.06;
    %% Split data and canculate FFT 
    v_Baseline = m_PreData(idx,v_time<(s_Baseline*60));  % Baseline segment 
    v_FFTBaseline = fft(v_Baseline);% Calculate fourier transform 
    v_FFTBaseline = v_FFTBaseline(b_cut);
    [v_EnvBase,~] = envelope(abs(v_FFTBaseline));
    % Cut the last minutes of the recording
    
    v_End = m_PreData(idx,b_End); %last segment of the recording 
    v_FFTEnd = fft(v_End);
    v_FFTEnd = v_FFTEnd(b_cut);
    [v_EnvEnd,~] = envelope(abs(v_FFTEnd));
    %% Plotting data
    fig = figure(); 
    fig.Position= [1,41,1920,963];
    v_axis = zeros (1,3);
    % Plot Beta Activity
    v_axis(1) = subplot(4,5,[1 2 3]);
    plot(v_time,m_PreData(idx,:))
    hold on 
    plot(v_time,v_Beta.*v_tuk)
    ylim([-50 50])
    xlabel('time (s)'); 
    ylabel('Amplitude (mV)')
    title('Beta band (10 - 30 Hz)')
    % Plot Gama Activity 
    v_axis(2) = subplot(4,5,[6 7 8]);
    plot(v_time,m_PreData(idx,:))
    hold on 
    plot(v_time,v_Gama.*v_tuk)
    xlabel('time (s)'); 
    ylabel('Amplitude (mV)')
    title('Gama Band (30 - 100 Hz)')
     ylim([-50 50])
    % Time Frequency 
    v_axis([3,4]) =subplot(4,5,[11 12 13 14 15 16 17 18 19 20]);
    f_ImageMatrix(m_MorseWT,v_TFTime,v_FreqAxis,pv_Limits)
    xlabel('time (s)'); 
    ylabel('Frequency (Hz)')
    cb = colorbar('Location','eastoutside');
    % cb.Position = [.925 0.1095 0.0381 0.3762];
    linkaxes(v_axis,'x');
    xlim([v_time(1) v_time(end)])
    % Fourier transform
    subplot(4,5,[4 5 9 10])
    plot(v_freqpl,v_EnvEnd)
    hold on 
    plot(v_freqpl,v_EnvBase)
    xlabel('Frequency (Hz)')
    legend({'End','Baseline'})
    %ylim([1500 20000])
    xlim([v_freqpl(1) v_freqpl(end)])
    %% Save Figure 
    % Define figure name 
    str_FigName = strcat(str_File(1:10),c_Channel{2},'.fig');
    % Define Full folder name
    str_FigSPath = fullfile(str_FigPath,str_FigName);
    fprintf('Damn saving the figure takes way too long ... \n')
    fprintf('What else can I write in this messages ... ?  \n')
    fprintf('Do you hear the fan?? this thing is going to explode!!!! \n')
    %savefig(str_FigSPath)
    fprintf('Wait... is it done now?  \n')
    fprintf('Awesome Next file!!!!  \n')
    fprintf('WAIT NO! save it as a png \n')
    % Seriously? reusing code? kinda lasy dont you think? 
    str_FigName = strcat(str_File(1:10),c_Channel{2},'.png');
    str_FigSPath = fullfile(str_FigPath,str_FigName);
    saveas(fig,str_FigSPath)
    fprintf('can we start with the new fine finally? \n')
    fprintf('thank you kinda ... \n')
    close all;
end


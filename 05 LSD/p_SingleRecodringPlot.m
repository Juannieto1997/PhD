%% Prepare workspace
clc; close all; clear
%% Load Varibles 
load("D:\LSD\MCS_RAW\20250206\20250206S1000RawData.mat")
c_baseline  = c_cut;
load("D:\LSD\MCS_RAW\20250206\20250206S1001RawData.mat")
c_LSD  = c_cut;
load("D:\LSD\MCS_RAW\20250206\20250206S1002RawData.mat")
c_Washout  = c_cut;

str_chann = 'G8'; 

c_chann = st_header.ChannOrder; 

%% Filter Design
% Variables for filter
str_typeF = 'bandpassiir'; % badstopiir to remove 50Hz noice, IIR (not good for phase analysis)
v_FreqBand1 = [300 310]; % First cut of the filter
v_FreqBand2 = [3000 3010]; % Second cut of the filter 
s_RippleBand = 1; % ripple in db used the same for both ripples
s_Attenuation = 70; %dB 
str_Des = 'ellip';
s_MinToCrop = 0; % Minutes of the signal to crop for the data

d = designfilt(str_typeF, ...       % Response type
       'StopbandFrequency1',v_FreqBand1(1), ...    % Frequency constraints
       'PassbandFrequency1',v_FreqBand1(2), ...
       'PassbandFrequency2',v_FreqBand2(1), ...
       'StopbandFrequency2',v_FreqBand2(2), ...
       'StopbandAttenuation1',s_Attenuation, ...   % Magnitude constraints
       'PassbandRipple',s_RippleBand, ...
       'StopbandAttenuation2',s_Attenuation, ...
       'DesignMethod',str_Des, ...      % Design method
       'MatchExactly','passband', ...   % Design method options
       'SampleRate',s_SampRate);               % Sample rate

str_folder = 'C:\Users\Juan\Desktop\LSD Presentation';

for idx = 1:length(c_chann)
str_chann = c_chann{idx};
idxChan = find(contains(c_chann,str_chann),1);

m_baseline = c_baseline{idxChan,1};
m_LSD = c_LSD{idxChan,1};
m_washout = c_Washout{idxChan,1};

v_time = ([0:size(m_washout,2)-1]/s_SRate)*1000;
v_time = v_time-100;

v_axis = zeros(1,3);
fig = figure(); 

v_axis(1) = subplot(3,1,1);
plot(v_time,m_baseline,'b')
hold on 
plot(v_time,mean(m_baseline,1),'LineWidth',2,'Color','k')
title (strcat('Baseline: 125\muM Li N:',num2str(size(m_baseline,1)),str_chann),FontSize=20)
ylabel('Amplitude (\muV)')
ax= gca; 
ax.FontSize = 18;
v_axis(2) = subplot(3,1,2);
plot(v_time,m_LSD,'b')
hold on 
plot(v_time,mean(m_LSD,1),'LineWidth',2,'Color','k')
title (strcat('LSD treatment: 125\muM Li + 1\muM LSD N:',num2str(size(m_LSD,1))),FontSize=20)
ylabel('Amplitude (\muV)')
ax= gca; 
ax.FontSize = 18;
v_axis(3) = subplot(3,1,3);
plot(v_time,m_washout,'b')
hold on 
plot(v_time,mean(m_washout,1),'LineWidth',2,'Color','k')
title (strcat('Washout: 125\muM Li N:',num2str(size(m_washout,1))),FontSize=20)
ylabel('Amplitude (\muV)')
xlabel('Duration (ms)')
ax= gca; 
ax.FontSize = 18;
linkaxes(v_axis,'xy')
xlim([-30 100])
fig.Position = [1,49,2560,1.319333333333333e+03];
str_FileName = strcat('Figure',num2str(idx),'LFP.png');
saveas(gcf,fullfile(str_folder,str_FileName))
close all 
end 
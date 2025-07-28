%% p_PlotChansTemp

clc; close all; clear

% load("D:\LSD\MCS_RAW\20250207\20250207S1000RawData.mat")
% m_baseline  = m_Data;
% s_BEnd = (size(m_baseline,2)-1)/s_SampRate;
% fprintf('Baseline Loaded \n')
% load("D:\LSD\MCS_RAW\20250207\20250207S1001RawData.mat")
% m_LSD  = m_Data;
% s_LEnd = s_BEnd + (size(m_LSD,2)/s_SampRate);
% fprintf('LSD Loaded \n')
% load("D:\LSD\MCS_RAW\20250207\20250207S1002RawData.mat")
% m_Washout  = m_Data;
% fprintf('Washout Loaded \n')

load("D:\LSD\MCS_RAW\20241210\20241210S20000RawData.mat")
m_baseline  = m_Data;
s_BEnd = (size(m_baseline,2)-1)/s_SampRate;
fprintf('Baseline Loaded \n')
load("D:\LSD\MCS_RAW\20241210\20241210S20001RawData.mat")
m_LSD  = m_Data;
s_LEnd = s_BEnd + (size(m_LSD,2)/s_SampRate);
fprintf('LSD Loaded \n')
load("D:\LSD\MCS_RAW\20241210\20241210S20002RawData.mat")
m_Washout  = m_Data;
fprintf('Washout Loaded \n')

c_channels = stru_Header.ChannOrder;
m_Data = [m_baseline(1,:),m_LSD(1,:),m_Washout(1,:)]; 
v_time = [0:size(m_Data,2)-1]/s_SampRate;


for idxChan = 1:length(c_channels)
    str_chan = split(c_channels{idxChan},'_');
    str_chan = str_chan{end}; 
    str_title = strcat('Chan No: ',num2str(idxChan),' Chan: ',str_chan);
    v_Data = [m_baseline(idxChan,:),m_LSD(idxChan,:),m_Washout(idxChan,:)];
    close all 
    figure()
    plot(v_time,v_Data)
    hold on 
    xline(s_BEnd,LineWidth=2,Color='g')
    xline(s_LEnd,LineWidth=2,Color='g')
    title(str_title) 
    pause(2)
    
end 



function [s_eFPTime,s_eFPAmp, s_delay,s_beg,s_end,s_Duration,s_Slope]= f_FindParamsMat (m_Data,v_time,v_Baseline,s_Artifact,s_eFP)
%% f_FindParamsMat 
% function to work identify the parameters of any channel of the MEA
% recording 
% INPUTS 
% m_Data        : Matrix containing the channel to analize 
% v_time        : time vector of the data in ms 
% v_Baseline    : points to take the baseline in ms
% s_Artifact    : estimated poinf of the stimulation artifact in ms 
% s_eFP         : estimated time of the eFP in ms 
% OUTPUTS 
% --------- eFP ----------------
% s_eFPTime     : when the event field potential was found 
% s_eFPAmp      : Amplitude of the eFP relatively to the baseline 
% --------- Delay ----------------
% s_delay       : Delay between the artifact and the eFP 
% --------- Duration ----------------
% s_beg         : Beggining of the eFP 
% s_end         : End of the eFP
% s_Duration    : Duration of the event
%--------- Slope ----------------
% s_Slope       : Slope of the data
%% initialize variables 
v_size = size(m_Data);
s_eFPTime = zeros([1,v_size(1)]);
s_eFPAmp = zeros([1,v_size(1)]);
s_delay = zeros([1,v_size(1)]);
s_beg = zeros([1,v_size(1)]);
s_end = zeros([1,v_size(1)]);
s_Duration = zeros([1,v_size(1)]);
s_Slope = zeros([1,v_size(1)]);
%% Analize individual channels
for idxEvt = 1:v_size(1)
    v_data = m_Data(idxEvt,:); % Current event of the data 
    [s_eFPTime(idxEvt),s_eFPAmp(idxEvt), s_delay(idxEvt),s_beg(idxEvt),s_end(idxEvt),s_Duration(idxEvt),s_Slope(idxEvt)]= ...
        f_FindParams(v_data,v_time,v_Baseline,s_Artifact,s_eFP);
end 
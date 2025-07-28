function [s_eFPTime,s_eFPAmp, s_delay,s_beg,s_end,s_Duration,s_Slope]= f_FindParams(v_Data,v_time,v_Baseline,s_Artifact,s_eFP)
%% f_FindParams 
% function to work identify the parameters of any channel of the MEA
% recording 
% INPUTS 
% v_Data        : vector containing the channel to analize 
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
%% Define initial variables 
s_Baseline = mean(v_Data(v_Baseline(1):v_Baseline(2))); % Baseline of the channel 
v_CentData = v_Data-s_Baseline; % Center data acording to the baseline
s_win = 7; % Window to anallize fot the eFP 
s_artCorr = 1.5; % Variable to remove the end of the artifact 
%% Calculate eFP 
if s_eFP-s_win < s_Artifact+s_artCorr
    s_lim = s_Artifact+s_artCorr;
else 
    s_lim =s_eFP-s_win;
end
[~,s_timeStart] = min(abs(v_time-(s_lim)));
[~,s_timeEnd] = min(abs(v_time-(s_eFP+s_win)));
v_Win = v_CentData (s_timeStart:s_timeEnd); % Window to anallize the for the maximal point  
v_Wintime = v_time(s_timeStart:s_timeEnd); % Time of the window to analize
[s_AmpMax,s_max] = max(abs(v_Win));
s_maxTime = v_Wintime(s_max); 
%% Calculate Delay 
s_Delay = s_maxTime-s_Artifact;
%% Calculate Start and end points 
[~,s_timeStart] = min(abs(v_time-(s_Artifact+s_artCorr)));
[~,s_timeEnd] = min(abs(v_time-(s_maxTime)));
v_Win = v_CentData (s_timeStart:s_timeEnd); % Window to anallize the for the begining of the Data  
v_Wintime = v_time(s_timeStart:s_timeEnd); % Time of the window to analize
[s_AmpMin,s_min] = min(abs(v_Win));
s_minTime = v_Wintime(s_min); 
%% Calculate slope 
s_Half = (s_AmpMax-s_AmpMin)/2;
v_Win = v_Win - s_Half;
[~,s_Half] = min(abs(v_Win));
s_halfTime = v_Wintime(s_Half);
s_Slp = (v_CentData(v_time==s_halfTime)-s_AmpMin)/(s_halfTime-s_minTime);
%% Calculate end points 
[~,s_timeStart] = min(abs(v_time-(s_eFP)));
v_Win = v_CentData (s_timeStart:end);
v_Wintime = v_time(s_timeStart:end);
v_Win = v_Win - s_AmpMin;
[~,s_End] = min(abs(v_Win));
s_End = v_Wintime(s_End);

%% Set the variables to send
s_eFPTime = s_maxTime;
s_eFPAmp = s_AmpMax;
s_delay = s_Delay;
s_beg = s_minTime;
s_end = s_End;
s_Duration = (s_end-s_beg);
s_Slope = s_Slp;








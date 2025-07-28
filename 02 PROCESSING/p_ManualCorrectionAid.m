%% p_ManualCorrectionAid 
% Temporal code to help witht the detection of the peaks when there is a
% lot of noise in the channels. 
%% Prepare workspace
% DO NOT clear all the variables since they are requiered to continue with
% the code (p_CutDataV2)
close all; 
s_IPI = 20; % Inter pulse interval in seconds. 
%% Plot initial detection for validation 
v_time = [0:length(v_data)-1]/s_SampRate; % Time in s 
figure
% Plot the channel and the found peaks for stimulation 
plot(v_time,v_data); 
hold on 
plot(v_time(v_locs),v_data(v_locs),'*r')
%% Manually update stimulation points 
s_first = find(v_time ==24.4488); % CHANGE to the time of first stim 
v_locs = [s_first:s_IPI*s_SampRate:length(v_data)]; % set all the peaks acording to IPI
plot(v_time(v_locs),v_data(v_locs),'*g') % Plpt the new locations

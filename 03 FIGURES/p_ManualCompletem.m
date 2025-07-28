%% p_ManualComplete
% Create the view of the a single file
clc; close all; clear; 
addpath C:\Users\Juan\Documents\GitHub\PhD\FUNCTIONS
str_file = "F:\Data\USEXP\MC_RAW\20240612\20240612S30004CutSegments.mat";
v_Rfreq = [45 55];

load(str_file)

% for idxc = 1:length(c_cut)
%     c_curr = c_cut{idxc}; 
%     v_size = size(c_curr); 
%     for idxs = 1:v_size(1)
%         c_curr(idxs,:) = f_FFTfilt (c_curr(idxs,:),s_SRate,v_Rfreq); 
%     end 
%     c_cut{idxc} = c_curr; 
% end 


f_MeanCompleteMatriz (c_cut,st_header.ChannOrder,'.\','20240426S2CutSegments',[-100 100])

idx = 1; 

c_CurrData = c_cut{idx,1};
v_time = (0:size(c_CurrData,2)-1)/s_SRate;
plot(v_time,c_CurrData,'Color',[95 15 64]/255)
hold on 
plot(v_time,mean(c_CurrData),'Color','k','LineWidth',3)

xlim([0 0.2])
%ylim([-20 120])

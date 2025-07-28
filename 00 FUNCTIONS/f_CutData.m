function [m_data,v_peaks] = f_CutData (v_Data,s_SampleRate,v_cut,str_type,v_add)
%% f_CutData
% Function by: Juan Nieto
% funtion to cut the data acording to a trigger or artifact in the signal 
% INPUT: 
%   v_Data = vector que contiene 
%   s_SampleRate = sample rate of the data
%   v_cut = time in samples to cut before and after format = [before after]
%   varaging
%       1 - trigger type 
%           - 'auto' - Automatic detection using a simple find 
%              peaks for the detection of the stimulation artifact 
%           - 'time' - sets a fixed time interval betweem stimuli 
%           - 'vector' - recieves a with the positions to cut 
%       2 - trigger details
%           - if auto = none validated
%           - if time
%               - vector containing the following data in order
%                 [first-element, time-interval]
%           - if vector
%               - vector of samples containing the positions of the stimuli
%% Validate variable inputs

if strcmp(str_type,'auto') % if there are less than 2 inputs the case is input
    v_peaks = f_detectPeaks(v_Data,s_SampleRate);
elseif strcmp(str_type,'time')  
    v_info = v_add;
    s_init = v_info(1);
    s_inter = v_info(2);
    v_peaks = (s_init:s_inter*s_SampleRate:length(v_Data));
elseif strcmp(str_type,'vector')  
    v_peaks = v_add;
end 
%% Cut the data acording to the index 
 [~,~,m_data] = f_cut (v_Data,v_peaks,v_cut);
clc; close all; clear; 
addpath C:\Users\Juan\Documents\GitHub\PhD\FUNCTIONS
%str_file = 'E:\Data\USEXP\MC_RAW\20250211\20250211S1016CutSegments.mat';
str_file = "F:\Data\USEXP\MC_RAW\20241204\20241204S10000CutSegments.mat";
str_chann = 'H6'; 

load(str_file)

c_chann = st_header.ChannOrder; 
idxChan = find(contains(c_chann,str_chann),1);

m_data = c_cut{idxChan,:}; 
v_time = (0:size(m_data,2)-1)/s_SRate;
v_time = (v_time *1000)-10; 
figure 

plot(v_time,m_data,Color='b')
title(str_chann); 
xlabel ('time (ms)');
ylabel ('amplitude (uV)') 
%xlim([-0.3 1])



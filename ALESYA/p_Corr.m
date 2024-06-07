%% Prepare workspace
% Clear workspace
clc; close all; clear;
% add functions to path 
addpath('C:\Users\Juan\Documents\GitHub\PhD\FUNCTIONS')
%% Load Data
load('E:\Data\ALESYA\Corr\20230623s1RawData.mat') 

%% Prepare Data

s_NoSegments = 10; % Number of segments to analize
% s_cut = 7531250; % point to cut the data cutting for normal data 

% Cut matriz 
%m_Data = m_Data(:,s_cut:end);
s_size = size(m_Data); % size of the matrix 
v_Rfreq = [45 55];
%% filter data 
for idxChan = 1:s_size(1)
    m_Data(idxChan,:) = v_FFTfilt (m_Data(idxChan,:),s_SampRate,v_Rfreq);
end 
%% Calculate correlation 
s_step = 10; % Steps to analize the data in seconds.
s_NoSteps = floor (s_size(2)/(s_step*s_SampRate));
s_start = 1; 
s_end = (s_step*s_SampRate);
m_corrAdd = zeros(s_size(1),s_size(1),s_NoSteps);
for idxStep = 1:s_NoSteps
    % Cut the matrix for analysis 
    m_temp= m_Data(:,s_start:s_end);
    % Calculate correlation 
    m_corr=corrcoef(m_temp');
    m_corrAdd(:,:,idxStep) = m_corr;
    H1=heatmap(m_corr);
    H1.ColorLimits= [0 1];
    str_name = strcat('Corretlation',num2str(idxStep),'.png');
    saveas(gcf,str_name)
    % Update values 
    s_start = s_end+1; 
    s_end = s_start+(s_step*s_SampRate);
end 
%% Processing for clustering 
m_











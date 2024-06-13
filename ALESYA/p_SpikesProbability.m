%% p_SpikesProbability 
%% Prepare workspace 
% clear workspace 
clc; close all; clear; 
% add functions to path 
addpath('C:\Users\Juan\Documents\GitHub\PhD\FUNCTIONS')
% Load File 
 load('E:\Data\ALESYA\Corr\20230623s1RawData.mat');
% Define General variables 
v_size = size(m_Data); % data size for channels and No Samps  
%% Preprocess Data
% Define Variables for filtering 
v_cutFreq = [45 55]; % frequencies to remove 
m_FilteredData = zeros(v_size);
% Define Variables for spike detection 
c_spikes = cell(1,v_size(1));
% parameters for detection 
v_FiltFreq = [50 5000]; % cur frequencies for filter 
s_int = 100;% Number of samples for validation of the maximum
s_PeakDistance = 0.1; % distance between peaks in secods 
s_Npeaks = []; % NUmber of peaks to detect 
% Variables for spike sorting 
v_Anot = [];
m_cutData = [];
% Processing Data  
for idx = 1:v_size(1)
    f_Promsg(idx,v_size(1))
    fprintf('Filtering data... \n')
    m_FilteredData(idx,:)= v_FFTfilt (m_Data(idx,:),s_SampRate,v_cutFreq);
    fprintf('Detecting spikes... \n')
    v_peaks = f_detectPeaks(m_FilteredData(idx,:),s_SampRate,v_FiltFreq,s_int,s_PeakDistance,s_Npeaks);
    v_peaks = unique(v_peaks); % Get only the unique values of the data 
    
    m_cut = f_CutData(m_FilteredData(idx,:),s_SampRate,[700 700],'vector',v_peaks);% Cut the matrix
    % add values to the storage variables 
    m_cutData = [m_cutData;m_cut];
    v_Anot = [v_Anot ones(1,length(v_peaks))*idx];
    % Calculate and show principal components of each channel
    m_coef = pca(m_cut');
    figure()
    scatter(m_coef(1,:),m_coef(2,:))
end 
%% Sorting 
m_Coeff = pca(m_cutData'); 
for idx=1:10
    figure()
    scatter3(m_Coeff(idx,:),m_Coeff(idx+1,:),m_Coeff(idx+2,:))
end

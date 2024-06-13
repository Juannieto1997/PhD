%% p_ICA
%% Prepare workspace 
% clear workspace 
clc; close all; clear; 
% add functions to path 
addpath('C:\Users\Juan\Documents\GitHub\PhD\FUNCTIONS')
% Load File 
 load('E:\Data\ALESYA\Corr\20230623s1RawData.mat');
% Define General variables 
v_size = size(m_Data); % data size for channels and No Samps 
% Define groups 
v_G1 = [1:6];
v_G2 = [7:12]; 
v_G3 = [12:19];
c_groups = {v_G1;v_G2;v_G3};
s_Groups = 3; 
%% Preprocess Data
% Define Variables for filtering 
v_cutFreq = [45 55]; % frequencies to remove 
m_FilteredData = zeros(v_size);
for idx = 1:v_size(1)
    f_Promsg(idx,v_size(1))
    fprintf('Filtering data... \n')
    m_FilteredData(idx,:)= v_FFTfilt (m_Data(idx,:),s_SampRate,v_cutFreq);
end 
m_FilteredData = m_FilteredData(:,170000:end);
%% ICA Analysiis 
fprintf('computing independent component analysis... \n')
for idxG = 1:s_Groups
    f_Promsg(idxG,length(s_Groups))
    m_currGroup = m_FilteredData(c_groups{idxG},:);
    m_Feat = kICA(m_currGroup,length(c_groups{idxG}));
    % plot Filtered data + IC
    figure();
    v_axis = zeros(2,length(c_groups{idxG}));
    for idx = 1:length(c_groups{idxG})
        v_axis(idx*2-1) = subplot (length(c_groups{idxG}),2,idx*2-1);
        plot(m_currGroup(idx,:))
        v_axis(idx*2) = subplot (length(c_groups{idxG}),2,idx*2);
        try
            plot(m_Feat(idx,:))
        catch
            continue
        end

    end

    linkaxes(v_axis,'x')
end
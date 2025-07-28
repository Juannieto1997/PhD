%% p_ManualMatriz
% Updated version for the movie. 
clc; close all; clear; 
addpath C:\Users\Juan\Documents\GitHub\PhD\FUNCTIONS
str_file = "F:\Data\USEXP\MC_RAW\20250702\20250702S40001CutSegments.mat";
str_Img = "F:\Data\USEXP\MC_MCD\20250702\20250702S4.png";
%str_filename = strcat('20250211S1-FUS-D7-N',num2str(size(m_data,1)),'-800us-700mV-IPI15s-Long');
str_filename = strcat('20250702S40001');
%str_filename = 'water-Elec-stim';

m_bgd = imread(str_Img);
% m_bgd = imrotate(m_bgd,-90);
m_bgd = m_bgd(355:1780,280:1710,:);
%imshow(m_bgd)
load(str_file)
v_caxis = [-80 60];
s_FPS = 25; 
v_Rfreq = [45 55];
for idxc = 1:length(c_cut)
    c_curr = c_cut{idxc}; 
    v_size = size(c_curr); 
    % for idxs = 1:v_size(1)
    %     c_curr(idxs,:) = f_FFTfilt (c_curr(idxs,:),s_SRate,v_Rfreq); 
    % end 
    c_cut{idxc} = c_curr; 
end 



for idx = 1:length(c_cut)
    m_data = c_cut{idx,1};
    %m_av(idx,:) = mean(m_data(:,1:1000));
    Temp = m_data(1,1:1000);%-m_data(2,71);
    s_STD = std(Temp(1:20)); 
    Temp(logical((Temp < 10).*(Temp > -10))) = 0;
    m_av(idx,:) = Temp;
end 
%m_av = m_av(:,1:180);
v_size = size (m_av);

%str_filename = fullfile('C:\Users\Juan\Desktop\Temp',str_filename);

f_mkVid(m_av,m_bgd,st_header.ChannOrder,v_caxis,s_SRate,s_FPS,str_filename)
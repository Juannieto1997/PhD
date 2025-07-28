%% clear workspace 
clc;close all;clear;
%% read data
strName = 'E:\Data\ALESYA\3Brain\20240502d1S1_02.brw';
h5disp(strName)



s_NoChannels = 64*64;

s_SRate = h5readatt(strName,'/','SamplingRate'); 
s_DigMin = h5readatt(strName,'/','MinDigitalValue');
s_DigMax = h5readatt(strName,'/','MaxDigitalValue');
s_AnaMin = h5readatt(strName,'/','MinAnalogValue');
s_AnaMax = h5readatt(strName,'/','MaxAnalogValue'); 

v_data = h5read(strName,'/Well_A1/Raw');

s_len = length(v_data)/s_NoChannels;
v_time = [0:s_len-1]/s_SRate;

m_Data = reshape(v_data,[s_NoChannels,s_len]);
v_keep = (v_time>40);
v_keep2 = v_time<100;
v_keep = logical(v_keep .* v_keep2);
m_Data = m_Data(:,v_keep);
%%
v_Row = [10:2:40];
v_Col = [54:2:64];

v_r = linspace(0,0,length(v_Col));
v_g = linspace(0,0,length(v_Col));
v_b = linspace(0,255,length(v_Col))/255;

m_Color = [v_r;v_g;v_b];

s_len = length(m_Data(1,:));
v_time = [0:s_len-1]/s_SRate;

for idxRow = 1:length(v_Row)
    s_Row = v_Row(idxRow);
    for idxCol = 1:length(v_Col)
        s_Col = v_Col(idxCol); 
        s_idx = s_Row*64+s_Col;   
        v_Data = m_Data(s_idx,:); 
        v_correctedData = s_AnaMin + double(v_Data) .* ((s_AnaMax-s_AnaMin) ./ (s_DigMax-s_DigMin));
        plot(v_time,v_correctedData,"Color",m_Color(:,idxCol))
        hold on
    end 
    hold off
    str_saveName = strcat('3BrainDataRow',num2str(s_Row),'.png');
    saveas(gcf,str_saveName)
    close all
end 


v_chan = m_Data(s_idx,:);
v_correctedData = s_AnaMin + double(v_chan) .* ((s_AnaMax-s_AnaMin) ./ (s_DigMax-s_DigMin));
plot(v_time, v_correctedData)

% 
% for idx = 1:length(a)
% plot(v_time, v_correctedData(idx,:),'Color',m_Color(:,idx))
% hold on 
% end 


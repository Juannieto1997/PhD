%% Code_Video_MEA 
% Ivan original code for generating movie

clc; close all; clear;
% cd('')
% 
% textFiles = dir('*.txt');
% textFiles_cell = struct2cell(textFiles);
% numfiles = length(textFiles);
% mydata = cell(1, numfiles);
% 
% conc_data = [];
% peaks_bis = [];
% peaks_bis2 = [];
% peaks_t = [];
% peaks_t_max = [];
% peaks_ampl = [];
% peaks_ampl2 = [];
% peaks_ampl_max = [];
% resp_dur = [];
% ampl_dur = [];
% max_peak = [];
% 
% t_start = 59.71; 
% 
% t_end = 59.9;
% 
% t_start = 62.97; 
% 
% t_end = 62.97+0.04;
% 
load('E:\Data\USEXP\MC_RAW\20240321\20240321S20003CutSegments.mat')

tic 
%% organize header 
c_Chan = st_header.ChannOrder; 
for idx = 1:length(c_Chan)
    str_temp = c_Chan{idx}; 
    c_part = split(str_temp,' '); 
    str_temp = c_part{end}; 
    if length(str_temp) == 5
        str_temp(6) = str_temp(5);
        str_temp(5) = '0';
    end 
c_Chan{idx} = str_temp;
end 

[a,b] = sortrows(c_Chan');
%%

t_start = 0; 
% 
 t_end = 1;
% 
% t_start = 64.71; 
% 
% t_end = 64.9;
% 
 sampling_rate = 1/10;%In milliseconds 
% PRF = 50000;   %In Hz
% Period = 1/PRF*(1/sampling_rate*1000);
% 
% 
% t_window = single([t_start:sampling_rate*0.001:t_end]*(1/sampling_rate*1000));
% data_window = zeros(length(t_window),numfiles);

% for k = 1:numfiles 
% %   mydata{k} = imread(textFiles(k).name); 
%  delimiterIn = ' ';
%  headerlinesIn = 4;
%  myStructure = importdata(textFiles(k).name,delimiterIn,headerlinesIn);
%  data = myStructure.data;
%  time = 0:sampling_rate:sampling_rate*size(data,1)-sampling_rate; 
% 
% k  %display iteration number
% 
% 
% data_window(:,k) = data(t_window',size(data,2));
% %data_window(:,1) = [];
% 
% end
m_av = [];

for idx = 1:length(c_cut)
    m_data = c_cut{b(idx),1};
    %m_av(idx,:) = mean(m_data);
    m_av(idx,:) = m_data(1,:);
end 
m_av = m_av(:,1:180);
data_window = m_av';

% Define matriz for the data. 

elec_matrix = zeros(12,12,size(data_window,1));
elec_matrix(1,1,:) = 0;
elec_matrix(1,12,:) = 0;
elec_matrix(12,1,:) = 0;
elec_matrix(12,12,:) = 0;

%

cont = 1;

for i = 1:1
    
    for j = 4:size(elec_matrix,2)-3
        
elec_matrix(j,i,:) = data_window(:,j-3);
cont = cont +1; 
    
    end
    
end

for i = 2:2
    
    for j = 3:size(elec_matrix,2)-2
        
elec_matrix(j,i,:) = data_window(:,4+j);    
    
    end
    
end

for i = 3:3
    
    for j = 2:size(elec_matrix,2)-1
        
elec_matrix(j,i,:) = data_window(:,13+j);    
    
    end
    
end

    

for i = 4:size(elec_matrix,1)-3
    
    for j = 1:size(elec_matrix,2)
        
elec_matrix(j,i,:) = data_window(:,size(elec_matrix,2)*i-(size(elec_matrix,2)-j)-12);    
   
    end
    
end


for i = 10:10
    
    for j = 2:size(elec_matrix,2)-1
        
elec_matrix(j,i,:) = data_window(:,j+95);  
    
    end
    
end


for i = 11:11
    
    for j = 3:size(elec_matrix,2)-2
        
elec_matrix(j,i,:) = data_window(:,104+j);    
    
    end
    
end

for i = 12:12
    
    for j = 4:size(elec_matrix,2)-3
        
elec_matrix(j,i,:) = data_window(:,111+j);    
    
    end
    
end






% back_img = imread('slice3+Tx1-78MHz.png');
% crop_back = back_img(:,341:1920+341,:);

back_img = imread('E:\Data\USEXP\MC_MCD\20241204\20241204S2.png');
crop_back = back_img(280:1117,97:898,:);
crop_back = imresize(crop_back,[500 500]);
%crop_back = back_img;

n_frames = size(elec_matrix,3);
v = VideoWriter(sprintf('Test3'));
movie_frames = 1:size(elec_matrix,3);
v.FrameRate = 5;
open(v);
t_window_adj = 0:sampling_rate:(t_end-t_start)*1000+sampling_rate;

fi = figure('units','normalized','outerposition',[0 0 1 1]);

for i = 1:size(elec_matrix,3)
    
hold off
ibg2 = image(crop_back);
daspect([1 1 1]);
axis off
hold on
iim2 = imagesc(elec_matrix(:,:,i),'XData',[0 0+60*8],'YData',[0 0+60*8]);
ax = get(iim2,'Parent');
% cmap = linspace(min(min(m_av)),max(max(m_av)),800);
% caxis([min(min(m_av)) max(max(m_av))])
% ax.CLim = [[min(min(m_av)) max(max(m_av))]];
cmap = linspace(-30,60,800);
caxis([-30 60]);
ax.CLim = [-30 60];
%ax.CLim = [-150 120];
set(iim2,'AlphaData',0.5);
try
title(sprintf('%d msec',t_window_adj(i)),'fontsize',22)
MyColorMap = colormap(ax,jet(800));
MyColorMap(321,:) = 1;
colormap(ax,MyColorMap);
c = colorbar(ax);
% set(c, 'ylim', [min(min(m_av)) max(max(m_av))])
% c.Limits = [[min(min(m_av)) max(max(m_av))]]
set(c, 'ylim', [-30 60]);
c.Limits = [[-30 60]];
title(c,'Voltage (\muV)','fontsize',22);

% if US_frames(i) == 1
% a = annotation('textbox',[0.1 0.5 0.5 0.5],'String','US','EdgeColor','none','Color','r','fontsize',36);
% end

xlabel('mm','fontsize',22);
ylabel('mm','fontsize',22);
set(gca,'fontsize',22,'XTick',0.5:1:8.5,'XTickLabel',str2mat('0','0.2','0.4','0.6','0.8','1.0','1.2','1.4','1.6'),'YTick',0.5:1:8.5,'YTickLabel',str2mat('1.6','1.4','1.2','1.0','0.8','0.6','0.4','0.2','0'));
opengl('software');
frame = getframe(gcf);
writeVideo(v,frame);

%close all
catch 
    break 
end 

end
toc 
close all
close(v)

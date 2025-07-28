cd('')

textFiles = dir('*.txt');
textFiles_cell = struct2cell(textFiles);
numfiles = length(textFiles);
mydata = cell(1, numfiles);

conc_data = [];
peaks_bis = [];
peaks_bis2 = [];
peaks_t = [];
peaks_t_max = [];
peaks_ampl = [];
peaks_ampl2 = [];
peaks_ampl_max = [];
resp_dur = [];
ampl_dur = [];
max_peak = [];

t_start = 59.71; 

t_end = 59.9;

t_start = 62.97; 

t_end = 62.97+0.04;

t_start = 59.71; 

t_end = 59.9;

t_start = 64.71; 

t_end = 64.9;

sampling_rate = 0.02;      %In milliseconds
PRF = 50000;   %In Hz
Period = 1/PRF*(1/sampling_rate*1000);


t_window = single([t_start:sampling_rate*0.001:t_end]*(1/sampling_rate*1000));
data_window = zeros(length(t_window),numfiles);

for k = 1:numfiles 
%   mydata{k} = imread(textFiles(k).name); 
 delimiterIn = ' ';
 headerlinesIn = 4;
 myStructure = importdata(textFiles(k).name,delimiterIn,headerlinesIn);
 data = myStructure.data;
 time = 0:sampling_rate:sampling_rate*size(data,1)-sampling_rate; 
  
k  %display iteration number


data_window(:,k) = data(t_window',size(data,2));
%data_window(:,1) = [];
  
end


elec_matrix = zeros(12,12,size(data_window,1));
elec_matrix(1,1,:) = 0;
elec_matrix(1,12,:) = 0;
elec_matrix(12,1,:) = 0;
elec_matrix(12,12,:) = 0;
%
for i = 1:1
    
    for j = 2:size(elec_matrix,2)-1
        
elec_matrix(j,i,:) = data_window(:,j-1);    
    
    end
    
end




for i = 2:size(elec_matrix,1)-1
    
    for j = 1:size(elec_matrix,2)
        
elec_matrix(j,i,:) = data_window(:,size(elec_matrix,2)*i-(size(elec_matrix,2)-j)-2);    
   
    end
    
end


for i = 8:8
    
    for j = 2:size(elec_matrix,2)-1
        
elec_matrix(j,i,:) = data_window(:,size(elec_matrix,2)^2-size(elec_matrix,2)-3+j);;    
    
    end
    
end



% back_img = imread('slice3+Tx1-78MHz.png');
% crop_back = back_img(:,341:1920+341,:);

back_img = imread('Slice_Pic_2.png');
crop_back = back_img(588:1188,801:1601,:);


n_frames = size(elec_matrix,3);
v = VideoWriter(sprintf('MEA AMPLITUDE MOVIE SLICE BACKGROUND 1HZ 2nd pulse - 0810020'));
movie_frames = 1:size(elec_matrix,3);
v.FrameRate = 50;
open(v);
t_window_adj = 0:sampling_rate:(t_end-t_start)*1000+sampling_rate;


for i = 1:size(elec_matrix,3)
    
    
figure('units','normalized','outerposition',[0 0 1 1])
ibg2 = image(crop_back);
daspect([1 1 1])
axis off
hold on
iim2 = imagesc(elec_matrix(:,:,i),'XData',[155 155+60*8],'YData',[60 60+60*8]);
ax = get(iim2,'Parent')
cmap = linspace(-200,300,800);
caxis(ax,[-200 300])
ax.CLim = [-200 300];
set(iim2,'AlphaData',0.5);
title(sprintf('%d msec',t_window_adj(i)),'fontsize',22)
MyColorMap = colormap(ax,jet(800));
MyColorMap(321,:) = 1;
colormap(ax,MyColorMap);
c = colorbar(ax);
set(c, 'ylim', [-200 300])
c.Limits = [-200 300]
title(c,'Voltage (\muV)','fontsize',22)

% if US_frames(i) == 1
% a = annotation('textbox',[0.1 0.5 0.5 0.5],'String','US','EdgeColor','none','Color','r','fontsize',36);
% end

xlabel('mm','fontsize',22)
ylabel('mm','fontsize',22)
set(gca,'fontsize',22,'XTick',0.5:1:8.5,'XTickLabel',str2mat('0','0.2','0.4','0.6','0.8','1.0','1.2','1.4','1.6'),'YTick',0.5:1:8.5,'YTickLabel',str2mat('1.6','1.4','1.2','1.0','0.8','0.6','0.4','0.2','0'))
opengl('software')
frame = getframe(gcf)
writeVideo(v,frame);

close all

end


close(v)
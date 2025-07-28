function f_mkVid(m_Data,m_bgd,c_channels,v_caxis,s_srate,s_FPS,str_filename)
%% f_mkVid 
% Function to make the video of the data over time based from Ivans
% Code_Video_MEA.m
% by Juan Nieto (2024)
% INPUTS:
%   m_Data      - Matrix with the data to show as a plot
%   m_bgd       - 3d matrix containing the background image to plot
%   c_channels  - cell with all the channels, used for gridd construction
%           and validation of the order. 
%   v_caxis     - lims of the axix in the c direction [min max]
%   s_state     - sample rate of the m_Data matrix 
%   s_FPS       - sampling rate of the video used to define the speed 
%   str_filename- Name of the video to save as a file, must include the
%           folder. 
% OUTPUTS
%   NULL    
%% Prepare Data
% organize the channels in alphanumerical order
[c_channels, v_idx] = f_organizeChannels (c_channels);
% Recognize unique characters in channels
c_Let = {};
c_Num = {};
for idx = 1:length(v_idx)
    str_curr = c_channels{idx}; 
    str_let = str_curr(1); 
    str_num = str_curr(2:3);
    if  isempty(find(contains(c_Let,str_let)))
        c_Let{end+1} = str_let; 
    end 
    if  isempty(find(contains(c_Num,str_num)))
        c_Num{end+1} = str_num; 
    end 
end 

[c_Let,~] = sortrows(c_Let');
[c_Num,~] = sortrows(c_Num');
% Organize the matriz acording to the previous index 
m_Data = m_Data(v_idx,:);
s_size = size(m_Data);

%% Create Variables
v_time = ([0:s_size(2)-1]/s_srate)*1000; % Create time vector

% Initialize Video 
vi_Rec = VideoWriter(sprintf(str_filename)); % Create Video file
vi_Rec.FrameRate=s_FPS; % Set Framerate of the video 
open(vi_Rec);
% Initialize Matrix 
m_ElecData = zeros(length(c_Num),length(c_Let),s_size(2));
%Background figure
m_bgd = imresize(m_bgd,[500 500]);
%% Organize data in matrix 
for idxNum = 1:length(c_Num)
    for idxLet = 1:length(c_Let)
        str_currChann = strcat(c_Let{idxLet},c_Num{idxNum}); % Channel to look for 
        idxChan = find(contains(c_channels,str_currChann),1); 
        % Channel not found skip 
        if isempty(idxChan)
            continue
        end 
        m_ElecData(idxNum,idxLet,:) = m_Data(idxChan,:);
    end
end

%% Prepare Figure
% Figure 
fi = figure('units','normalized','outerposition',[0 0 1 1]);
for i = 1:size(m_ElecData,3)
    %% Plot background figure. 
    hold off
    ibg2 = image(m_bgd);
    daspect([1 1 1]);
    axis off
    hold on
    %% Plot data 
    iim2 = imagesc(m_ElecData(:,:,i),'XData',[25 475],'YData',[25 475]);
    ax = get(iim2,'Parent');
    % Adjust c axis
    cmap = linspace(v_caxis(1),v_caxis(2),800);
    caxis(v_caxis);
    ax.CLim = v_caxis;
    set(iim2,'AlphaData',0.5); % Set Transparency 
    %% Title and axis 
    title(sprintf('%d msec',v_time(i)),'fontsize',22) % Display time 
    MyColorMap = colormap(ax,jet(800));
    MyColorMap(321,:) = 1;
    colormap(ax,MyColorMap);
    c = colorbar(ax);
    set(c, 'ylim', v_caxis);
    c.Limits = [v_caxis];
    title(c,'Voltage (\muV)','fontsize',22);
    xlabel('mm','fontsize',22);
    ylabel('mm','fontsize',22);
    set(gca,'fontsize',22,'XTick',0.5:1:8.5,'XTickLabel',str2mat('0','0.2','0.4','0.6','0.8','1.0','1.2','1.4','1.6'),'YTick',0.5:1:8.5,'YTickLabel',str2mat('1.6','1.4','1.2','1.0','0.8','0.6','0.4','0.2','0'));
    %% Wirte Frame in video
    opengl('software');
    frame = getframe(gcf);
    writeVideo(vi_Rec,frame);
end

close all
close(vi_Rec)

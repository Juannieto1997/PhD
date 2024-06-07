%% p_plotMean
%generate the figures for the mean activity of each of the files 
%% Initialize workspace 
clc;close all; clear
addpath("C:\Users\Juan\Documents\CODES\FUNCTIONS")
%path fo save the figures
str_saveFolder = uigetdir("D:\Figures\");
% Define extention to look for 
str_ext = 'CutSegments.mat';
% Find all files
[str_path]=uigetdir("D:\Data\");
[c_Files] = f_SearchFiles (str_path,str_ext);

%% Open the files
for idx = 1:length(c_Files)
    % progress message
    f_Promsg(idx,length(c_Files))
    % Define the file names 
    str_folder = c_Files(idx).folder;
    str_file = c_Files(idx).name; 
    str_fileS = split(str_file,str_ext);
    str_fileS = str_fileS{1};
    fprintf(str_file)
    fprintf('\n')
    %load Data
    st_data = load(fullfile(str_folder,str_file));
    c_data = st_data.c_cut; 
    s_SRate = st_data.s_SRate; 
    c_chan = st_data.st_header.ChannOrder; 
    % Define time vector 
    v_time = ([0:size(c_data{1,1},2)-1]/s_SRate)*1000;
    %process each matrix 
    for idxChan = 1:length(c_chan)
        m_data = c_data{idxChan,1};
        str_chan = split(c_chan{idxChan},'_');
        str_chan = str_chan{2};
        strtitle = strcat(str_fileS,str_chan);
        % Calculate mean 
        v_mean = mean(m_data);
        % if isempty(m_data)
        %     continue 
        % end 
        % plot the data 
        try
        plot(v_time,m_data,'Color',[95 15 64]/255);
        hold on 
        plot(v_time,v_mean,'k','LineWidth',1);
        catch 
            continue 
        end 
        % add title and labels 
        ylabel('Amplitude (mV)');
        xlabel('Time (ms)')
        legend({'Raw','Average'})
        ylim([-40 90])
        xlim([-0 200])
        title(strtitle)
        %str_saveName = fullfile (str_saveFolder,strcat(strtitle,'.fig'));
        %savefig(str_saveName)
        % str_saveName = fullfile (str_saveFolder,strcat(strtitle,'.svg'));
        % print(gcf,'-vector','-dsvg',str_saveName)
        str_saveName = fullfile (str_saveFolder,strcat(strtitle,'.png'));
        saveas(gcf,str_saveName)
        close all
    end
end 
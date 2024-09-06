%% p_ProcessLTP
% code to plot the amplitude of the eFP of the different files that makes
% an LTP experimet. 
%% Preoare workspace
% Clear work space 
clc;close all; clear; 
% add functions to path 
addpath('C:\Users\Juan\Documents\GitHub\PhD\FUNCTIONS')
% look for folder 
str_DataPath = uigetdir('E:\Data','Select folder containing the data'); % Folder containing the data to read 
str_FigPath = uigetdir('E:\Figures','Select folder to store the figures'); % Folger to store the data into
%--------------------------------------------------------------------------
% String variables to identify the documents 
%--------------------------------------------------------------------------
str_ext = '*eFP*'; % Extension to look files 
str_savExt = 'ampvtime'; % Extension to save the figures 
%--------------------------------------------------------------------------
% Numerical Values of the experiment
%--------------------------------------------------------------------------
s_IPI = 20; % Inter Pulse interval in seconds 


%% Arrage files to work with 
% look for all the files containing eFP features. 
st_files = dir(fullfile(str_DataPath,str_ext));
c_names = {st_files.name}';
c_split = split(c_names,'.mat');
c_slices = {c_split{:,1}};
c_channels = {c_split{:,2}};

% Get the unique values of the slices and the channels
%NOTE: it is assumed that the recording is starting in the first file, esle
%it will not detect the first file of the recording and skip the whole
%slice. 
[c_unSlices,~,~] = unique(c_slices);
c_unSlices =c_unSlices(strlength(c_unSlices)<=13);
[c_unChann,~,~] = unique (c_channels);

%% Data processing 

for idxSlice = 1:length(c_unSlices)
    % Find all the files corresponding to the current slice
    fprintf('Slice Progress')
    f_Promsg(idxSlice,length(c_unSlices))
    str_currSlice = c_unSlices{idxSlice}; % Get current slice to work with 
    str_currSlice = str_currSlice(1:10);
    [b_Slice] = contains(c_slices,str_currSlice); % find all files corresponding to this slice
    for idxChannel = 1:length(c_unChann)
        fprintf('Channel Progress')
        f_Promsg(idxChannel,length(c_unChann))
        str_currChann = c_unChann{idxChannel}; 
        b_chann = contains(c_channels,str_currChann);
        % Define the files of interest (same slice, same channel)
        b_int = b_Slice.*b_chann;
        s_Nofiles = sum(b_int); 
        % if there is only one file skip the recording 
        if s_Nofiles <= 1
            fprintf(strcat('Current file, ',str_currSlice,' Does not have enough Recordings \n',...
                'Continue to next file'))
            break
        end 
        c_currNames = {c_names{logical(b_int)}};
        c_currNames = circshift(c_currNames,[1]);
        v_Features = f_loadFeatures(str_DataPath,c_currNames);
        v_time = [0:length(v_Features)-1]*s_IPI;
        % PLOT FIGURE 
        fig = figure (); 
        scatter(v_time,v_Features)
        hold on 
        xline(600,'Color','g')
        % save figure 
        str_saveName = strcat (str_currSlice,str_currChann,str_savExt,'.png'); 
        str_saveName = fullfile(str_FigPath,str_saveName);
        saveas(fig,str_saveName)
        close all;
    end 

end 
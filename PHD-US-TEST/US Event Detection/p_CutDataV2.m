%% p_CutDataV2
% code to find the peaks of the Ultrasound an electrical stimulation then
% cut the data into segments for análisis. 
%% Initialize Variables
% Variables for files
% Select folder using explorer
[str_path]=uigetdir("D:");
% get all folders ending with "RawData.mat"
% Extetion to look for
str_ext = 'RawData.mat';
[c_Files] = f_SearchFiles (str_path,str_ext);
% New extention of save files 
str_svExt = 'CutSegments.mat';
% Variables for cutting
%Triger type
str_type = 'vector';
% length of the segment to cut in seconds 
v_cutl = [0.01 0.25];
% Inter Pulse Interval in Seconds 
s_IPI = 15; % Seconds 
% Parameters for peak detection 
v_FiltFreq = [100 4000]; % Frequency to filter for peak detection
s_int = 100;  % Window size to validate the peaks 
s_int2 = 500; 
s_PeakDistance = 5; % Minimum distance in seconds between peaks
%% Check all the files in the folder 
for idx = 1:length(c_Files)
    % Progress bar 
    f_Promsg (idx,length(c_Files))
    str_folder = c_Files(idx).folder; % folder of the file 
    str_filename = c_Files(idx).name; % file name
     % Generate a new name to store the cut data
    str_saveName = split(str_filename,str_ext);
     % Validate if the file is Electrical stimulation or US stimulation 
    % Electrical files will end in S1 
    % US Files will end in S1000x
    ch_valIdx = str_saveName{1}(end-1); 
    num = str2num(ch_valIdx); 
    if isempty(num)
        b_Elec = 1; 
    else 
        b_Elec = 0;
    end 
    str_saveName = strcat(str_saveName{1},str_svExt);
    % load the data
    str_FFile = fullfile(str_folder,str_filename); % Full name 
    st_Data  = load(str_FFile); % Loading the data into the workspace
    % Defininf Varaibles to work with 
    m_data = st_Data.m_Data;
    s_SRate = st_Data.s_SampRate;
    st_header = st_Data.stru_Header;
    % Deleting the generated structure 
    clear st_Data;
    % Print update message 
    fprintf('Removing all variables that are no longer usefull \n (see you can do it! remove the negativity) \n')
    s_size = size(m_data);
    % create cell to store the result matrix
    c_cut = cell(s_size(1),2);
    % Estimate the number of pulses 
    s_NP = round((s_size(2)/s_SRate)/s_IPI);
   % Select a channel to perform the peak detection 
   fprintf ('Selecting a channel for the peak detection, this really takes a long time ...  \n I hope it is workth the wait, I mean according to math it should \n')
   [~,I]= max(max(m_data'));
   fprintf('OMG, cutting the channel also takes a while, I really hope this takes minutes of the code... \n')
   v_DetetChan = m_data(I,:);
   fprintf('Should be done now \n')
   % Peak detection in the channel with the most activity.
   v_peaks = f_detectPeaks(v_DetetChan,s_SRate,v_FiltFreq,s_int,s_PeakDistance,s_NP);
   if b_Elec
   else
        v_height = v_DetetChan(v_peaks); 
        [~,I] = max(v_height.^2);
        I = v_peaks(I);
        v_low = [I:-s_IPI*s_SRate:0];
        v_hi = [I:s_IPI*s_SRate:length(v_DetetChan)];
        v_peaks = [v_low v_hi];
        v_locs = unique(v_peaks);
        v_peaks = zeros(size(v_locs));
        % peak correction
        for idxP = 1:length(v_locs)
            s_peak = v_locs(idxP);
            if ((s_peak-s_int2)<0) + (s_peak+s_int2>length(v_DetetChan))
                continue
            end
            v_cut = v_DetetChan(s_peak-s_int2:s_peak+s_int2);
            % Find the first peak of the data
            v_max = abs(v_cut);
            [~,v_LocPeak] = max(v_max);
            % s_LocM = mean(v_max);
            % s_LocS = std(v_max);
            % [~,v_LocPeak] = findpeaks(v_max,'NPeaks',1,'MinPeakHeight',s_LocM+s_LocS);
            % Correct positon in the location vector
            v_peaks(idxP) = v_locs(idxP) + (v_LocPeak-s_int2)+1;
        end

   end
   % Initialize cell to store data
   c_cut = cell(s_size(1),2);
   % Cut all the data according to the peaks found in the channel with the
   % most activity
   fprintf('Cutting all channels (Yes Finally although you guess it it is going to be long ) \n')
    for idxChan = 1:s_size(1)
        % Progress bar 
        f_Promsg (idxChan,s_size(1))
        v_data = m_data(idxChan,:);
        [c_cut{idxChan,1},c_cut{idxChan,2}] = f_CutData (v_data,s_SRate,v_cutl.*s_SRate,str_type,v_peaks);
    end 
    fprintf('Saving data \n')
    save(fullfile(str_folder,str_saveName),'c_cut','s_SRate','st_header','-v7.3')
end 

function v_peaks = f_detectPeaks(v_data,s_SampRate,varargin)
%% f_detectPeaks
% Function to automaticalle detect the peaks from the gived data. 
%    INPUTS: 
%       v_data = vector containing the data to analize it will detect artifacts
%       in the elements 
%       s_SampRate = sample rate of the data 
%       varargin 
%           1 -  v_FiltFreq = Frequencies to filter the data to remove
%           artifacts 
%           2 - s_int = interval in samples for validation 
%           3 - s_PeakDistance = peak distances in seconds in seconds  
%           4 - s_Npeaks = Number of peaks to detect 
%% Define preset variables 
% Cut frequencies for the FFT filters 
if nargin == 2
    v_FiltFreq = [50 4000];
    s_int = 100;
    s_Npeaks = 40;
    s_PeakDistance = 1*s_SampRate;
elseif nargin == 3
    v_FiltFreq = varargin{1};
    s_int = 100;
    s_Npeaks = 40;
    s_PeakDistance = 1*s_SampRate;
elseif nargin == 4
    v_FiltFreq = varargin{1};
    s_int = varargin{2};
    s_Npeaks = 40;
    s_PeakDistance = 1*s_SampRate;
elseif nargin == 5
    v_FiltFreq = varargin{1};
    s_int = varargin{2};
    s_PeakDistance = varargin{3}*s_SampRate;
    s_Npeaks = 40;
elseif nargin == 6
    v_FiltFreq = varargin{1};
    s_int = varargin{2};
    s_PeakDistance = varargin{3}*s_SampRate;
    s_Npeaks = varargin{4};
end

%% preprocessing data
%fast fourier filtering 
% Keep only the higgest frequencies (should corresponds to the artifacts)
v_Filter = f_FFTfilt (v_data,s_SampRate,v_FiltFreq);
% plot(v_Filter)
% Difference of the data;
%v_diff = diff(v_Filter);
% apply a cuadratic function to the difference 
v_cuad = v_Filter.^2;
% calculate standard deviation 
s_STD = std(v_cuad)*3;
% calculate mean deviation 
s_mean = mean(mean(v_cuad));
%% Peak detection 
[~,v_locs] = findpeaks(v_cuad,'MinPeakHeight',s_mean+s_STD,'MinPeakDistance',s_PeakDistance,'NPeaks',s_Npeaks,'SortStr','descend');
% initialize final position for locations
v_peaks = zeros(size(v_locs));
% peak correction
for idx = 1:length(v_locs)
    s_peak = v_locs(idx);
    if ((s_peak-s_int)<0) + (s_peak+s_int>length(v_data)) 
        continue 
    end 
    v_cut = v_data(s_peak-s_int:s_peak+s_int);
    % Find the first peak of the data
    v_max = (v_cut).^2; 
    [~,v_LocPeak] = max(v_max);
    % s_LocM = mean(v_max);
    % s_LocS = std(v_max);
    % [~,v_LocPeak] = findpeaks(v_max,'NPeaks',1,'MinPeakHeight',s_LocM+s_LocS);
    % Correct positon in the location vector 
    v_peaks(idx) = v_locs(idx) + (v_LocPeak-s_int)-1;
end 




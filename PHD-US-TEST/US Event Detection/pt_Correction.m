idx = 3; 
v_Data = v_DetetChan;
plot(v_Data);
s_x = s_x.Position(1);
s_time = 15;
v_peaks = [s_x:(s_time*s_SRate)+32:length(v_Data)];
plot(v_Data)
hold on 
plot(v_peaks,v_Data(v_peaks),'r*')

v_locs = v_peaks;
s_int = 100;

v_peaks = zeros(size(v_locs));
% peak correction
for idx = 1:length(v_locs)
    s_peak = v_locs(idx);
    if ((s_peak-s_int)<0) + (s_peak+s_int>length(v_Data)) 
        continue 
    end 
    v_cut = v_Data(s_peak-s_int:s_peak+s_int);
    % Find the first peak of the data
    v_max = diff(v_cut).^2; 
    [~,v_LocPeak] = max(abs(v_cut));
    % s_LocM = mean(v_max);
    % s_LocS = std(v_max);
    % [~,v_LocPeak] = findpeaks(v_max,'NPeaks',1,'MinPeakHeight',s_LocM+s_LocS);
    % Correct positon in the location vector 
    v_peaks(idx) = v_locs(idx) + (v_LocPeak-s_int)+1;
end 
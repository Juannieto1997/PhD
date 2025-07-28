function [v_Peaks] = f_WinDetect(v_Data,s_nSTD,s_Win,s_dis,s_width)
v_Peaks = [];
s_nWin = ceil(length(v_Data)/s_Win); % Number of windows to use 
s_start = 1; 
s_end = s_start + s_Win;


for idx = 1:s_nWin
    if s_end > length(v_Data)
        s_end = length(v_Data); 
    end 
    v_cWin = v_Data(s_start:s_end);
    s_mean = mean(v_cWin); 
    s_STD = std(v_cWin); 
    s_thr = s_mean+(s_nSTD*s_STD); 
    try
        [~,v_locs] = findpeaks(v_cWin,'MinPeakHeight',s_thr,'MinPeakDistance',s_dis,'MinPeakWidth',s_width);
    catch
        continue 
    end 
    v_Peaks = [v_Peaks v_locs+s_start+1];
    % UPDATE THE CUTPOINTS! you almost forgot dude
    s_start = s_end;
    s_end = s_end+s_Win;
end 

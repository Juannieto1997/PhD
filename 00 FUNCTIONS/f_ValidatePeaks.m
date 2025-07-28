function v_peaks = f_ValidatePeaks(v_Data,v_peaks,s_Ini)
v_peaksF = zeros(length(v_peaks)+1);
for idx = 1:length(v_peaks)
    s_peak = v_peaks(idx);
    % Remove peaks that will cause problems
    if ((s_peak-s_Ini)<=0) + (s_peak+s_Ini>length(v_Data)) 
        continue
    end 
    s_start = s_peak-s_Ini; 
    s_end = s_peak+s_Ini; 
    v_WinVal = v_Data(s_start:s_end); 
    v_cuad = v_WinVal.^2; 
    [~,v_LocPeak] = max(v_cuad);
    v_peaksF(idx) = s_peak + (v_LocPeak-s_Ini)-1;
end 
v_peaks = unique(v_peaksF); 
v_peaks(1)=[];
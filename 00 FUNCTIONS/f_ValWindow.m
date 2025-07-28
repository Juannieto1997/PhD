function v_peaks = f_ValWindow(v_data,m_wins,s_off)
v_peaks = [];
cont = 1; 
for idx = 1:size(m_wins,2)
    try
        v_currWin = v_data(m_wins(1,idx)-s_off:m_wins(2,idx)+s_off);
    catch 
        continue 
    end 
    [~,s_max] = max((v_currWin.^2));
    v_peaks(cont) = s_max + m_wins(1,idx) - 1;

    cont = cont +1;
end 
v_peaks = unique (v_peaks);
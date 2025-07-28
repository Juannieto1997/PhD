function m_cross=f_zeroCross(v_Data)
v_PosPeaks = find(diff(v_Data>0));
v_PosPeaks2 = find(diff(v_Data(length(v_Data):-1:1)>0));
v_PosPeaks2 = length(v_Data)-v_PosPeaks2+1;
v_PosPeaks2= v_PosPeaks2(length(v_PosPeaks2):-1:1);
m_cross = [v_PosPeaks(1:end-1);v_PosPeaks2(2:end)];


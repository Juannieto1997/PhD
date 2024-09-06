v_time = [0:length(v_data)-1]/s_SampRate;
figure
plot(v_time,v_data); 
hold on 
plot(v_time(v_locs),v_data(v_locs),'*r')
s_first = find(v_time ==1.2580);
v_locs = [s_first:5*s_SampRate:length(v_data)];
plot(v_time(v_locs),v_data(v_locs),'*g')
plot(v_time(v_peaks),v_data(v_peaks),'*g')
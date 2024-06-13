v_time = [0:length(v_data)-1]/s_SampRate;
figure
plot(v_time,v_data); 
hold on 
plot(v_time(v_locs),v_data(v_locs),'*r')
s_first = find(v_time == 11.5802);
v_locs = [s_first:15*s_SampRate:length(v_data)];
plot(v_time(v_locs),v_data(v_locs),'*g')
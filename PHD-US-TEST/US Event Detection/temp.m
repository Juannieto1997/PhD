%% Energy plot
v_time = linspace (0,9,10000); % Time in ms
s_222MHz_700mV = 26; % (W) Power
v_222MHz_700mV = (v_time/1000)*s_222MHz_700mV; % (J) Energy
plot(v_time,v_222MHz_700mV,'--')
xlabel('Time (ms)',FontSize=20)
ylabel('Energy (J)',FontSize=20)
hold on 

s_523MHz_700mV = 7.91;
s_maxtime = v_222MHz_700mV(end)/s_523MHz_700mV;
v_time = linspace(0,s_maxtime*1000,10000); 
v_523MHz_700mV = s_523MHz_700mV * (v_time/1000);
plot(v_time,v_523MHz_700mV)

s_523MHz_800mV = 9.18;
s_maxtime = v_222MHz_700mV(end)/s_523MHz_800mV;
v_time = linspace(0,s_maxtime*1000,10000); 
v_523MHz_800mV = s_523MHz_800mV * (v_time/1000);
plot(v_time,v_523MHz_800mV)

s_822MHz_700mV = 3.5;
s_maxtime = v_222MHz_700mV(end)/s_822MHz_700mV;
v_time = linspace(0,s_maxtime*1000,10000); 
v_822MHz_700mV = s_822MHz_700mV * (v_time/1000);
plot(v_time,v_822MHz_700mV)

s_822MHz_700mV = 4.09;
s_maxtime = v_222MHz_700mV(end)/s_822MHz_700mV;
v_time = linspace(0,s_maxtime*1000,10000); 
v_822MHz_700mV = s_822MHz_700mV * (v_time/1000);
plot(v_time,v_822MHz_700mV)



legend({'2.22MHz 700mV','5.23MHz 700mV','5.23MHz 800mV','8.22MHz 700mV','8.22MHz 800mV'})

str_name = 'E:\Data\ALESYA\3Brain\20240416s1_08.brw';

h5disp (str_name)
s_SRate = h5readatt(str_name,'/','SamplingRate');
s_MinDig = h5readatt(str_name,'/','MinDigitalValue');
s_MinAna = h5readatt(str_name,'/','MinAnalogValue');
s_MaxDig = h5readatt(str_name,'/','MaxDigitalValue');
s_MaxAna = h5readatt(str_name,'/','MaxAnalogValue');


v_data = h5read(str_name,'/Well_A1/Raw');
v_data = s_MinAna+(double(v_data)*((s_MaxAna-s_MinAna)/(s_MaxDig -s_MinDig)));
s_Nchan = 64*64;
s_l = length(v_data)/s_Nchan;
m_data = reshape(v_data,[s_Nchan s_l]);

v_temp = m_data(1098:1128,:);
v_size = size(v_temp);
v_b = round(linspace(0,255,v_size(1)));
v_r = zeros (1,v_size(1));
v_g = zeros (1,v_size(1));
m_colors = single([v_r' v_g' v_b']/255);
v_time = [0:s_l-1]/s_SRate; 

ax = plot(v_time,v_temp); 
ax = gca;
ax.ColorOrder = m_colors;
xlabel('time (s)')
ylabel('Amplitude (uV)')
clc; close all; clear; 
str_dir = 'C:\Users\Juan\Desktop\Temp';
st_files = dir(str_dir);
st_files = st_files(3:end);
v_amp = [];
for idx = 1:length(st_files)
    str_currFile = st_files(idx).name;
    v_temp = load(fullfile(str_dir,str_currFile));
    v_amp = [v_amp v_temp.v_eFPAmp];
end 
v_time = linspace(0,length(v_amp)*20,length(v_amp));
scatter(v_time,v_amp)
xlabel ('time (s)')
ylabel('Amplitude (uV)')
xline(1800,'Color','g')

fig = gcf;
ax = gca;
ax.FontSize = 15;
set(fig, 'color', 'none');    
set(ax, 'color', 'none');
exportgraphics(fig,fullfile(str_dir,'LTP.eps'),...   % since R2020a
    'ContentType','vector',...
    'BackgroundColor','none')
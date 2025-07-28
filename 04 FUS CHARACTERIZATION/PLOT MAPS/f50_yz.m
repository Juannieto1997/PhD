clear all;
close all;

%% Paramètres
cd('D:\MUTATION\Champs de Pression\TETE1\jeudi 21\focale 50\f50 yz\');
fileformat = 'scanTestZ%03dY%03d.txt';
focal = 50;

% nombre de lignes et colonnes acquises
n_row = 101; % <-- deuxieme variable
n_col = 51; % <-- premiere variable

hydrophoneSensitivity = 4.2e-7; % V/Pa

% 1ere ligne de données du fichier texte
dataline = 81;

% plan Yx
z0 = 0;
y0 = 50;
dz = 0.2;
dy = 0.2;

%% lecture des fichiers

n = 1;

min_array = zeros(n_row, n_col);
max_array = zeros(n_row, n_col);
rms_array = zeros(n_row, n_col);

r_array = zeros(n_col,1);
z_array = zeros(n_row,1);
pmin_array = zeros(n_col*n_row,3);
pmax_array = zeros(n_col*n_row,3);
prms_array = zeros(n_col*n_row,3);

z = z0;
y = y0;

for i=0:1:n_row-1
    for j=0:1:n_col-1
      z = z+dz;
      n_in = sprintf(fileformat,i,j); 
      [d1,d2,d3] = textread(n_in, '%f %f %f', 'headerlines', dataline);
      pmin = min(d2);
      pmax = max(d2);
      prms = sqrt(mean(d2.^2));
      Vmin_array(i+1,j+1) = min(d2);
      Vmax_array(i+1,j+1) = max(d2);
      Vrms_array(i+1,j+1) = prms;
      pmin_array(n,:) = [pmin,z,y];
      pmax_array(n,:) = [pmax,z,y];
      prms_array(n,:) = [prms,z,y];
      n = n+1;
    end
    100*n/(n_row*n_col)
    z = 0;
    y = y+dy;
end

% Conversion en MPa
min_array = Vmin_array ./ hydrophoneSensitivity * 1e-6;
max_array = Vmax_array ./ hydrophoneSensitivity * 1e-6;
rms_array = Vrms_array ./ hydrophoneSensitivity * 1e-6;

%% affichage

% axes
z = (-floor(n_col/2)*dz:dz:floor(n_col/2)*dz) + z0;
y = (-floor(n_row/2)*dy:dy:floor(n_row/2)*dy) + y0;

width = 908;     % Width in pixels
height = 747;    % Height in pixels
alw = 0.75;    % AxesLineWidth
fsz = 20;      % Fontsize
lw = 1.5;      % LineWidth
msz = 8;       % MarkerSize
fig1 = figure(13);
fig1.Position = [520 51 908 747];
% pos=get(gcf,'Position');
% set(gcf, 'Position', [pos(1) pos(2) (width) (height)]); %<- Set size
set(gca, 'FontSize', fsz); %<- Set properties
set(gcf,'InvertHardcopy','on');
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- (width))/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left-1, bottom, width+1, height];

% affichage max
fig2 = figure(2);
fig2.Position = [520 51 908 747];
pos=get(gcf,'Position');
set(gcf, 'Position', [pos(1) pos(2) (width) (height)]);
set(gca, 'FontSize', fsz); %<- Set properties
set(gcf,'InvertHardcopy','on');
[cv,ch]=contourf(z,y,(max_array),80);
set(ch,'edgecolor','none');
axis equal
axis([-floor(n_col/2)*dz+z0 floor(n_col/2)*dz+z0 -floor(n_row/2)*dy+y0 floor(n_row/2)*dy+y0])
axisLabel = get(gca,'YTickLabel');
set(gca,'YTickLabel',axisLabel(end:-1:1));
xlabel('Transducer length, \ity\rm (mm)','fontsize',20)
ylabel('US propagation axis, \itz\rm (mm)','fontsize',20)
set(gca,'fontsize',20)
title(' ')
c = colorbar;
ylabel(c,'Pressure, \itp_{pk} \rm(MPa)')
xticks([-5 -4 -3 -2 -1 0 1 2 3 4 5])
xticklabels({'-5','-4','-3','-2', '-1','0','1', '2','3','4', '5'})

frame = getframe(fig2)
im = frame2im(frame);
[A,map] = rgb2ind(im,256);
imwrite(A,map,'testPMax.tif','tiff');
savefig(fig2,'testPMax.fig')




% affichage min
fig3 = figure(14);
fig3.Position = [520 51 908 747];
pos=get(gcf,'Position');
set(gcf, 'Position', [pos(1) pos(2) (width) (height)]);
set(gca, 'FontSize', fsz); %<- Set properties
set(gcf,'InvertHardcopy','on');
[cv,ch]=contourf(z,y,(-min_array),80);
set(ch,'edgecolor','none');
axis equal
axis([-floor(n_col/2)*dz+z0 floor(n_col/2)*dz+z0 -floor(n_row/2)*dy+y0 floor(n_row/2)*dy+y0])
axisLabel = get(gca,'YTickLabel');
set(gca,'YTickLabel',axisLabel(end:-1:1));
xlabel('Transducer length, \ity\rm (mm)','fontsize',20)
ylabel('US propagation axis, \itz\rm (mm)','fontsize',20)
set(gca,'fontsize',20)
xticks([-5 -4 -3 -2 -1 0 1 2 3 4 5])
xticklabels({'-5','-4','-3','-2', '-1','0','1', '2','3','4', '5'})
title(' ')
c = colorbar;
ylabel(c,'Pressure, \itp_{pk-} \rm(MPa)')
[row,col] = find((min_array) == min(min_array(:)));

frame = getframe(fig3)
im = frame2im(frame);
[A,map] = rgb2ind(im,256);
imwrite(A,map,'testPMin.tif','tiff');
savefig(fig3,'testPMin.fig')




% affichage rms
fig4 = figure(15);
fig4.Position = [520 51 908 747];
pos=get(gcf,'Position');
set(gcf, 'Position', [pos(1) pos(2) (width) (height)]);
set(gca, 'FontSize', fsz); %<- Set properties
set(gcf,'InvertHardcopy','on');
[cv,ch]=contourf(z,y,(rms_array),80);
set(ch,'edgecolor','none');
axis equal
axis([-floor(n_col/2)*dz+z0 floor(n_col/2)*dz+z0 -floor(n_row/2)*dy+y0 floor(n_row/2)*dy+y0])
axisLabel = get(gca,'YTickLabel');
set(gca,'YTickLabel',axisLabel(end:-1:1));
xlabel('Transducer length, \ity\rm (mm)','fontsize',20)
ylabel('US propagation axis, \itz\rm (mm)','fontsize',20)
set(gca,'fontsize',20)
xticks([-5 -4 -3 -2 -1 0 1 2 3 4 5])
xticklabels({'-5','-4','-3','-2', '-1','0','1', '2','3','4', '5'})
title(' ')
c = colorbar;
ylabel(c,'Pressure, \itp_{RMS}\rm(MPa)','Fontsize', 20)
set(gca,'fontsize',20)

frame = getframe(fig4)
im = frame2im(frame);
[A,map] = rgb2ind(im,256);
imwrite(A,map,'testPrms.tif','tiff');
savefig(fig4,'testPrms.fig')

%% Plot Max en Echelle Log

max_array_norm = max_array/(max(max(max_array)));
max_array_dB = 20*log10(max_array_norm);

% axes
y = (-floor(n_col/2)*dz:dz:floor(n_col/2)*dz) + z0;
x = (-floor(n_row/2)*dy:dy:floor(n_row/2)*dy) + y0;

width = 908;     % Width in pixels
height = 747;    % Height in pixels
alw = 0.75;    % AxesLineWidth
fsz = 24;      % Fontsize
lw = 1.5;      % LineWidth
msz = 8;       % MarkerSize
fig1 = figure(13);
fig1.Position = [520 51 908 747];
% pos=get(gcf,'Position');
% set(gcf, 'Position', [pos(1) pos(2) (width) (height)]); %<- Set size
set(gca, 'FontSize', fsz); %<- Set properties
set(gcf,'InvertHardcopy','on');
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- (width))/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left-1, bottom, width+1, height];

% affichage max
fig5 = figure(2);
fig5.Position = [520 51 908 747];
pos=get(gcf,'Position');
set(gcf, 'Position', [pos(1) pos(2) (width) (height)]);
set(gca, 'FontSize', fsz); %<- Set properties
set(gcf,'InvertHardcopy','on');
[cv,ch]=contourf(y,x,(max_array_dB),80);
set(ch,'edgecolor','none');
hold on
[con_v con_h]  = contour(y,x,(max_array_dB),[-3,-3]);
set(con_h,'LineStyle','-')
set(con_h,'LineWidth',1)
set(con_h,'LineColor','k')
hold on
[con_v con_h]  = contour(y,x,(max_array_dB),[-6,-6]);
set(con_h,'LineStyle',':')
set(con_h,'LineWidth',1)
set(con_h,'LineColor','k')
axis equal
axis([-floor(n_col/2)*dz+z0 floor(n_col/2)*dz+z0 -floor(n_row/2)*dy+y0 floor(n_row/2)*dy+y0])
axisLabel = get(gca,'YTickLabel');
set(gca,'YTickLabel',axisLabel(end:-1:1));
xlabel('Transducer length, \ity\rm (mm)','fontsize',20)
ylabel('US propagation axis, \itz\rm (mm)','fontsize',20)
set(gca,'fontsize',20)
xticks([-5 -4 -3 -2 -1 0 1 2 3 4 5])
xticklabels({'-5','-4','-3','-2', '-1','0','1', '2','3','4', '5'})
max_p = max(max(max_array));
str = sprintf('p_{sp} = %.2f MPa', max_p)
title(str)
c = colorbar;
ylabel(c,'Pressure, \itp_{pk} \rm(dB)')

% set(gca, 'YDir', 'normal'); %<- Set properties
% if ~exist('PressureField','var')
% save('PressureField','max_array')
% else
% save('PressureField','max_array','-append')
% end

frame = getframe(fig5)
im = frame2im(frame);
[A,map] = rgb2ind(im,256);
imwrite(A,map,'testPMax_dB.tif','tiff');
savefig(fig5,'testPMax_dB.fig')


%% 2D Plot of Pressure at focus

max_pixel_index = find(max_array == max(max(max_array)),1,'first');
max_col = rem(max_pixel_index,size(max_array,2));
max_line = rem(max_pixel_index,size(max_array,1));
figure
plot(y,max_array(max_line,:))

figure
plot(x,max_array(:,max_col))


%% 1D Plot of Pressure at focus
max_array_norm = max_array/(max(max(max_array)));
max_array_dB = 20*log10(max_array_norm);

[max_p, max_pixel_index] =max(max_array(:))
[max_line max_col]=ind2sub(size(max_array),max_pixel_index)

fig8 = figure(8)
plot(y,max_array_dB(max_line,:))
fig8.Position = [520 51 908 747];
pos=get(gcf,'Position');
set(gcf, 'Position', [pos(1) pos(2) (width) (height)]);
str = sprintf('p_{sp} = %.2f MPa', max_p)
title(str)
xlabel('Transducer length, \ity\rm (mm)')
ylabel('Pressure, \itp_{pk} \rm(dB)')%(c,'Pressure, \itp_{pk} \rm(MPa)')
frame = getframe(fig8)
im = frame2im(frame);
[A,map] = rgb2ind(im,256);
imwrite(A,map,'Y.tif','tiff');
savefig(fig8,'Y.fig')

fig7 = figure(7)
plot(x,max_array_dB(:,max_col))
fig7.Position = [520 51 908 747];
pos=get(gcf,'Position');
set(gcf, 'Position', [pos(1) pos(2) (width) (height)]);
str = sprintf('p_{sp} = %.2f MPa', max_p)
title(str)
xlabel('US propagation axis, \itz \rm(mm)')
ylabel('Pressure, \itp_{pk} \rm(dB)')
frame = getframe(fig7)
im = frame2im(frame);
[A,map] = rgb2ind(im,256);
imwrite(A,map,'Z.tif','tiff');
savefig(fig7,'Z.fig')

%% Contour Plot
% 
% max_array_norm = max_array/(max(max(max_array)));
% max_array_dB = 20*log10(max_array_norm);
% 
% % axes
% y = (-floor(n_col/2)*dz:dz:floor(n_col/2)*dz) + z0;
% x = (-floor(n_row/2)*dy:dy:floor(n_row/2)*dy) + y0;
% 
% width = 908;     % Width in pixels
% height = 747;    % Height in pixels
% alw = 0.75;    % AxesLineWidth
% fsz = 24;      % Fontsize
% lw = 1.5;      % LineWidth
% msz = 8;       % MarkerSize
% %fig1 = figure(13);
% fig1.Position = [520 51 908 747];
% % pos=get(gcf,'Position');
% % set(gcf, 'Position', [pos(1) pos(2) (width) (height)]); %<- Set size
% set(gca, 'FontSize', fsz); %<- Set properties
% set(gcf,'InvertHardcopy','on');
% papersize = get(gcf, 'PaperSize');
% left = (papersize(1)- (width))/2;
% bottom = (papersize(2)- height)/2;
% myfiguresize = [left-1, bottom, width+1, height];



fig6 = figure(6)
fig6.Position = [520 51 908 747];
pos=get(gcf,'Position');
set(gcf, 'Position', [pos(1) pos(2) (width) (height)]);
[X,Y] = meshgrid(y,x)
[M,c] = contourf(X,Y,max_array_dB,[-12 -9 -6 -3],'Fill','on')
cmap = [ 239,255,255;  255,223,255; 191,191,255; 255,0,0]/255;
colormap(cmap);
colorbar()
set(gca, 'clim', [-12 0])
max_p = max(max(max_array));
str = sprintf('p_{sp} = %.2f MPa', max_p)
title(str)
xlabel('Transducer length, \ity\rm (mm)','fontsize',20)
ylabel('US propagation axis, \itz\rm (mm)')
xticks([-5 -4 -3 -2 -1 0 1 2 3 4 5])
xticklabels({'-5','-4','-3','-2', '-1','0','1', '2','3','4','5'})
yticks([40:1:60])
yticklabels({'60','','58','','56','','54','','52','','50','','48','','46','','44','','42','','40'})
c = colorbar;
c.YTick = [-12, -9, -6, -3, -0];
ylabel(c,'Pressure, \itp_{pk} \rm(dB)','fontsize',20)
set(gca,'fontsize',20)
str = sprintf('p_{sp} = %.2f MPa', max_p)
title(str)
daspect ([1,1,1])

frame = getframe(fig6)
im = frame2im(frame);
[A,map] = rgb2ind(im,256);
imwrite(A,map,'testPMax_dB_contour.tif','tiff');
savefig(fig6,'testPMax_dB_contour.fig')
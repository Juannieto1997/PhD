clear all;
close all;

%% Paramètres
cd('E:\Data\PRESURE MAPS\20240911\82MHz\xy');
fileformat = ['ScanY%03dX%03d.txt'];
focal = 50;

% nombre de lignes et colonnes acquises
n_row = 41; % <-- deuxieme variable
n_col = 41; % <-- premiere variable

hydrophoneSensitivity = 4.2e-7; % V/Pa

% 1ere ligne de données du fichier texte
dataline = 82;

% plan Yx
y0 = 0;
x0 = 0;
dy = 0.05;
dx = 0.05;

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

y = y0;
x = x0;

for i=0:1:n_row-1
    for j=0:1:n_col-1
      y = y+dy;
      n_in = sprintf(fileformat,i,j); 
      [d1,d2,d3] = textread(n_in, '%f %f %f', 'headerlines', dataline);
      %d2 = d2+0.0022;
      d2 = abs(d2);
      pmin = min(d2);
      pmax = max(d2);
      prms = sqrt(mean(d2.^2));
      Vmin_array(i+1,j+1) = min(d2);
      Vmax_array(i+1,j+1) = max(d2);
      Vrms_array(i+1,j+1) = prms;
      pmin_array(n,:) = [pmin,y,x];
      pmax_array(n,:) = [pmax,y,x];
      prms_array(n,:) = [prms,y,x];
      n = n+1;
    end
    100*n/(n_row*n_col)
    y = 0;
    x = x+dx;
end

% Conversion en MPa
min_array = Vmin_array ./ hydrophoneSensitivity * 1e-6;

max_array = Vmax_array ./ hydrophoneSensitivity * 1e-6;
rms_array = Vrms_array ./ hydrophoneSensitivity * 1e-6;

%% affichage

% axes
y = (-floor(n_col/2)*dy:dy:floor(n_col/2)*dy) + y0;
x = (-floor(n_row/2)*dx:dx:floor(n_row/2)*dx) + x0;

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
[cv,ch]=contourf(y,x,(max_array),80);
set(ch,'edgecolor','none');
axis equal
axis([-floor(n_col/2)*dy+y0 floor(n_col/2)*dy+y0 -floor(n_row/2)*dx+x0 floor(n_row/2)*dx+x0])
axisLabel = get(gca,'YTickLabel');
set(gca,'YTickLabel',axisLabel(end:-1:1));
xlabel('Transducer width, \itx\rm (mm)','fontsize',20)
ylabel('Transducer length, \ity\rm (mm)','fontsize',20)
set(gca,'fontsize',20)
title(' ')
c = colorbar;
ylabel(c,'Pressure, \itp_{pk} \rm(MPa)','fontsize',20)
set(gca,'fontsize',20)
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
[cv,ch]=contourf(y,x,(-min_array),80);
set(ch,'edgecolor','none');
axis equal
axis([-floor(n_col/2)*dy+y0 floor(n_col/2)*dy+y0 -floor(n_row/2)*dx+x0 floor(n_row/2)*dx+x0])
axisLabel = get(gca,'YTickLabel');
set(gca,'YTickLabel',axisLabel(end:-1:1));
xlabel('Transducer width, \itx\rm (mm)','fontsize',20)
ylabel('Transducer length, \ity\rm (mm)','fontsize',20)
set(gca,'fontsize',20)
xticks([-5 -4 -3 -2 -1 0 1 2 3 4 5])
xticklabels({'-5','-4','-3','-2', '-1','0','1', '2','3','4', '5'})
title(' ')
c = colorbar;
ylabel(c,'Pressure, \itp_{pk-} \rm(MPa)','fontsize',20)
set(gca,'fontsize',20)
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
[cv,ch]=contourf(y,x,(rms_array),80);
set(ch,'edgecolor','none');
axis equal
axis([-floor(n_col/2)*dy+y0 floor(n_col/2)*dy+y0 -floor(n_row/2)*dx+x0 floor(n_row/2)*dx+x0])
axisLabel = get(gca,'YTickLabel');
set(gca,'YTickLabel',axisLabel(end:-1:1));
xlabel('Transducer width, \itx\rm (mm)','fontsize',20)
ylabel('Transducer length, \ity\rm (mm)','fontsize',20)
set(gca,'fontsize',20)
xticks([-5 -4 -3 -2 -1 0 1 2 3 4 5])
xticklabels({'-5','-4','-3','-2', '-1','0','1', '2','3','4', '5'})
title(' ')
c = colorbar;
ylabel(c,'Pressure, \itp_{RMS} \rm(MPa)','fontsize',20)
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
y = (-floor(n_col/2)*dy:dy:floor(n_col/2)*dy) + y0;
x = (-floor(n_row/2)*dx:dx:floor(n_row/2)*dx) + x0;

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
axis([-floor(n_col/2)*dy+y0 floor(n_col/2)*dy+y0 -floor(n_row/2)*dx+x0 floor(n_row/2)*dx+x0])
axisLabel = get(gca,'YTickLabel');
set(gca,'YTickLabel',axisLabel(end:-1:1));
xlabel('Transducer width, \itx\rm (mm)','fontsize',20)
ylabel('Transducer length, \ity\rm (mm)','fontsize',20)
set(gca,'fontsize',20)
xticks([-5 -4 -3 -2 -1 0 1 2 3 4 5])
xticklabels({'-5','-4','-3','-2', '-1','0','1', '2','3','4', '5'})
yticks([-5 -4 -3 -2 -1 0 1 2 3 4 5])
yticklabels({'-5','-4','-3','-2', '-1','0','1', '2','3','4', '5'})
max_p = max(max(max_array));
str = sprintf('P_{sp} = %.2f MPa', max_p)
title(str)
c = colorbar;
ylabel(c,'Pressure, \itp \rm(dB)','fontsize',20)
set(gca,'fontsize',20); 

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
fig8.Position = [520 51 908 747];
pos=get(gcf,'Position');
set(gcf, 'Position', [pos(1) pos(2) (width) (height)]);
plot(y,max_array_dB(max_line,:),'k','LineWidth',2)
str = sprintf('p_{sp} = %.2f MPa', max_p)
title(str)
xlabel('Transducer width, \itx\rm (mm)')
ylabel('Pressure, \itp_{pk} \rm(dB)')
set(gca,'fontsize',20)
frame = getframe(fig8)
im = frame2im(frame);
[A,map] = rgb2ind(im,256);
imwrite(A,map,'X.tif','tiff');
savefig(fig8,'X.fig')

fig9 = figure(9)
fig9.Position = [520 51 908 747];
pos=get(gcf,'Position');
set(gcf, 'Position', [pos(1) pos(2) (width) (height)]);
plot(x,max_array_dB(:,max_col),'k','LineWidth',2)
str = sprintf('p_{spp} = %.2f MPa', max_p)
title(str)
xlabel('Transducer length, \ity\rm (mm)')
ylabel('Pressure, \itp_{spp} \rm(dB)')
set(gca,'fontsize',20)
frame = getframe(fig9)
im = frame2im(frame);
[A,map] = rgb2ind(im,256);
imwrite(A,map,'Y.tif','tiff');
savefig(fig9,'Y.fig')
ylim([-30 0])
yticks([-30 -25 -20 -15 -10 -5 0])

fig7 = figure(7)
plot(x,max_array_dB(:,max_col))
fig7.Position = [520 51 908 747];
pos=get(gcf,'Position');
set(gcf, 'Position', [pos(1) pos(2) (width) (height)]);
str = sprintf('p_{sp} = %.2f MPa', max_p)
title(str)
xlabel('Transducer length, \ity \rm(mm)')
ylabel('Pressure, \itp_{pk} \rm(dB)')
frame = getframe(fig7)
im = frame2im(frame);
[A,map] = rgb2ind(im,256);
imwrite(A,map,'Y.tif','tiff');
savefig(fig7,'Y.fig')


%% Contour Plot

% max_array_norm = max_array/(max(max(max_array)));
% max_array_dB = 20*log10(max_array_norm);
% 
% % axes
% y = (-floor(n_row/2)*dy:dy:floor(n_row/2)*dy) + y0;
% x = (-floor(n_row/2)*dx:dx:floor(n_row/2)*dx) + x0;
% 
% width = 908;     % Width in pixels
% height = 747;    % Height in pixels
% alw = 0.75;    % AxesLineWidth
% fsz = 24;      % Fontsize
% lw = 1.5;      % LineWidth
% msz = 8;       % MarkerSize
% fig1 = figure(13);
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
xlabel('Transducer width, \itx\rm (mm)','fontsize',20)
ylabel('Transducer length, \ity\rm (mm)','fontsize',20)
set(gca,'fontsize',20)
xticks([-5 -4 -3 -2 -1 0 1 2 3 4 5])
xticklabels({'-5','-4','-3','-2', '-1','0','1', '2','3','4', '5'})
yticks([-5 -4 -3 -2 -1 0 1 2 3 4 5])
yticklabels({'-5','-4','-3','-2', '-1','0','1', '2','3','4', '5'})
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

%% point (0,0,0) is supposed to be 50mm focus from the z or (0,0,50)
%
% For the scans located in 
% C:\Users\cbawiec\Documents\INSERM\Project MUTATION\MUTATION_Prostate\Hydrophone Measurements\SondeCMUT_3D\March2017\42mm
% The phasing of the active elements was set to be 42 mm in depth
%

cd('D:\MUTATION\Champs de Pression\TETE1\jeudi 21\focale 50\f50_3D')

close all

path = pwd;
%% set transducer and measurement parameters
inner_rad = 7; % [mm] inner most transducer ring radius
outer_rad = 28.75; % [mm] outer most transducer ring radius
speed_of_sound = 1480; % [m/s]

%%Set Scan Parameters
focal_depth = 50; % [mm] (0,0,0) location for the scan files

z = -10:1:9; % z scan start:step size:stop [mm]
actual_z = 40:1:60; % actual z distance from transducer face
x = -5:0.5:5; % x scan start:step size:stop [mm]
y = -5:0.5:5; % y scan start:step size:stop [mm]

[X,Y,Z] = meshgrid(x,y,actual_z);

%%preallocate large arrays
field3d = zeros(length(x),length(y),length(z));
scanMapMax = zeros(length(x),length(y),length(z));


%%input data from already organized file into field variable
for ii=1:length(z)
    slice2dfile = fullfile(path,num2str(z(ii)),'Planar_UMSmap.txt');
    field3d(:,:,ii) = dlmread(slice2dfile,'',5,1);
end

%% go through each file and find maximum within the desired window where all the signals have converged
for jj = 1:length(z)
    for mm = 1:length(y)
        for kk = 1:length(x)   
            pointScanFile = fullfile(path,num2str(z(jj)), ...
                ['PlanarY' num2str(mm-1,'%03d') 'X' num2str(kk-1,'%03d') '.txt']);  
            pointScan = dlmread(pointScanFile,'',81,0); %time, voltage data, uknown data
            pointScan(:,3) = []; % get rid of unknown data
            
            % isolate the  converging part of the signal based
            % on calculated time of flight 
%             window_start = sqrt(actual_z(jj)^2 + (outer_rad + sqrt(x(kk)^2 + y(mm)^2))^2) ...
%                             / speed_of_sound / 1000 + 1e-6
            window_start = pointScan(1,1);
            window_end = window_start + 2e-6; 
            
            %find indices of voltage array for desired window of values
%             [~,start_bin]=histc(window_start,pointScan(:,1));
%             [~,end_bin]=histc(window_end,pointScan(:,1));
%             
%             pointScanMax = max(pointScan(start_bin:end_bin,2)); %find max value of voltage data
% %             pointScanMax = max(pointScan(:,2)); %find max value of voltage data
            pointScanMax = max(pointScan(1:length(pointScan),2))
            scanMapMax(kk,mm,jj) = pointScanMax; %save max value in matrix
        end
    end
end


jj = 21;
z = -10:1:10;

for mm = 1:length(y)
        for kk = 1:length(x)   
            pointScanFile = fullfile(path,num2str(z(jj)), ...
                ['scanTestY' num2str(mm-1,'%03d') 'X' num2str(kk-1,'%03d') '.txt']);  
            pointScan = dlmread(pointScanFile,'',81,0); %time, voltage data, uknown data
            pointScan(:,3) = []; % get rid of unknown data
            
            % isolate the  converging part of the signal based
            % on calculated time of flight 
%             window_start = sqrt(actual_z(jj)^2 + (outer_rad + sqrt(x(kk)^2 + y(mm)^2))^2) ...
%                             / speed_of_sound / 1000 + 1e-6
            window_start = pointScan(1,1);
            window_end = window_start + 2e-6; 
            
            %find indices of voltage array for desired window of values
%             [~,start_bin]=histc(window_start,pointScan(:,1));
%             [~,end_bin]=histc(window_end,pointScan(:,1));
%             
%             pointScanMax = max(pointScan(start_bin:end_bin,2)); %find max value of voltage data
% %             pointScanMax = max(pointScan(:,2)); %find max value of voltage data
            pointScanMax = max(pointScan(1:length(pointScan),2))
            scanMapMax(kk,mm,jj) = pointScanMax; %save max value in matrix
        end
    end

save('3dScanMapMax-F32mm.mat','scanMapMax');
 
% figure(1)

image_num = 7;
% imagesc(X,Y,field3d(:,:,image_num))
% title(['Planar UMSmap - 3D Field of image slice at ' num2str(z(image_num)) 'mm'])
% c=colorbar;
% ylabel(c,'Unknown - Volts?')
% xlabel('X-Distance [mm]')
% ylabel('Y-Distance [mm]')
% zlabel('Z-Distance [mm]')
% axis image

% figure(2)
% surfaceValue = max(field3d(:))/2;
% p = patch(isosurface(X,Y,Z,field3d,surfaceValue));
% isonormals(X,Y,Z,field3d,p)
% p.FaceColor = 'red';
% p.EdgeColor = 'none';
% daspect([1,1,1])
% view(3); 
% axis ([-5 5 -5 5 -15 15])
% title(['Planar UMSmap - IsoSurface of Half Max Intensity: ' num2str(surfaceValue)])

% xlabel('X-Distance [mm]')
% ylabel('Y-Distance [mm]')
% zlabel('Z-Distance [mm]')
% camlight 
% lighting gouraud

figure(3)
scanMapSurfaceMax = max(scanMapMax(:));
p = patch(isosurface(X,Y,Z,scanMapMax,scanMapSurfaceMax/sqrt(2))); % /sqrt(2) for half intensity
isonormals(X,Y,Z,scanMapMax,p)
p.Parent.ZDir = 'reverse';
p.FaceColor = 'red';
p.EdgeColor = 'none';

q = patch(isosurface(X,Y,Z,scanMapMax,scanMapSurfaceMax/sqrt(4))); % /sqrt(4) for quarter intensity
isonormals(X,Y,Z,scanMapMax,q)
q.FaceAlpha = 0.25;
q.FaceColor = 'blue';
q.EdgeColor = 'none';

r = patch(isosurface(X,Y,Z,scanMapMax,scanMapSurfaceMax/sqrt(10))); % /sqrt(10) for 5x less intensity
isonormals(X,Y,Z,scanMapMax,r)
r.FaceAlpha = 0.125;
r.FaceColor = 'magenta';
r.EdgeColor = 'none';

s = patch(isosurface(X,Y,Z,scanMapMax,scanMapSurfaceMax/sqrt(20))); % /sqrt(20) for 10x less intensity
isonormals(X,Y,Z,scanMapMax,s)
s.FaceAlpha = 0.0625;
s.FaceColor = 'cyan';
s.EdgeColor = 'none';

daspect([1,1,1])
view(3); 
axis ([-5 5 -5 5 40 60]);
title(sprintf('p_{sp} = %.2f MPa', max(scanMapMax(:))*2.38))
xlabel('Width, \itx \rm(mm)')%'fontsize',20)
ylabel('Length, \ity \rm(mm)')%'fontsize',20)
zlabel('US propagation axis, \itz \rm(mm)')%'fontsize',20)
set(gca,'fontsize',20)
hcb=legend([p q r s],'-3 dB', '-6 dB', '-9 dB', '-12 dB')
camlight 
lighting gouraud

cb_handle = get(hcb,'Title');
titleString = 'Acoustic pressure, \itp \rm (dB)';
set(cb_handle ,'String',titleString);

figure(4)
imagesc(x,y,scanMapMax(:,:,image_num))
c=colorbar;
title(['scanMapMax - 3D Field of image slice at ' num2str(actual_z(image_num)) 'mm'])
ylabel(c,'Volts')
xlabel('X-Distance [mm]')
ylabel('Y-Distance [mm]')
zlabel('Z-Distance [mm]')
axis image
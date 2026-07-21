% run caiman_register_tiff_stack.py before this

%% Load registered image stacks
% Load AO/NoAO volumes and choose the motion-correction workflow.
clear; clc; close all;
pathfolder = 'K:\AO_FACED2026\20260327_RBead_M507_tdtomato\20260327_RBead_M507\M507_data1_3D\'; cd(pathfolder);
load([pathfolder 'FACED_pro.mat']);
str1 = 'sampleAO'; str2 = 'sysAO';

oldprocessing = 0;

if oldprocessing==1
    %% Register the AO stack
    % Apply the legacy within-stack motion correction.
    load([pathfolder 'Img_sampleAO_avg_noInterCorr.mat']);
    Img_AO = img_avg; clear img_avg

    options.overwrite = true;
    registerSeg = 0; registerSegSize = []; % check first file if used
    [datareg,yxShiftAll] = Motion_correction_Mod(cat(3,max(Img_AO,[],3),Img_AO),...
        registerSeg,registerSegSize,options);

    Img_AO = datareg(:,:,2:end);

    %% Register the NoAO stack
    % Apply the same legacy correction to the NoAO volume.
    load([pathfolder 'Img_sysAO_avg_noInterCorr.mat']);
    Img_NoAO_reg = img_avg; clear img_avg

    options.overwrite = true;
    registerSeg = 0; registerSegSize = []; % check first file if used
    [datareg,yxShiftAll] = Motion_correction_Mod(cat(3,max(Img_NoAO_reg,[],3),Img_NoAO_reg),...
        registerSeg,registerSegSize,options);

    Img_NoAO_reg = datareg(:,:,2:end);

else
    Img_AO_ori = tiffreadVolume([pathfolder 'MotionCorr_Caiman\Img_sampleAO_avg_noInterCorr.tif']);
    Img_NoAO_reg_ori = tiffreadVolume([pathfolder 'MotionCorr_Caiman\Img_sysAO_avg_noInterCorr.tif']);

    aostr = 'Img_sampleAO_avg_noInterCorr_caiman_reg.tif';
    noaostr = 'Img_sysAO_avg_noInterCorr_caiman_reg.tif';
    Img_AO = tiffreadVolume([pathfolder 'MotionCorr_Caiman\' aostr]);
    Img_NoAO_reg = tiffreadVolume([pathfolder 'MotionCorr_Caiman\' noaostr]);

    binN = 50;
    minval = min([min(Img_AO_ori(:)) min(Img_AO(:))]);
    maxval = max([max(Img_AO_ori(:)) max(Img_AO(:))]);
    figure; histogram(Img_AO_ori(:),linspace(minval,maxval,binN)); hold on; histogram(Img_AO(:),linspace(minval,maxval,binN),'FaceAlpha', 0.4, 'EdgeAlpha', 0.3);
    legend('AO, Ori','AO, Caiman'); set(gca,'Yscale','log')

    binN = 50;
    minval = min([min(Img_NoAO_reg_ori(:)) min(Img_NoAO_reg(:))]);
    maxval = max([max(Img_NoAO_reg_ori(:)) max(Img_NoAO_reg(:))]);
    figure; histogram(Img_NoAO_reg_ori(:),linspace(minval,maxval,binN)); hold on; histogram(Img_NoAO_reg(:),linspace(minval,maxval,binN),'FaceAlpha', 0.4, 'EdgeAlpha', 0.3);
    legend('NoAO, Ori','NoAO, Caiman'); set(gca,'Yscale','log')
end

%% Align NoAO to AO
% Estimate the bulk XY offset and apply it to every NoAO plane.
registerSeg = 0; registerSegSize = []; % check first file if used
[datareg,yxShiftAll] = Motion_correction_Mod(cat(3,max(Img_AO,[],3),max(Img_NoAO_reg,[],3)),...
    registerSeg,registerSegSize,options);

% apply image shift in XY to all slices
[Img_NoAO_reg] = fixyxshift(Img_NoAO_reg,repmat(yxShiftAll(2,:),size(Img_NoAO_reg,3),1));

options.overwrite = true;
saveastiff(Img_NoAO_reg,[noaostr(1:end-4) 'bulkReg2sample.tif'],options);

%% Crop and save aligned volumes
% Keep the shared spatial and depth range used for comparison.
d1crop = 1:size(Img_AO,1);
d2crop = 201:600;
zCrop_Extra = 1:55;

dep_ao = dep_ao(zCrop_Extra);
dep_noao = dep_noao(zCrop_Extra);
Img_AO = Img_AO(d1crop,d2crop,zCrop_Extra);
Img_NoAO_reg = Img_NoAO_reg(d1crop,d2crop,zCrop_Extra);

saveastiff(Img_AO,[aostr(1:end-4) '_USE.tif'],options);
saveastiff(Img_NoAO_reg,[noaostr(1:end-4) 'bulkReg2sample_USE.tif'],options);

%% image display
% Display matched maximum-intensity projections of both volumes.
NFACED = 100;
dX_FACED = FOV_x/NFACED; %um
dZ_FACED = abs(mean(diff(dep_ao))); % um

x_dim = dX_FACED*size(Img_AO,2);
y_dim = FOV_y*size(Img_AO,1);

displayLim = [0 60000];
close all;
figure(100); imagesc(max(Img_AO,[],3)); colormap gray
pbaspect([x_dim y_dim 1]);
colorbar; clim(displayLim); axis off; title(['Max IP: ' str1]);
Colormap_create(); colormap(greenHotMap);
set(gcf,'Color',[1,1,1]); set(gca,'FontSize',16);
scalebar_mod2(x_dim,y_dim,50,'\mum','w','w');
% set(gcf,'Position',[-958 42 958 954]);

figure(200); imagesc(max(Img_NoAO_reg,[],3)); colormap gray
pbaspect([x_dim y_dim 1]);
colorbar; clim(displayLim); axis off; title(['Max IP: ' str2]);
Colormap_create(); colormap(greenHotMap);
set(gcf,'Color',[1,1,1]); set(gca,'FontSize',16);
scalebar_mod2(x_dim,y_dim,50,'\mum','w','w');
% set(gcf,'Position',[-1918 42 958 954]);

% saveFigs_mod(pplotXY,1,'.fig','Img');
% saveFigs_mod(pplotXY,1,'.pdf','Img');

%% Select XY line-profile ROIs
% Launch interactive ROI selection on the aligned projections.
ROI_XY();

%% MIP, single plane, XY
% Compare AO and NoAO projections over a selected depth range.
% close all;
displayLim = [1000 60000];

% croprg1 = [1:500]; croprg2 = 1:400; cropSample = 33; cropSys = 33; %M507_data1_3D
croprg1 = [1:500]; croprg2 = 1:400; cropSample = 4:55; cropSys = cropSample-3; %M507_data1_3D
sampleZoom = max(Img_AO(croprg1,croprg2,cropSample),[],3);
sysZoom = max(Img_NoAO_reg(croprg1,croprg2,cropSys),[],3);

x_dim = dX_FACED*size(sampleZoom,2);
y_dim = FOV_y*size(sampleZoom,1);

figure(500); imagesc(sampleZoom); colormap gray
pbaspect([x_dim y_dim 1]);
colorbar; clim(displayLim); axis off; title(['Max IP: ' str1])
set(gcf,'Color',[1,1,1]); set(gca,'FontSize',16);
scalebar_mod2(x_dim,y_dim,50,'\mum','w','w',1);
get(gca,'Position'); set(gca,'Position',ans);
h = colorbar; h.Position = [0.85 0.1095 0.02 0.2]; % For zoom-in data
h.Ticks = displayLim;
Colormap_create(); colormap(greenHotMap);

figure(600); imagesc(sysZoom); colormap gray
pbaspect([x_dim y_dim 1]);
colorbar; clim(displayLim); axis off; title(['Max IP: ' str2])
set(gcf,'Color',[1,1,1]); set(gca,'FontSize',16);
scalebar_mod2(x_dim,y_dim,50,'\mum','w','w',1);
get(gca,'Position'); set(gca,'Position',ans);
h = colorbar; h.Position = [0.85 0.1095 0.02 0.2]; % For zoom-in data
clim(displayLim);
h.Ticks = displayLim;
Colormap_create(); colormap(greenHotMap);

% saveFigs_mod(pplotXY,1,'.fig',['MIP_' num2str(dep_ao(cropSample(1)),'%0.1f') '_' num2str(dep_ao(cropSample(end)),'%0.1f') 'um']);
% saveFigs_mod(pplotXY,1,'.pdf',['MIP_' num2str(dep_ao(cropSample(1)),'%0.1f') '_' num2str(dep_ao(cropSample(end)),'%0.1f') 'um']);

%% Zoom region, single plane
% Generate matched zoomed projections for the manuscript figure.
% close all;
displayLim = [2000 32000];

croprg1 = [225:350]; croprg2 = 220:292; cropSample = 1:55; cropSys = 1:55; % M507_data1_3D, zoom 1
% croprg1 = [215:350]; croprg2 = 207:292; cropSample = 33; cropSys = 33; % M507_data1_3D, zoom 1
sampleZoom = max(Img_AO(croprg1,croprg2,cropSample),[],3);
croprg2 = croprg2+1;
sysZoom = max(Img_NoAO_reg(croprg1,croprg2,cropSys),[],3);

x_dim = dX_FACED*size(sampleZoom,2);
y_dim = FOV_y*size(sampleZoom,1);

figure(500); imagesc(sampleZoom); colormap gray
pbaspect([x_dim y_dim 1]);
colorbar; clim(displayLim); axis off; title(['Max IP: ' str1])
set(gcf,'Color',[1,1,1]); set(gca,'FontSize',16);
scalebar_mod2(x_dim,y_dim,5,'\mum','w','w',1);
get(gca,'Position'); set(gca,'Position',ans);
Colormap_create(); colormap(greenHotMap);
h = colorbar; h.Position = [0.85 0.1095 0.02 0.2]; % For zoom-in data
h.Ticks = displayLim;

figure(600); imagesc(sysZoom); colormap gray
pbaspect([x_dim y_dim 1]);
colorbar; clim(displayLim); axis off; title(['Max IP: ' str2])
set(gcf,'Color',[1,1,1]); set(gca,'FontSize',16);
scalebar_mod2(x_dim,y_dim,5,'\mum','w','w',1);
get(gca,'Position'); set(gca,'Position',ans);
Colormap_create(); colormap(greenHotMap);
h = colorbar; h.Position = [0.85 0.1095 0.02 0.2]; % For zoom-in data
clim(displayLim);
h.Ticks = displayLim;

%% Overlay saved line ROIs
% Draw the selected AO and NoAO profile paths on the zoomed images.
% Zoom1
load([pplotXY 'ROI_3.mat']);
% open sampleAO zoom
hold on; plot(x_line_ao-croprg2(1)+1,y_line_ao-croprg1(1)+1)
% open sysAO zoom
hold on; plot(x_line_noao-croprg2(1)+1,y_line_noao-croprg1(1)+1)

% saveFigs_mod(pplotXY,1,'.fig',['Zoom1_MIP_' num2str(dep_ao(cropSample(1)),'%0.1f') '_' num2str(dep_ao(cropSample(end)),'%0.1f') 'um']);
% saveFigs_mod(pplotXY,1,'.pdf',['Zoom1_MIP_' num2str(dep_ao(cropSample(1)),'%0.1f') '_' num2str(dep_ao(cropSample(end)),'%0.1f') 'um']);

%% Mark the zoom region
% Outline the crop on the full-field image.
hold on; line([croprg2(1) croprg2(end) croprg2(end) croprg2(1) croprg2(1)],...
    [croprg1(1) croprg1(1) croprg1(end) croprg1(end) croprg1(1)])

%% format, Line profiles
% Apply the final line-profile axes layout.
set(gca,'Position',[0.1585 0.1559 0.7465 0.7575]);

%% format, colorscale adjust
% Match the colormap and color scale across image panels.
displayLim = [2000 28000];
Colormap_create(); colormap(greenHotMap);
get(gca,'Position'); set(gca,'Position',ans);
colorbar off;
h = colorbar; h.Position = [0.82 0.1095 0.02 0.2]; % For zoom-in data
clim(displayLim);
h.Ticks = displayLim;

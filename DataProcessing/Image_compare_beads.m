% Compare fluorescent-bead stacks acquired with system AO and without AO.
clear; clc; close all;
pathfolder = 'K:\AO_FACED2026\20260416_RBead_vessel\FACED\2umRbead\'; cd(pathfolder);
psave = 'K:\AO_FACED2026\20260416_RBead_vessel\FACED\2umRbead\';
namestr = '2umRBead';
Xrg = [5 28-5]; %FACED-AO, 1035 nm, 03/2026

loaddata = 0;
% Identify the system-AO and NoAO acquisition folders from their name suffixes.
folderPro = dir([pathfolder namestr '*']);
AA = [folderPro.isdir];
folderPro(AA==0)=[]; clear AA;

for iii = 1:length(folderPro)
    if contains(folderPro(iii).name,'_sysAO')
        AOpath = [folderPro(iii).folder '\' folderPro(iii).name '\']
    elseif contains(folderPro(iii).name,'_NoAO')
        NoAOpath = [folderPro(iii).folder '\' folderPro(iii).name '\']
    end
end

framedig = 4;
reshapeData = 0;
Offset = nan;

% Load cached volumes when available; otherwise read raw FACED data and save it.
if loaddata==1
    Img_AO = tiffreadVolume([pathfolder 'Img_AO_avg.tif']);
    Img_NoAO = tiffreadVolume([pathfolder 'Img_NoAO_avg.tif']);
    load([pathfolder 'FACED_pro.mat']);
else
    [Img_AO,Nframe,dT,FOV_x,FOV_xNum,FOV_y,FOV_yNum,LineAvg,dep_ao] = read_SWZ(AOpath,reshapeData,Offset,Xrg,framedig);
    saveastiff(Img_AO,[pathfolder 'Img_AO_avg.tif']);
    [Img_NoAO,Nframe,dT,FOV_x,FOV_xNum,FOV_y,FOV_yNum,LineAvg,dep_noao] = read_SWZ(NoAOpath,reshapeData,Offset,Xrg,framedig);
    saveastiff(Img_NoAO,[pathfolder 'Img_NoAO_avg.tif']);
    save([pathfolder 'FACED_pro.mat'],'Offset','Xrg','Nframe','dT','FOV_x','FOV_xNum','FOV_y','FOV_yNum','LineAvg','dep_ao','dep_noao');
end
dep = round(dep_ao,2);
dep2 = round(dep_noao,2);
if unique(Xrg)==1
    NFACED = 128;
else
    NFACED = 100;
end
dX_FACED = FOV_x/100; %um

% Register the NoAO volume to the AO maximum-intensity projection using a global XY shift.
%% image shift
registerSeg = 0; registerSegSize = []; % check first file if used
options.maxLoop = 100; % maximum iterations allowed to register images;
options.previousRefFlag = 0;
options.extRefImage = max(Img_AO,[],3);
options.extRef = 1;
[datareg,yxShiftAll] = Motion_correction_Mod(cat(3,max(Img_AO,[],3),max(Img_NoAO,[],3)),...
        registerSeg,registerSegSize,options);
save([pathfolder 'FACED_pro.mat'],'registerSeg','registerSegSize','options','yxShiftAll','-append')

% apply image shift in XY to all slices
Img_NoAO_reg = zeros(size(Img_NoAO));
for iii = 1:size(Img_NoAO,3)
    Img_NoAO_reg(:,:,iii) = circshift(Img_NoAO(:,:,iii),yxShiftAll(2,:));
end
figure; imagesc_auto(max(Img_AO,[],3)); clim([0 5000]);
figure; imagesc_auto(max(Img_NoAO_reg,[],3)); clim([0 5000]);

% Display matched maximum-intensity projections using physical dimensions and common limits.
%% image display
display_lim = [0 8000];
close all;
x_dim = dX_FACED*size(Img_AO,2);
y_dim = FOV_y*size(Img_AO,1);
figure; imagesc(max(Img_AO,[],3)); colormap gray
pbaspect([x_dim y_dim 1]);
colorbar; clim(display_lim); axis off; title('Max IP: system AO');
set(gcf,'Color',[1,1,1]); set(gca,'FontSize',16);
scalebar_mod2(x_dim,y_dim,50,'\mum','w','w');
Colormap_create(); colormap(redHotMap);
get(gca,'Position'); set(gca,'Position',ans);
h = colorbar; h.Position = [0.85 0.1095 0.02 0.2]; % For large FOV data
h.Ticks = display_lim;

figure; imagesc(max(Img_NoAO_reg,[],3)); colormap gray
pbaspect([x_dim y_dim 1]);
colorbar; clim(display_lim); axis off; title('Max IP: no system AO')
set(gcf,'Color',[1,1,1]); set(gca,'FontSize',16);
scalebar_mod2(x_dim,y_dim,50,'\mum','w','w');
Colormap_create(); colormap(redHotMap);
get(gca,'Position'); set(gca,'Position',ans);
h = colorbar; h.Position = [0.85 0.1095 0.02 0.2]; % For large FOV data
h.Ticks = display_lim;

% Calculate a local peak-intensity ratio between the AO and registered NoAO volumes.
%% signal ratio
close all;
winsz = [3,1];
sig_ratio = zeros(size(Img_AO,1),size(Img_AO,2));
for iii = 1:(size(Img_AO,1)-winsz(1)+1)
    d1rg = iii:iii+winsz(1)-1;
    for jjj = 1:(size(Img_AO,2)-winsz(2)+1)
        d2rg = jjj:jjj+winsz(2)-1;
        sig_AO_curr = Img_AO(d1rg,d2rg,:);
        sig_NoAO_curr = Img_NoAO_reg(d1rg,d2rg,:);
        sig_ratio(iii,jjj) = max(sig_AO_curr(:))./max(sig_NoAO_curr(:));
    end
end

figure; imagesc(sig_ratio); colormap jet; colorbar;
pbaspect([FOV_xNum*FOV_x FOV_yNum*FOV_y 1]);
colorbar; clim([0 5]); axis off; title('I_{AO}/I_{NoAO}');
set(gcf,'Color',[1,1,1]); set(gca,'FontSize',16);
scalebar_mod2(x_dim,y_dim,20,'\mum','y','y');

% Segment beads from the mean AO/NoAO projection; roi_method selects the algorithm.
%% get ROIs
close all;
roi_method = 4;
I1 = max(Img_AO,[],3);
I2 = max(Img_NoAO_reg,[],3);
I = (I1+I2)/2;

if roi_method==1
    blursize = [2,2];
    Ic = medfilt2(I,blursize);
    exethre = mean(Ic(:))*5.1;
    BW = imbinarize(Ic,exethre);

    se = strel('rectangle',[5,5]);
    J = imerode(BW,se);
    se = strel('rectangle',[2,2]);
    J = imdilate(J,se);
    [imgd,N] = bwlabel(J);
    
    close all;
    figure; imagesc(I);
    figure; imagesc(BW);
    figure; imagesc(J);
    figure; imagesc(imgd);

elseif roi_method==2 % using cellpose
    blursize = [3,3];
    Ic = medfilt2(I,blursize);
    AA = double(Ic);
    Lim1 = 0; Lim2 = 1000; % small FOV, Lim2 = 2000; large FOV, Lim2 = 1000;
    AA(AA<Lim1) = Lim1; AA(AA>Lim2) = Lim2;
    AA = (AA-min(AA(:)))/(max(AA(:))-min(AA(:)));
    imwrite(AA,[psave 'I_all.png']);
    %load to cellpose

    II = imread([psave 'I_all_seg_0_1000.png']); II = II(:,:,1);
    [imgd,N] = bwlabel(II);
    close all;
    figure; imagesc(I);
    figure; imagesc(imgd); colormap jet; clim([-10 N]);
    figure; imagesc(max(Img_NoAO_reg,[],3));

% Process FACED tiles independently to compensate the lateral illumination envelope.
elseif roi_method==3  % using cellpose, with additional processing
    pseg = [psave 'seg\']; mkdir(pseg);
    blursize = [3,3];
    Ic = medfilt2(I,blursize);
    Lim1 = 5; Lim2 = 210; % small FOV, Lim2 = ; large FOV, Lim2 = 100;

    for iii = 1:FOV_xNum
        AA = Ic(:,(iii-1)*NFACED+1:iii*NFACED);
        curvem = 1./sum(AA); curvem = curvem/max(curvem);
        AA = AA.*repmat(curvem,size(Ic,1),1);
        AA(AA<Lim1) = Lim1; AA(AA>Lim2) = Lim2;
        AA = (AA-min(AA(:)))/(max(AA(:))-min(AA(:)));
        imwrite(AA,[psave 'seg\I_' num2str(iii) '.png']);
    end

    % reading
    masks = dir([pseg '*masks.png']);
    masks = masks(~[masks.isdir]);
    [~,idx] = sort([masks.datenum]);
    masks = masks(idx);
    IIlarge = [];
    for iii = 1:FOV_xNum
        if length(masks)==FOV_xNum
            movefile([masks(iii).folder '\' masks(iii).name],[pseg 'I_' num2str(iii) '_seg.png']);
        end
        II = imread([pseg 'I_' num2str(iii) '_seg.png']); II = II(:,:,1);
        IIlarge = cat(2,IIlarge,II);
    end
    II = IIlarge;
    [imgd,N] = bwlabel(II);
    close all;
    figure; imagesc(I); clim([0 2000]);
    figure; imagesc(imgd); colormap jet; clim([-10 N]);
    figure; imagesc(max(Img_NoAO_reg,[],3)); clim([0 2000]);
    
% Normalize each FACED tile and segment beads directly with Cellpose on the GPU.
elseif roi_method==4 % cellpose package in Matlab
    blursize = [3,3];
    Ic = medfilt2(I,blursize);
    % Lim1 = 5; Lim2 = 210; averageCellDiameter = 5; %small FOV, Lim2 = ; large FOV, Lim2 = 100;
    Lim1 = 5; Lim2 = 50; averageCellDiameter = 10; % 750 & 850 mm lenses to the DM

    IIlarge = [];
    for iii = 1:ceil(FOV_xNum)
        iii
        AA = Ic(:,(iii-1)*NFACED+1:min(iii*NFACED,size(Img_AO,2)));
        curvem = 1./sum(AA); curvem = curvem/max(curvem);
        AA = AA.*repmat(curvem,size(Ic,1),1);
        AA(AA<Lim1) = Lim1; AA(AA>Lim2) = Lim2;
        AA = (AA-min(AA(:)))/(max(AA(:))-min(AA(:)));

        Img = AA;
        cp = cellpose(Model="nuclei",ExecutionEnvironment="gpu");
        labels = segmentCells2D(cp,Img,ImageCellDiameter=averageCellDiameter);
        % B = labeloverlay(Img,labels); figure; imshow(B);
        IIlarge = cat(2,IIlarge,labels);
    end
    [imgd,N] = bwlabel(IIlarge);
    close all;
    figure; imagesc(I); clim([0 2000]);
    figure; imagesc(imgd); colormap jet; clim([-10 N]);
    figure; imagesc(max(Img_NoAO_reg,[],3)); clim([0 2000]);
end

save([psave 'para.mat'],'roi_method','blursize','Lim1','Lim2','imgd','N','IIlarge');

% Extract an axial intensity profile for every bead under both optical conditions.
%% get Z plots from ROIs
% For cal_method = 4, max or sum give similar results

cal_method = 4; %1-mean; 2-max; 3-top N brightness pixels; 4-daw a line
fig_display = 1;
removebkg = 0;
loadprofile = 0;
getwidth = 1;

if removebkg==1
    pplot = [psave 'ProfileNew_' num2str(cal_method) '_removebkg\']; mkdir(pplot);
else
    pplot = [psave 'ProfileNew_' num2str(cal_method) '\']; mkdir(pplot);
end

if loadprofile==1
    load([pplot 'profile.mat']);
else
    roimask = zeros(size(I,1),size(I,2),N);
    profile_all = zeros(size(Img_AO,3),N,2);
    maskcenter_all = zeros(2,N);
end

if getwidth==1
    w_all = nan(N,2);
end

% Build each bead mask, record its centroid, and measure its AO and NoAO profiles.
for roiind = 1:N
    roiind
    if loadprofile==1
        [];
    else
        roi_current = imgd;
        roi_current(roi_current~=roiind) = nan;
        roi_current(roi_current==roiind) = 1;
        roimask(:,:,roiind) = roi_current;
        temp = find(roi_current==1);
        [rr,cc] = ind2sub([size(I,1),size(I,2)],temp);
        maskcenter_all(:,roiind) = [nanmean(rr);nanmean(cc)];

        sig_AO_curr = Img_AO.*repmat(roi_current,[1,1,size(Img_AO,3)]);
        sig_NoAO_curr = Img_NoAO_reg.*repmat(roi_current,[1,1,size(Img_AO,3)]);

        if cal_method==1
            Nroi = nansum(nansum(roi_current));
            profile_all(:,roiind,1) = squeeze(nansum(nansum(sig_AO_curr,1),2))/Nroi;
            profile_all(:,roiind,2) = squeeze(nansum(nansum(sig_NoAO_curr,1),2))/Nroi;
        elseif cal_method==2
            profile_all(:,roiind,1) = squeeze(max(max(sig_AO_curr,[],1),[],2));
            profile_all(:,roiind,2) = squeeze(max(max(sig_NoAO_curr,[],1),[],2));
        elseif cal_method==3
            Nroi = nansum(nansum(roi_current));
            ratiouse = 0.15;
            for iii = 1:size(sig_AO_curr,3)
                temp = sig_AO_curr(:,:,iii);
                profile_all(iii,roiind,1) = (nansum(maxk(temp(:),round(Nroi*ratiouse))))/round(Nroi*ratiouse);
                temp = sig_NoAO_curr(:,:,iii);
                profile_all(iii,roiind,2) = (nansum(maxk(temp(:),round(Nroi*ratiouse))))/round(Nroi*ratiouse);
            end
        elseif cal_method==4
            temp = max(sig_AO_curr,[],3);
            [~,loc] = max(temp(:));
            [rr,cc] = ind2sub([size(I,1),size(I,2)],loc);
            profile_all(:,roiind,1) = sig_AO_curr(rr,cc,:);
            temp = max(sig_NoAO_curr,[],3);
            [~,loc] = max(temp(:));
            [rr,cc] = ind2sub([size(I,1),size(I,2)],loc);
            profile_all(:,roiind,2) = sig_NoAO_curr(rr,cc,:);
        end

        if removebkg==1
            profile_all(:,roiind,:) = profile_all(:,roiind,:)-repmat(min(profile_all(:,roiind,:)),[size(Img_AO,3),1,1]);
        end
    end

    % Reject boundary-truncated profiles, then calculate axial FWHM in micrometers.
    if getwidth==1
        maxval = squeeze(max(profile_all(:,roiind,:)));
        edgeval = [squeeze(profile_all(1,roiind,:)),squeeze(profile_all(end,roiind,:))];
        if any((edgeval(1,:)/maxval(1,1))>0.5)
            [];
        else
            w_all(roiind,1) = calc_width(profile_all(:,roiind,1),1,0.5)*abs(mean(diff(dep)));
            if w_all(roiind,1)<2
                w_all(roiind,1) = nan;
            end
        end
        if any((edgeval(2,:)/maxval(2,1))>0.5)
            [];
        else
            w_all(roiind,2) = calc_width(profile_all(:,roiind,2),1,0.5)*abs(mean(diff(dep)));
            if w_all(roiind,2)<2
                w_all(roiind,2) = nan;
            end
        end
    end
    
    if fig_display==1
        if (rem(roiind,20)==0)||(any(w_all(roiind,:)<3))
            close all;
            figure; plot(dep,profile_all(:,roiind,1),'r','LineWidth',1.5); hold on;
            plot(dep,profile_all(:,roiind,2),'b','LineWidth',1.5);
            legend(['AO ', num2str(w_all(roiind,1),'%.1f') ' \mum'],['No AO ', num2str(w_all(roiind,2),'%.1f') ' \mum'],'location','best');
            xlabel('Depth (\mum)');
            if cal_method==1
                ylabel('Average signal (A.U.)');
            elseif cal_method==2
                ylabel('Max signal (A.U.)');
            elseif ismember(cal_method,[3,4])
                ylabel('Signal (A.U.)');
            end
            title(['ROI ' num2str(roiind) ', Xavg = ' num2str(rem(maskcenter_all(2,roiind),NFACED),'%.0f') ', Yavg = ' num2str(maskcenter_all(1,roiind),'%.0f')]);
            set(gcf,'Color',[1,1,1]); xlim([min(dep) max(dep)]); set(gca,'FontSize',16);
            grid on; set(gca,'LineWidth',1.5);
            if mean(diff(dep))<0
                set(gca, 'XDir','reverse');
            end
            set(gca,'Position',[0.1585    0.1559    0.7465    0.7575]);
            saveFigs_mod(pplot,1,'.fig',['ROI_' num2str(roiind)]);
            saveFigs_mod(pplot,1,'.tif',['ROI_' num2str(roiind)]);
        end
    end
end
save([pplot 'profile_.mat'],'profile_all','maskcenter_all','removebkg','roimask','imgd','w_all','-v7.3');

%%
close all;
plotLim = 5;
figure; set(gca,'Position',[0.2   0.1100    0.6693    0.8150]);
scatter((maskcenter_all(2,:)-1)*dX_FACED,(maskcenter_all(1,:)-1)*FOV_y,25,w_all(:,1),'filled');
pbaspect([x_dim y_dim 1]);
xlabel('X location (\mum)'); ylabel('Y location (\mum)'); 
title(['\deltaz_{sysAO}']);
colormap jet; colorbar; clim([0 plotLim]);
set(gcf,'Color',[1,1,1]); xlim([0 (size(I,2)-1)*dX_FACED]); set(gca,'FontSize',16);
ylim([0 (size(I,1)-1)*FOV_y]); set(gca, 'YDir','reverse');
set(gca,'LineWidth',1.5); box on;
h = colorbar; h.Position = [0.90 0.1095 0.02 0.2];
h.Ticks = [0 plotLim]; h.Title.String = '\mum';

edges = 2:0.2:max(w_all(:));
edges = [edges edges(end)+0.2];
figure; set(gca,'Position',[0.2    0.15    0.6693    0.75]);
histogram(w_all(:,1),'BinEdges', edges);
xlabel('FWHM (\mum)'); ylabel('#'); 
title(['\deltaz_{sysAO}, mean = ' num2str(nanmean(w_all(:,1)),'%0.2f') ' , median = ' num2str(nanmedian(w_all(:,1)),'%0.2f')]);
set(gcf,'Color',[1,1,1]); set(gca,'FontSize',16);
set(gca,'LineWidth',1.5); box on;

figure; set(gca,'Position',[0.2    0.1100    0.6693    0.8150]);
scatter((maskcenter_all(2,:)-1)*dX_FACED,(maskcenter_all(1,:)-1)*FOV_y,25,w_all(:,2),'filled');
pbaspect([x_dim y_dim 1]);
xlabel('X location (\mum)'); ylabel('Y location (\mum)'); 
title('\deltaz_{NoAO}');
colormap jet; colorbar; clim([0 plotLim]);
set(gcf,'Color',[1,1,1]); xlim([0 (size(I,2)-1)*dX_FACED]); set(gca,'FontSize',16);
ylim([0 (size(I,1)-1)*FOV_y]); set(gca, 'YDir','reverse');
set(gca,'LineWidth',1.5); box on;
h = colorbar; h.Position = [0.90 0.1095 0.02 0.2];
h.Ticks = [0 plotLim]; h.Title.String = '\mum';

figure; set(gca,'Position',[0.2    0.15    0.6693    0.75]);
histogram(w_all(:,2),'BinEdges', edges);
xlabel('FWHM (\mum)'); ylabel('#'); 
title(['\deltaz_{NoAO}, mean = ' num2str(nanmean(w_all(:,2)),'%0.2f') ' , median = ' num2str(nanmedian(w_all(:,2)),'%0.2f')]);
set(gcf,'Color',[1,1,1]); set(gca,'FontSize',16);
set(gca,'LineWidth',1.5); box on;

% Display matched orthogonal planes through a representative, dataset-specific bead region.
%% Display SliceXZ image
close all;
displayLim = [0 10000];

% Bead 439
% figure; imagesc(roimask(:,:,439)); 
% figure; imagesc_auto(max(Img_NoAO_reg,[],3)); clim(displayLim);
% figure; imagesc_auto(max(Img_AO,[],3)); clim(displayLim);
% linkaxes_cus();

d1use = 181:192; d1_z = FOV_y;
d2use = 370:395; d2_z = dX_FACED;
d3_z = abs(mean(diff(dep)));

% XZ
figure; imagesc(squeeze(Img_NoAO_reg(187,d2use,:)).'); clim(displayLim);
x_dim = d2_z*length(d2use); y_dim = d3_z*size(Img_NoAO_reg,3); Beads_plot_format(); title(['XZ, NoAO']);
figure; imagesc(squeeze(Img_AO(188,d2use,:)).');  clim(displayLim); Beads_plot_format(); title(['XZ, sysAO']);

% YZ
figure; imagesc(squeeze(Img_NoAO_reg(d1use,374,:)).');  clim(displayLim);
x_dim = d1_z*length(d1use); y_dim = d3_z*size(Img_NoAO_reg,3); Beads_plot_format(); title(['YZ, NoAO']);
figure; imagesc(squeeze(Img_AO(d1use,375,:)).');  clim(displayLim); Beads_plot_format(); title(['YZ, sysAO']);

% XY
figure; imagesc(squeeze(Img_NoAO_reg(d1use,d2use,35)));  clim(displayLim);
x_dim = d2_z*length(d2use); y_dim = d1_z*length(d1use); Beads_plot_format(); title(['XY, NoAO']);
figure; imagesc(squeeze(Img_AO(d1use,d2use,36)));  clim(displayLim); Beads_plot_format(); title(['XY, sysAO']);

saveFigs_mod(pplot,0,'.fig','Beads_SinglePlane');
saveFigs_mod(pplot,0,'.pdf','Beads_SinglePlane');

% Optionally reject beads by profile-edge signal, ROI area, or peak-depth disagreement.
%% Execlude out of focus beads
excludeBeads = 0;
if excludeBeads==0
    [];
else
    if excludeBeads==1
        exethre = 0.6;
        roi_select = ones(1,N);
        % Build each bead mask, record its centroid, and measure its AO and NoAO profiles.
for roiind = 1:N
            temp = squeeze(profile_all(:,roiind,:));
            maxval = max(temp);
            if (temp(1,1)/maxval(1,1)>exethre)||(temp(end,1)/maxval(1,1)>exethre)||(temp(1,2)/maxval(1,2)>exethre)||(temp(end,2)/maxval(1,2)>exethre)
                roi_select(roiind) = 0;
            end          
        end
    elseif excludeBeads==2
        roi_select = squeeze(nansum(nansum(roimask,1),2));
        exethre = 20;
        roi_select = roi_select>exethre; % larger than 20 pixels
    elseif excludeBeads==3
        [vao,lao] = max(profile_all(:,:,1));
        [vnoao,lnoao] = max(profile_all(:,:,2));
        roi_select = ones(1,N);
        roi_select(abs(lao-lnoao)>25) = 0;
    end

    profile_all(:,find(roi_select==0),:) = [];
    maskcenter_all(:,find(roi_select==0)) = [];
    roimask(:,:,find(roi_select==0)) = [];
    save([pplot 'profile_execlude_' num2str(excludeBeads) '.mat'],'profile_all','maskcenter_all','roi_select','exethre','-v7.3');
end

% Inspect the axial profiles, mean images, and mask for one selected bead.
%% plot for a specific bead
[vao,lao] = max(profile_all(:,:,1));
[vnoao,lnoao] = max(profile_all(:,:,2));

roiind = 200;
figure; plot(profile_all(:,roiind,1),'r'); hold on; plot(profile_all(:,roiind,2),'b');

sig_AO_curr = Img_AO.*repmat(roimask(:,:,roiind),[1,1,size(Img_AO,3)]);
sig_NoAO_curr = Img_NoAO_reg.*repmat(roimask(:,:,roiind),[1,1,size(Img_AO,3)]);
figure; imagesc_auto(mean(sig_AO_curr,3));
figure; imagesc_auto(mean(sig_NoAO_curr,3));
figure; imagesc(roimask(:,:,roiind));

% Relate AO/NoAO signal gain to global position, within-tile position, and ROI size.
%% Location based analysis
close all;
sig_ratio = max(profile_all(:,:,1))./max(profile_all(:,:,2));
figure; plot((maskcenter_all(2,:)-1)*dX_FACED,sig_ratio,'k*','MarkerSize',8);
xlabel('X location (\mum)'); ylabel('I_{AO}/I_{NoAO}');
set(gcf,'Color',[1,1,1]); xlim([0 (size(I,2)-1)*dX_FACED]); set(gca,'FontSize',16);
grid on; set(gca,'LineWidth',1.5);

figure; plot((maskcenter_all(1,:)-1)*FOV_y,sig_ratio,'k*','MarkerSize',8);
xlabel('Y location (\mum)'); ylabel('I_{AO}/I_{NoAO}');
set(gcf,'Color',[1,1,1]); xlim([0 (size(I,1)-1)*FOV_y]); set(gca,'FontSize',16);
grid on; set(gca,'LineWidth',1.5);

figure; plot((rem(maskcenter_all(2,:),NFACED)-1)*dX_FACED,sig_ratio,'k*','MarkerSize',8);
xlabel('X location (\mum)'); ylabel('I_{AO}/I_{NoAO}');
set(gcf,'Color',[1,1,1]); xlim([0 NFACED*dX_FACED]); set(gca,'FontSize',16);
grid on; set(gca,'LineWidth',1.5);

figure; scatter((maskcenter_all(2,:)-1)*dX_FACED,(maskcenter_all(1,:)-1)*FOV_y,[],sig_ratio,'filled');
xlabel('X location (\mum)'); ylabel('Y location (\mum)'); title('I_{AO}/I_{NoAO}');
colormap jet; colorbar; clim([0 3])
set(gcf,'Color',[1,1,1]); xlim([0 (size(I,2)-1)*dX_FACED]); set(gca,'FontSize',16);
ylim([0 (size(I,1)-1)*FOV_y]); set(gca, 'YDir','reverse');
grid on; set(gca,'LineWidth',1.5); box on;

figure; scatter(maskcenter_all(2,:),maskcenter_all(1,:),[],sig_ratio,'filled');
xlabel('X location (#)'); ylabel('Y location (#)'); title('I_{AO}/I_{NoAO}');
colormap jet; colorbar; clim([0 3])
set(gcf,'Color',[1,1,1]); xlim([0 size(I,2)]); set(gca,'FontSize',16);
ylim([0 size(I,1)]); set(gca, 'YDir','reverse');
grid on; set(gca,'LineWidth',1.5); box on;

figure; imagesc(imgd);
xlabel('X location (#)'); ylabel('Y location (#)'); title('ROI');
colormap jet; colorbar; clim([-10 N])
set(gcf,'Color',[1,1,1]); xlim([0 size(I,2)]); set(gca,'FontSize',16);
ylim([0 size(I,1)]); set(gca, 'YDir','reverse');
set(gca,'LineWidth',1.5); box on;

figure; scatter((rem(maskcenter_all(2,:),NFACED)-1)*dX_FACED,(maskcenter_all(1,:)-1)*FOV_y,[],sig_ratio,'filled');
xlabel('X location (\mum)'); ylabel('Y location (\mum)'); title('I_{AO}/I_{NoAO}');
colormap jet; colorbar; clim([0 3])
set(gcf,'Color',[1,1,1]); xlim([0 NFACED*dX_FACED]); set(gca,'FontSize',16);
ylim([0 (size(I,1)-1)*FOV_y]); set(gca, 'YDir','reverse');
grid on; set(gca,'LineWidth',1.5); box on;

sz = squeeze(nansum(nansum(roimask,1),2));
[~,loci] = sort(sz);
figure; plot(sz(loci),sig_ratio(loci),'k*');
xlabel('ROI size (#)'); ylabel('I_{AO}/I_{NoAO}');
set(gcf,'Color',[1,1,1]); set(gca,'FontSize',16);
grid on; set(gca,'LineWidth',1.5);

%%
close all;
figure; imagesc(([1:size(I,2)]-1)*dX_FACED,([1:size(I,1)]-1)*FOV_y,I); colormap gray; clim([0 3000]);
xlabel('X location (\mum)'); ylabel('Y location (\mum)'); title('Beads signal (A.U.)'); colorbar;
set(gcf,'Color',[1,1,1]); xlim([0 (size(I,2)-1)*dX_FACED]); set(gca,'FontSize',16);
ylim([0 (size(I,1)-1)*FOV_y]);
set(gca,'LineWidth',1.5);
figure; imagesc(([1:size(I,2)]-1)*dX_FACED,([1:size(I,1)]-1)*FOV_y,imgd); colormap jet; clim([-10 N]);
xlabel('X location (\mum)'); ylabel('Y location (\mum)'); title('ROI index'); colorbar;
set(gcf,'Color',[1,1,1]); xlim([0 (size(I,2)-1)*dX_FACED]); set(gca,'FontSize',16);
ylim([0 (size(I,1)-1)*FOV_y]);
set(gca,'LineWidth',1.5);
% saveFigs_mod(pplot,0,'.fig','Profile_All');

%%
imagesc(Img_AO(:,:,9)); colormap gray
pbaspect([FOV_x/100*NFACED*FOV_xNum FOV_yNum*FOV_y 1]);
% colorbar; clim([50 3500]); axis off; title('Max IP: system AO')
% set(gcf,'Color',[1,1,1]); set(gca,'FontSize',16);

imagesc(Img_NoAO(:,:,15)); colormap gray
pbaspect([FOV_x/100*NFACED*FOV_xNum FOV_yNum*FOV_y 1]);
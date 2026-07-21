%% Load data
% Load registered AO/NoAO image stacks and calculate mean images.
clear; clc; close all;
foldernum = 9;
pathfolder = ['K:\AO_FACED2026\20260416_RBead_vessel\FACED\data' num2str(foldernum) '_VS\']; cd(pathfolder);
aofileInd = 0; noaofileInd = 1;
aopath = [pathfolder 'data' num2str(foldernum) '_sampleAO\Olympus25X_data' num2str(foldernum) '_sampleAO_TS_CH0_CHA_F' num2str(aofileInd) '_Pkg\'];
aostr = [aopath 'F' num2str(aofileInd) '_RM_ref1.tif'];
noaopath = [pathfolder 'data' num2str(foldernum) '_sysAO\Olympus25X_data' num2str(foldernum) '_sysAO_TS_CH0_CHA_F' num2str(noaofileInd) '_Pkg\'];
noaostr = [noaopath 'F' num2str(noaofileInd) '_RM_ref1.tif'];

pplotT = [pathfolder 'ProfileT_bicubic1_202604_compare_AOF' num2str(aofileInd) '_NoAOF' num2str(noaofileInd) '\']; mkdir(pplotT);

Img_AO = tiffreadVolume(aostr);
Img_NoAO_reg = tiffreadVolume(noaostr);
str1 = 'sampleAO';
str2 = 'sysAO';

save([pathfolder 'pro.mat'],'aopath','noaopath','aofileInd','noaofileInd');

IAO = mean(Img_AO,3);
INoAO = mean(Img_NoAO_reg,3);

%% Display mean images
% Show the paired fields of view for ROI selection.
load([pathfolder 'FACED_pro.mat']);
NFACED = 100;
dX_FACED = FOV_x/100; %um

x_dim = dX_FACED*size(IAO,2);
y_dim = FOV_y*size(IAO,1);

displayLim = [0 4500]; close all;
figure; imagesc(IAO); colormap gray;
pbaspect([x_dim y_dim 1]);
colorbar; clim(displayLim); axis off; title(['Mean IP: ' str1])
set(gcf,'Color',[1,1,1]); set(gca,'FontSize',16);
scalebar_mod2(x_dim,y_dim,50,'\mum','w','w');
get(gca,'Position'); set(gca,'Position',ans);
h = colorbar; h.Position = [0.85 0.1095 0.02 0.2]; % For large FOV data
h.Ticks = displayLim;
Colormap_create(); colormap(redHotMap);

figure; imagesc(INoAO); colormap gray;
pbaspect([x_dim y_dim 1]);
colorbar; clim(displayLim); axis off; title(['Mean IP: ' str2])
set(gcf,'Color',[1,1,1]); set(gca,'FontSize',16);
scalebar_mod2(x_dim,y_dim,50,'\mum','w','w');
get(gca,'Position'); set(gca,'Position',ans);
h = colorbar; h.Position = [0.85 0.1095 0.02 0.2]; % For large FOV data
h.Ticks = displayLim;
Colormap_create(); colormap(redHotMap);

%% ROI selection on the whole FOV
% Load vessel ROIs and generate their AO/NoAO kymographs.
[sROI_AO] = ReadImageJROI([aopath,'SampleAO_ROI_seg_F' num2str(aofileInd) '.zip']);
[sROI_NoAO] = ReadImageJROI([noaopath,'SysAO_ROI_seg_F' num2str(noaofileInd) '.zip']);

colorall = 'yrbcmgk';
for iii = 1:length(sROI_AO)
    close all;
    clear ptsi total_lengthAO XZSlice_AO XZSlice_NoAO

    roiAO = sROI_AO{1,iii};
    roiNoAO = sROI_NoAO{1,iii};
    isOvalAO = localIsOvalROI(roiAO);
    isOvalNoAO = localIsOvalROI(roiNoAO);
    isOvalROI = isOvalAO || isOvalNoAO;

    if xor(isOvalAO, isOvalNoAO)
        warning('ROI %d:  AO ROI and NoAO ROI types do not match. Treating this pair as oval ROIs.', iii);
    end

    Lim1 = displayLim(1);
    Lim2 = displayLim(2);

    if isOvalROI
        maskAO = localBuildOvalMask(roiAO, size(Img_AO,1), size(Img_AO,2));
        maskNoAO = localBuildOvalMask(roiNoAO, size(Img_NoAO_reg,1), size(Img_NoAO_reg,2));

        ROI_AO_line = localMeanTraceFromMask(double(Img_AO), maskAO);
        ROI_NoAO_line = localMeanTraceFromMask(double(Img_NoAO_reg), maskNoAO);

        contrast_cal = 2; Nratio = 0.15;
        [RMS_C_AO,RMS_C_NoAO] = VesselLineCst_v2(ROI_AO_line,ROI_NoAO_line,contrast_cal,Nratio);
        tt_XZ = [0:(length(ROI_AO_line)-1)]*dT;

        f = figure(300);
        set(gcf,'Position',[30 250 1200 500]);

        subplot('Position',[0.05 0.18 0.2 0.64]);
        imagesc(max(double(Img_NoAO_reg),[],3)); colormap gray; clim([Lim1 Lim2]); axis image off;
        title([str2 ' oval ROI']); hold on;
        localPlotMaskOutline(maskNoAO, colorall(1));

        subplot('Position',[0.3 0.18 0.2 0.64]);
        imagesc(max(double(Img_AO),[],3)); colormap gray; clim([Lim1 Lim2]); axis image off;
        title([str1 ' oval ROI']); hold on;
        localPlotMaskOutline(maskAO, colorall(1));

        subplot('Position',[0.58 0.25 0.37 0.5]);
        plot(tt_XZ,ROI_AO_line,'r','LineWidth',1.5); hold on;
        plot(tt_XZ,ROI_NoAO_line,'b','LineWidth',1.5);
        legend(str1,str2,'Location','NorthOutside','NumColumns',2);
        ylim([0 1.1*max(cat(1,ROI_AO_line(:),ROI_NoAO_line(:)))]);
        xlabel('t (s)'); ylabel('Signal (A.U.)');
        set(gcf,'Color',[1,1,1]); set(gca,'FontSize',16);
        grid off; set(gca,'LineWidth',1.5);

        resp_bht = menu('Select','Brightness only','Brightness & contrast');
        if resp_bht==1
            labelN = 1;
        elseif resp_bht==2
            labelN = 2;
        else
            labelN = 0;
        end
        flatmethod = [];
        up_factor = 1;
        d1_size = [];
        d2_size = [];
        XZ_AO = [];
        XZ_NoAO = [];
        XZSlice_AO = [];
        XZSlice_NoAO = [];
        s_final_AO = [];
        s_final_NoAO = [];
        ROI_AO = [];
        ROI_AO_up = [];
        ROI_AO_up_flat = [];
        ROI_NoAO = [];
        ROI_NoAO_up = [];
        ROI_NoAO_up_flat = [];
        dXY_FACED = [];
        total_lengthAO = [];
        win_paral = [];
        roiType = 'oval';

        save([pplotT 'ROI_' num2str(iii) '.mat'], ...
            'maskAO','maskNoAO','ROI_AO_line','ROI_NoAO_line','tt_XZ', ...
            'RMS_C_AO','RMS_C_NoAO','labelN','roiType');
        save([pplotT 'ROI_' num2str(iii) '_1.mat'], ...
            'XZ_NoAO','XZ_AO','s_final_AO','flatmethod','up_factor','ROI_AO','ROI_AO_up','ROI_AO_up_flat','ROI_AO_line', ...
            's_final_NoAO','ROI_NoAO','ROI_NoAO_up','ROI_NoAO_up_flat','ROI_NoAO_line', ...
            'tt_XZ','RMS_C_AO','RMS_C_NoAO','d1_size','d2_size','labelN', ...
            'maskAO','maskNoAO','roiType');
        saveas(gcf, [pplotT 'ROI_' num2str(iii) '_1_LP.fig']);
        saveas(gcf, [pplotT 'ROI_' num2str(iii) '.fig']);
        saveas(gcf, [pplotT 'ROI_' num2str(iii) '.tif']);
        close all;
        continue;
    end

    win_paral = 0; %linspace(-0.6,0.6,3); % shift, with the unit of um
    ptsi = roiAO.mnCoordinates+1;
    [total_lengthAO] = getSegmentLinelength(ptsi,dX_FACED,FOV_y); % physical distance
    [XZSlice_AO] = SliceXZ(double(Img_AO),ptsi(:,1),ptsi(:,2),round(total_lengthAO*10),'bicubic',win_paral,dX_FACED,FOV_y);

    ptsi = roiNoAO.mnCoordinates+1;
    [total_lengthNoAO] = getSegmentLinelength(ptsi,dX_FACED,FOV_y);
    if total_lengthNoAO~=total_lengthAO % should be the same!
        total_lengthNoAO = total_lengthAO;
    end
    [XZSlice_NoAO] = SliceXZ(double(Img_NoAO_reg),ptsi(:,1),ptsi(:,2),round(total_lengthNoAO*10),'bicubic',win_paral,dX_FACED,FOV_y);

    dXY_FACED = total_lengthAO/round(total_lengthAO*10);
    save([pplotT 'ROI_' num2str(iii) '.mat'],'XZSlice_AO','XZSlice_NoAO','total_lengthAO','dXY_FACED','win_paral');

    %% Select kymograph ROIs and extract line profiles
    % Select vessel bands, flatten them, and average across width.
    close all;
    f = figure; set(gcf,'unit','centimeter','Position',[2 2 total_lengthAO/5 24]);
    sp1 = subplot(1,2,1); set(sp1,'Position',[0.05 0.05 0.4 0.9]);
    imagesc(XZSlice_NoAO); colormap gray; clim([Lim1 Lim2]); hold on; title([str2]);
    sp2 = subplot(1,2,2); set(sp2,'Position',[0.55 0.05 0.4 0.9]);
    imagesc(XZSlice_AO); colormap gray; clim([Lim1 Lim2]); hold on; title([str1]);
    linkaxes([sp1,sp2]);

    %
    MM = 0;
    while(1)
        while(1)
            resp = menu('Select a RBC line in XZ (the dark stripe)?','Yes','No','Increase brightness','Decrease brightness');
            if resp==1
                % draw the line along the RBC traveling direction
                % pick similar lines in NoAO & AO
                set(f,'currentaxes',sp2); [x2,y2,c] = improfile; %#ok<NASGU>
                x2 = [x2(1) x2(end)].'; y2 = [y2(1) y2(end)].'; % two points on the line
                x2 = round(x2); y2 = round(y2);
                hold on; hh2 = plot(x2,y2,'y-','LineWidth',2);
                % set(f,'currentaxes',sp1); hh1 = images.roi.Polygon(gca,'Position',[x2,y2]);
                set(f,'currentaxes',sp1); hh1 = images.roi.Line(gca,'Position',[x2,y2]);
                pos = customWait(hh1); %#ok<NASGU>
                x1 = hh1.Position(:,1); y1 = hh1.Position(:,2);
                x1 = round(x1); y1 = round(y1);
                delete(hh1); hh1 = plot(x1,y1,'y-','LineWidth',2);
                break;
            elseif resp==2
                break;
            elseif resp==3
                Lim2 = Lim2/1.25;
                set(f,'currentaxes',sp2); clim([Lim1 Lim2]);
                set(f,'currentaxes',sp1); clim([Lim1 Lim2]);
            elseif resp==4
                Lim2 = Lim2*1.25;
                set(f,'currentaxes',sp2); clim([Lim1 Lim2]);
                set(f,'currentaxes',sp1); clim([Lim1 Lim2]);
            end
        end
        if resp==2
            break;
        end

        resp = menu('Process the current ROI?','Yes','No');
        if resp==1
            MM = MM+1;
            XZ_NoAO = [x1,y1]; XZ_AO = [x2,y2];
            hh1.Color = colorall(mod(MM-1,length(colorall))+1);
            hh2.Color = colorall(mod(MM-1,length(colorall))+1);
            set(f,'currentaxes',sp2);
            txth = text(mean(x2),mean(y2),num2str(MM),'FontSize',50,'Color',colorall(mod(MM-1,length(colorall))+1));

            % Image flattening
            s_final_AO = linspace(XZ_AO(1,2),XZ_AO(2,2),length(XZ_AO(1,1):XZ_AO(2,1)));
            s_final_NoAO = linspace(XZ_NoAO(1,2),XZ_NoAO(2,2),length(XZ_NoAO(1,1):XZ_NoAO(2,1)));
            flatmethod = 3; up_factor = 8;

            d1_size = round(0.8/2/dT); % half of the size, 201*dT
            d2_size = 50; % half of the size
            while(1)
                [ROI_AO,ROI_AO_up,ROI_AO_up_flat,ROI_AO_line] = imgflat(XZSlice_AO,XZ_AO,d1_size,d2_size,s_final_AO,flatmethod,up_factor);
                [ROI_NoAO,ROI_NoAO_up,ROI_NoAO_up_flat,ROI_NoAO_line] = imgflat(XZSlice_NoAO,XZ_NoAO,d1_size,d2_size,s_final_NoAO,flatmethod,up_factor);

                contrast_cal = 2;
                Nratio = 0.15;
                [RMS_C_AO,RMS_C_NoAO] = VesselLineCst_v2(ROI_AO_line,ROI_NoAO_line,contrast_cal,Nratio);

                slidingWindowSec = 0.1;
                adaptiveFwhmThreshold = 0.1;
                [PVD_AO_curr,~,PVDMethod_AO_curr] = localPeakValleyDistinguishability_v2( ...
                    ROI_AO_line,Nratio,'adaptive',dT/up_factor,slidingWindowSec,adaptiveFwhmThreshold);
                if PVDMethod_AO_curr=="sliding"
                    NoAOMethod_curr = 'slidingaverage';
                else
                    NoAOMethod_curr = char(PVDMethod_AO_curr);
                end
                [PVD_NoAO_curr,~,PVDMethod_NoAO_curr] = localPeakValleyDistinguishability_v2( ...
                    ROI_NoAO_line,Nratio,NoAOMethod_curr,dT/up_factor,slidingWindowSec,adaptiveFwhmThreshold);

                tt_XZ = [0:(length(ROI_AO_line)-1)]*dT/up_factor;
                close(figure(300));
                figure(300); set(gcf,'Position',[30 250 1200 500]);
                subplot('Position',[0.05 0.58 0.2 0.3]);
                imagesc(ROI_AO_up); clim([Lim1 Lim2]); axis off; title('ROI AO');
                set(gcf,'Color',[1,1,1]); set(gca,'FontSize',16);
                subplot('Position',[0.3 0.58 0.2 0.3]);
                imagesc(ROI_AO_up_flat); clim([Lim1 Lim2]); axis off; title('ROI AO flat');
                set(gcf,'Color',[1,1,1]); set(gca,'FontSize',16);
                subplot('Position',[0.05 0.1 0.2 0.3]);
                imagesc(ROI_NoAO_up); clim([Lim1 Lim2]); axis off; title('ROI NoAO');
                set(gcf,'Color',[1,1,1]); set(gca,'FontSize',16);
                subplot('Position',[0.3 0.1 0.2 0.3]);
                imagesc(ROI_NoAO_up_flat); clim([Lim1 Lim2]); axis off; title('ROI NoAO flat, 0.1 s');
                get(gca,'Position'); set(gca,'Position',ans);
                h = colorbar; h.Position = [0.52 0.1095 0.02 0.2]; % For large FOV data
                h.Ticks = [Lim1 Lim2];
                Colormap_create(); colormap(redHotMap);
                set(gcf,'Color',[1,1,1]); set(gca,'FontSize',16);
                scalebar_mod2(0,max(tt_XZ),0.1,'s','w','w');
                scalebar_mod2(size(ROI_NoAO_up_flat,2)*dXY_FACED,0,2,'\mum','w','w',2);

                subplot('Position',[0.63 0.4 0.35 0.4]);
                plot(tt_XZ,ROI_AO_line,'r','LineWidth',1.5); hold on;
                plot(tt_XZ,ROI_NoAO_line,'b','LineWidth',1.5);
                legend(str1,str2,'Location','NorthOutside','NumColumns',2);
                ylim([0 1.1*max(cat(1,ROI_AO_line,ROI_NoAO_line))]);
                text(0,max(cat(1,ROI_AO_line,ROI_NoAO_line)),['AO:' num2str(PVD_AO_curr,'%0.2f') ' NoAO: ' num2str(PVD_NoAO_curr,'%0.2f')...
                    ' , ' num2str(PVD_AO_curr./PVD_NoAO_curr)]);
                xlabel('t (s)'); ylabel('Signal (A.U.)');
                set(gcf,'Color',[1,1,1]); set(gca,'FontSize',16);
                grid off; set(gca,'LineWidth',1.5);

                resp_sz = menu('Keep the current d1_size, d2_size?','Save','d1_size+','d1_size-','d2_size+','d2_size-','Remove');
                if resp_sz==1
                    resp_bht = menu('Select','Brightness only','Brightness & contrast');
                    if resp_bht==1
                        labelN = 1;
                    elseif resp_bht==2
                        labelN = 2;
                    end
                    save([pplotT 'ROI_' num2str(iii) '_' num2str(MM) '.mat'],'XZ_NoAO','XZ_AO',...
                        's_final_AO','flatmethod','up_factor','ROI_AO','ROI_AO_up','ROI_AO_up_flat','ROI_AO_line',...
                        's_final_NoAO','ROI_NoAO','ROI_NoAO_up','ROI_NoAO_up_flat','ROI_NoAO_line',...
                        'tt_XZ','RMS_C_AO','RMS_C_NoAO','d1_size','d2_size','labelN');
                    for fig_num = 300
                        figure(fig_num);
                        saveas(gcf, [pplotT 'ROI_' num2str(iii) '_' num2str(MM) '_LP.fig']);
                    end
                    close(figure(300));
                    break;
                elseif resp_sz==2
                    d1_size = d1_size*1.1;
                elseif resp_sz==3
                    d1_size = d1_size/1.1;
                elseif resp_sz==4
                    d2_size = d2_size*1.1;
                elseif resp_sz==5
                    d2_size = d2_size/1.1;
                elseif resp_sz==6
                    MM = MM-1;
                    delete(hh2);
                    delete(hh1);
                    delete(txth);
                    close(figure(300));
                    break;
                end
            end

        else
            delete(hh2);
            delete(hh1);
        end
    end
    saveas(gcf, [pplotT 'ROI_' num2str(iii) '.fig']);
    close all;

    scaleInd = round(total_lengthAO/5/5)*5;
    if scaleInd>0
        [];
    else
        scaleInd = 5;
    end

    figure(100); set(gcf,'unit','centimeter','Position',[2 2 total_lengthAO/5 24]);
    imagesc(XZSlice_AO); colormap gray;
    colorbar; clim([Lim1 Lim2]); axis off; title([str1])
    set(gcf,'Color',[1,1,1]); set(gca,'FontSize',16);
    scalebar_mod2(total_lengthAO,0,scaleInd,'\mum','y','y');

    figure(200); set(gcf,'unit','centimeter','Position',[2 2 total_lengthAO/5 24]);
    imagesc(XZSlice_NoAO); colormap gray;
    colorbar; clim([Lim1 Lim2]); axis off; title([str2])
    set(gcf,'Color',[1,1,1]); set(gca,'FontSize',16);
    scalebar_mod2(total_lengthAO,0,scaleInd,'\mum','y','y');

    saveFigs_mod(pplotT,0,'.fig',['ROI_' num2str(iii)]);
    saveFigs_mod(pplotT,0,'.tif',['ROI_' num2str(iii)]);

    close all;
end

function isOval = localIsOvalROI(roi)
isOval = false;
typeFields = {'strType','Type','type','roiType'};
for kk = 1:numel(typeFields)
    if isfield(roi, typeFields{kk})
        roiTypeValue = roi.(typeFields{kk});
        if isstring(roiTypeValue) || ischar(roiTypeValue)
            roiTypeText = lower(char(roiTypeValue));
            if contains(roiTypeText,'oval') || contains(roiTypeText,'ellipse') || contains(roiTypeText,'ellipt')
                isOval = true;
                return;
            end
        elseif isnumeric(roiTypeValue)
            if any(roiTypeValue == [2 3])
                isOval = true;
                return;
            end
        end
    end
end
end

function mask = localBuildOvalMask(roi, nRows, nCols)
xy = [];
if isfield(roi,'mnCoordinates') && ~isempty(roi.mnCoordinates)
    xy = double(roi.mnCoordinates) + 1;
end

if size(xy,1) >= 3
    mask = poly2mask(xy(:,1), xy(:,2), nRows, nCols);
elseif isfield(roi,'vnRectBounds') && numel(roi.vnRectBounds) >= 4
    rectBounds = double(roi.vnRectBounds(:).');
    top = rectBounds(1) + 1;
    left = rectBounds(2) + 1;
    bottom = rectBounds(3) + 1;
    right = rectBounds(4) + 1;

    cx = (left + right) / 2;
    cy = (top + bottom) / 2;
    rx = max((right - left) / 2, 0.5);
    ry = max((bottom - top) / 2, 0.5);

    [xx,yy] = meshgrid(1:nCols,1:nRows);
    mask = ((xx - cx).^2 ./ rx.^2 + (yy - cy).^2 ./ ry.^2) <= 1;
elseif size(xy,1) == 2
    xMin = min(xy(:,1));
    xMax = max(xy(:,1));
    yMin = min(xy(:,2));
    yMax = max(xy(:,2));
    cx = (xMin + xMax) / 2;
    cy = (yMin + yMax) / 2;
    rx = max((xMax - xMin) / 2, 0.5);
    ry = max((yMax - yMin) / 2, 0.5);

    [xx,yy] = meshgrid(1:nCols,1:nRows);
    mask = ((xx - cx).^2 ./ rx.^2 + (yy - cy).^2 ./ ry.^2) <= 1;
else
    error('Unable to build an oval mask from the ROI. Inspect ROI fields from ReadImageJROI.');
end
end

function lineTrace = localMeanTraceFromMask(imgStack, mask)
mask = logical(mask);
if ~any(mask(:))
    error('Oval ROI mask is empty.');
end

nFrames = size(imgStack,3);
lineTrace = zeros(1,nFrames);
for frameIdx = 1:nFrames
    frameImg = imgStack(:,:,frameIdx);
    lineTrace(frameIdx) = mean(frameImg(mask), 'omitnan');
end
end

function localPlotMaskOutline(mask, lineColor)
B = bwboundaries(mask);
for bb = 1:numel(B)
    plot(B{bb}(:,2), B{bb}(:,1), '-', 'Color', lineColor, 'LineWidth', 2);
end
end

%% Export AO and NoAO videos
% Save separate videos for visual comparison.
croprg1 = 1:size(Img_AO,1); croprg2 = 1:size(Img_AO,2);
displayLim = [0 4500];

for strsave = {'AO','No AO'}
    close all;

    if strcmp(strsave,'AO')
        Img_volume = Img_AO;
    else strcmp(strsave,'No AO')
        Img_volume = Img_NoAO_reg;
    end
    Img_volume = Img_volume(croprg1,croprg2,:);
    FS = 1/dT; tx = [0:(Nframe-1)]*dT; tx = tx(:);
    x_dim = dX_FACED*size(Img_volume,2);
    y_dim = FOV_y*size(Img_volume,1);

    % Video, p2, scalebar size 50 um
    v = VideoWriter([pplotT 'vessel_' strsave{1} '.mp4'],'MPEG-4');
    v.Quality = 100;
    v.FrameRate = 80; %fps
    open(v); figure;
    for iii = 1:length(tx)
        iii
        imagesc(Img_volume(croprg1,croprg2,iii)); colormap gray
        pbaspect([x_dim y_dim 1]);
        colorbar; clim(displayLim); axis off;
        title([strsave{1} ', ' num2str(FS,'%0.1f') ' Hz, ' num2str(dep,'%0.1f') ' \mum, t = ' num2str(tx(iii),'%0.2f') ' s' ]);
        set(gcf,'Color',[1,1,1]); set(gca,'FontSize',16);
        scalebar_mod2(x_dim,0,50,'\mum','w','w',1,10,0);
        % Colormap_create(); colormap(redHotMap);
        % h.Color = [1 1 1]; h.Label.Color = [1 1 1];

        ax = gca;
        axPos = ax.Position;
        h = colorbar;
        ax.Position = axPos;

        imageRatio = x_dim / y_dim;
        axRatio = axPos(3) / axPos(4);
        if axRatio > imageRatio
            imageWidth = axPos(4) * imageRatio;
            imageRight = axPos(1) + (axPos(3) + imageWidth) / 2;
        else
            imageRight = axPos(1) + axPos(3);
        end

        colorbarGap = 0.02;
        colorbarWidth = 0.02;
        colorbarHeight = 0.2;
        h.Position = [imageRight + colorbarGap, 0.4, colorbarWidth, colorbarHeight]; % For large FOV data
        h.Ticks = displayLim; ylabel(h,'(A.U.)');
        h.Label.Units = 'normalized';
        h.Label.Position(1) = 1.8;

        mov = getframe(gcf);
        writeVideo(v,mov);
    end
    close(v);
end
save([pplotT 'video_para.mat'],'croprg1','croprg2','FS','tx','x_dim','y_dim');

%% Generate manuscript figures
% Plot selected ROIs, kymographs, and averaged line profiles.
pfigsave = [pathfolder 'Fig2\']; mkdir(pfigsave);
pdata = pplotT;
dy_FACED = FOV_y;
contrast_cal = 2;
Nratio = 0.2;

for ROIind = 14
    %% Prepare the selected ROI
    % Load its saved profiles and define the display region.
    close all;
    displayLim_XY = [0 1000]; displayLim_XZ = [0 5000];
    savethis = 1; savethismap = 1;
    ptsi = sROI_AO{1,ROIind}.mnCoordinates+1;

    load([pdata 'ROI_' num2str(ROIind) '.mat']);
    [cx,cy,~] = improfile(IAO,ptsi(:,1),ptsi(:,2),round(total_lengthAO*10),'bicubic');
    displayd1 = 1680:2350; displayd2 = 1:size(XZSlice_AO,2); % if showing part of the data

    ptsi_rg = abs(ptsi(end,:)-ptsi(1,:));
    dimen_rg = [ptsi_rg(1)*dX_FACED ptsi_rg(2)*dy_FACED 20*dX_FACED];
    [~,loc] = max(dimen_rg);
    if loc==1
        displaySz = round(ptsi_rg(1));
    elseif loc==2
        displaySz = round(ptsi_rg(2)/(dX_FACED/dy_FACED));
    elseif loc==3
        displaySz = 15;
    end
    ptsi_c = nanmean(ptsi);
    zoomd1 = round(ptsi_c(2)-displaySz*(dX_FACED/dy_FACED)):round(ptsi_c(2)+displaySz*(dX_FACED/dy_FACED));
    zoomd1 = intersect(zoomd1,1:size(IAO,1));
    zoomd2 = round(ptsi_c(1)-displaySz):round(ptsi_c(1)+displaySz);
    zoomd2 = intersect(zoomd2,1:size(IAO,2));

    sampleZoom = IAO(zoomd1,zoomd2);
    if savethismap==1
        if exist([pfigsave 'FOV_0002_ROImap.fig'], 'file') == 2
            uiopen([pfigsave 'FOV_0002_ROImap.fig'],1);
        else
            figure; imagesc(IAO); colormap gray
            x_dim = size(IAO,2)*dX_FACED; y_dim = size(IAO,1)*FOV_y;
            pbaspect([x_dim y_dim 1]);
            colorbar; clim(displayLim); axis off; title(['Mean IP: ' str1])
            set(gcf,'Color',[1,1,1]); set(gca,'FontSize',16);
            scalebar_mod2(x_dim,y_dim,50,'\mum','w','w');
            Colormap_create(); colormap(redHotMap);
            get(gca,'Position'); set(gca,'Position',ans);
            h = colorbar; h.Position = [0.85 0.1095 0.02 0.2]; % For large FOV data
            h.Ticks = displayLim;
        end
        drawROI_iGlu4(zoomd1,zoomd2,'w',num2str(ROIind));
        saveas(gcf,[pfigsave 'FOV_0002_ROImap.fig']);
    end
    close all;

    ptsi = sROI_NoAO{1,ROIind}.mnCoordinates+1;
    ptsi_c = nanmean(ptsi);
    zoomd1 = round(ptsi_c(2)-displaySz*(dX_FACED/dy_FACED)):round(ptsi_c(2)+displaySz*(dX_FACED/dy_FACED));
    zoomd1 = intersect(zoomd1,1:size(INoAO,1));
    zoomd2 = round(ptsi_c(1)-displaySz):round(ptsi_c(1)+displaySz);
    zoomd2 = intersect(zoomd2,1:size(INoAO,2));

    sysZoom = INoAO(zoomd1,zoomd2);
    if savethismap==1
        if exist([pfigsave 'FOV_0001_ROImap.fig'], 'file') == 2
            uiopen([pfigsave 'FOV_0001_ROImap.fig'],1);
        else
            figure; imagesc(INoAO); colormap gray
            pbaspect([x_dim y_dim 1]);
            colorbar; clim(displayLim); axis off; title(['Mean IP: ' str2])
            set(gcf,'Color',[1,1,1]); set(gca,'FontSize',16);
            scalebar_mod2(size(INoAO,2)*dX_FACED,size(INoAO,1)*FOV_y,50,'\mum','w','w');
            Colormap_create(); colormap(redHotMap);
            get(gca,'Position'); set(gca,'Position',ans);
            h = colorbar; h.Position = [0.85 0.1095 0.02 0.2]; % For large FOV data
            h.Ticks = displayLim;
        end
        drawROI_iGlu4(zoomd1,zoomd2,'w',num2str(ROIind));
        saveas(gcf,[pfigsave 'FOV_0001_ROImap.fig']);
    end
    close all;

    figure(100); imagesc(sampleZoom); colormap gray
    pbaspect([dX_FACED*size(sampleZoom,2) dy_FACED*size(sampleZoom,1) 1]);
    colorbar; clim(displayLim_XY); axis off; title(['Mean IP: ' str1])
    set(gcf,'Color',[1,1,1]); set(gca,'FontSize',16);
    x_dim = dX_FACED*size(sampleZoom,2);
    y_dim = dy_FACED*size(sampleZoom,1);
    scalebar_mod2(x_dim,y_dim,10,'\mum','w','w');
    Colormap_create(); colormap(redHotMap);
    get(gca,'Position'); set(gca,'Position',ans);
    h = colorbar; h.Position = [0.82 0.1095 0.02 0.2]; % For large FOV data
    h.Ticks = displayLim_XY;

    figure(200); imagesc(sysZoom); colormap gray
    pbaspect([dX_FACED*size(sysZoom,2) dy_FACED*size(sysZoom,1) 1]);
    hold on; plot(cx(displayd2)-zoomd2(1)+1,cy(displayd2)-zoomd1(1)+1,'c')
    colorbar; clim(displayLim_XY); axis off; title(['Mean IP: ' str2])
    set(gcf,'Color',[1,1,1]); set(gca,'FontSize',16);
    scalebar_mod2(x_dim,y_dim,10,'\mum','w','w');
    Colormap_create(); colormap(redHotMap);
    get(gca,'Position'); set(gca,'Position',ans);
    h = colorbar; h.Position = [0.82 0.1095 0.02 0.2]; % For large FOV data
    h.Ticks = displayLim_XY;

    %
    scaleInd = round(total_lengthAO/5/5)*5;
    if scaleInd>0
        [];
    else
        scaleInd = 5;
    end

    % scaleInd = 20
    figure(300);
    set(gcf,'unit','centimeter','Position',[2 2 8 24]);
    imagesc(XZSlice_AO(displayd1,displayd2)); colormap gray;
    colorbar; clim(displayLim_XZ); axis off; title([str1])
    set(gcf,'Color',[1,1,1]); set(gca,'FontSize',16);
    scalebar_mod2(length(displayd2)/size(XZSlice_AO,2)*total_lengthAO,0,scaleInd,'\mum','w','w');
    scalebar_mod2(0,length(displayd1)*dT,0.1,'s','w','w',0);
    Colormap_create(); colormap(redHotMap);
    set(gca,'Position',[0.1434    0.1100    0.5797    0.8150]);
    get(gca,'Position'); set(gca,'Position',ans);
    h = colorbar; h.Position = [0.75 0.1095 0.02 0.2]; % For large FOV data
    h.Ticks = displayLim_XZ;

    figure(400);
    set(gcf,'unit','centimeter','Position',[2 2 8 24]);
    imagesc(XZSlice_NoAO(displayd1,displayd2)); colormap gray;
    colorbar; clim(displayLim_XZ); axis off; title([str2])
    set(gcf,'Color',[1,1,1]); set(gca,'FontSize',16);
    scalebar_mod2(length(displayd2)/size(XZSlice_AO,2)*total_lengthAO,0,scaleInd,'\mum','w','w');
    scalebar_mod2(0,length(displayd1)*dT,0.1,'s','w','w',0);
    Colormap_create(); colormap(redHotMap);
    set(gca,'Position',[0.1434    0.1100    0.5797    0.8150]);
    get(gca,'Position'); set(gca,'Position',ans);
    h = colorbar; h.Position = [0.75 0.1095 0.02 0.2]; % For large FOV data
    h.Ticks = displayLim_XZ;

    if savethis
        save([pfigsave,'ROI' num2str(ROIind) '_display.mat'],'displayd1','displayd2','cx','cy');
    end

    ROI_line_file = dir([pdata 'ROI_' num2str(ROIind) '_*.mat']);
    XZ_AO_crop_prev = cell(length(ROI_line_file),1);
    XZ_NoAO_crop_prev = cell(length(ROI_line_file),1);
    d1_size_prev = nan(length(ROI_line_file),1);
    for iii = 1:length(ROI_line_file)
        clear ROI_AO ROI_AO_original ROI_NoAO ROI_NoAO_original ROI_AO_up ROI_AO_up_original ROI_NoAO_up ROI_NoAO_up_original
        load([pdata ROI_line_file(iii).name]);
        d1_size_original = d1_size;
        d1_size_plot = d1_size;
        ROI_AO_line_original = ROI_AO_line;
        ROI_NoAO_line_original = ROI_NoAO_line;
        if exist('ROI_AO','var')
            ROI_AO_original = ROI_AO;
        end
        if exist('ROI_NoAO','var')
            ROI_NoAO_original = ROI_NoAO;
        end
        if exist('ROI_AO_up','var')
            ROI_AO_up_original = ROI_AO_up;
        end
        if exist('ROI_NoAO_up','var')
            ROI_NoAO_up_original = ROI_NoAO_up;
        end

        update_d1_size_answer = 1;
        while update_d1_size_answer~=3
            ROI_line_N0 = length(ROI_AO_line_original);
            ROI_line_N = max(1,round(ROI_line_N0*d1_size_plot/d1_size_original));
            ROI_line_i1 = floor((ROI_line_N0-ROI_line_N)/2)+1;
            ROI_line_i2 = ROI_line_i1+ROI_line_N-1;
            ROI_AO_line = ROI_AO_line_original(ROI_line_i1:ROI_line_i2);
            ROI_NoAO_line = ROI_NoAO_line_original(ROI_line_i1:ROI_line_i2);
            ROI_raw_N = 2*d1_size_plot+1;
            if exist('ROI_AO_original','var')
                ROI_AO_N0 = size(ROI_AO_original,1);
                ROI_AO_i1 = max(1,floor((ROI_AO_N0-ROI_raw_N)/2)+1);
                ROI_AO_i2 = min(ROI_AO_N0,ROI_AO_i1+ROI_raw_N-1);
                ROI_AO = ROI_AO_original(ROI_AO_i1:ROI_AO_i2,:);
            end
            if exist('ROI_NoAO_original','var')
                ROI_NoAO_N0 = size(ROI_NoAO_original,1);
                ROI_NoAO_i1 = max(1,floor((ROI_NoAO_N0-ROI_raw_N)/2)+1);
                ROI_NoAO_i2 = min(ROI_NoAO_N0,ROI_NoAO_i1+ROI_raw_N-1);
                ROI_NoAO = ROI_NoAO_original(ROI_NoAO_i1:ROI_NoAO_i2,:);
            end
            if exist('ROI_AO_up_original','var')
                ROI_AO_up = ROI_AO_up_original(ROI_line_i1:ROI_line_i2,:);
            end
            if exist('ROI_NoAO_up_original','var')
                ROI_NoAO_up = ROI_NoAO_up_original(ROI_line_i1:ROI_line_i2,:);
            end

            XZ_AO_crop = [XZ_AO(:,1)-displayd2(1)+1,XZ_AO(:,2)-displayd1(1)+1];
            figure(300); clf;
            set(gcf,'unit','centimeter','Position',[2 2 8 24]);
            imagesc(XZSlice_AO(displayd1,displayd2)); colormap(redHotMap);
            colorbar; clim(displayLim_XZ); axis off; title([str1])
            set(gcf,'Color',[1,1,1]); set(gca,'FontSize',16);
            scalebar_mod2(length(displayd2)/size(XZSlice_AO,2)*total_lengthAO,0,scaleInd,'\mum','w','w');
            scalebar_mod2(0,length(displayd1)*dT,0.1,'s','w','w',0);
            set(gca,'Position',[0.1434    0.1100    0.5797    0.8150]);
            get(gca,'Position'); set(gca,'Position',ans);
            h = colorbar; h.Position = [0.75 0.1095 0.02 0.2];
            h.Ticks = displayLim_XZ;
            hold on;
            for jjj = 1:iii-1
                if ~isempty(XZ_AO_crop_prev{jjj})
                    plot(XZ_AO_crop_prev{jjj}(:,1),XZ_AO_crop_prev{jjj}(:,2),'w');
                    draw_rect(XZ_AO_crop_prev{jjj},d1_size_prev(jjj),d2_size,'w',num2str(jjj));
                end
            end
            plot(XZ_AO_crop(:,1),XZ_AO_crop(:,2),'w');
            draw_rect(XZ_AO_crop,d1_size_plot,d2_size,'w',num2str(iii));

            XZ_NoAO_crop = [XZ_NoAO(:,1)-displayd2(1)+1,XZ_NoAO(:,2)-displayd1(1)+1];
            figure(400); clf;
            set(gcf,'unit','centimeter','Position',[2 2 8 24]);
            imagesc(XZSlice_NoAO(displayd1,displayd2)); colormap(redHotMap);
            colorbar; clim(displayLim_XZ); axis off; title([str2])
            set(gcf,'Color',[1,1,1]); set(gca,'FontSize',16);
            scalebar_mod2(length(displayd2)/size(XZSlice_AO,2)*total_lengthAO,0,scaleInd,'\mum','w','w');
            scalebar_mod2(0,length(displayd1)*dT,0.1,'s','w','w',0);
            set(gca,'Position',[0.1434    0.1100    0.5797    0.8150]);
            get(gca,'Position'); set(gca,'Position',ans);
            h = colorbar; h.Position = [0.75 0.1095 0.02 0.2];
            h.Ticks = displayLim_XZ;
            hold on;
            for jjj = 1:iii-1
                if ~isempty(XZ_NoAO_crop_prev{jjj})
                    plot(XZ_NoAO_crop_prev{jjj}(:,1),XZ_NoAO_crop_prev{jjj}(:,2),'w');
                    draw_rect(XZ_NoAO_crop_prev{jjj},d1_size_prev(jjj),d2_size,'w',num2str(jjj));
                end
            end
            plot(XZ_NoAO_crop(:,1),XZ_NoAO_crop(:,2),'w');
            draw_rect(XZ_NoAO_crop,d1_size_plot,d2_size,'w',num2str(iii));

            [PVD_AO_curr,~,PVDMethod_AO_curr] = localPeakValleyDistinguishability_v2( ...
                ROI_AO_line,Nratio,'adaptive',dT/up_factor,slidingWindowSec,adaptiveFwhmThreshold);
            if PVDMethod_AO_curr=="sliding"
                NoAOMethod_curr = 'slidingaverage';
            else
                NoAOMethod_curr = char(PVDMethod_AO_curr);
            end
            [PVD_NoAO_curr,~,PVDMethod_NoAO_curr] = localPeakValleyDistinguishability_v2( ...
                ROI_NoAO_line,Nratio,NoAOMethod_curr,dT/up_factor,slidingWindowSec,adaptiveFwhmThreshold);
            figure(5000); clf; set(gcf,'Position',[100 390 560 260]);
            set(gca,'Position',[0.16 0.25 0.75 0.72]);
            up_factor = 8; tt_XZ = (0:(length(ROI_AO_line)-1))*dT/up_factor;
            plot(tt_XZ,ROI_AO_line,'r','LineWidth',1.5); hold on;
            plot(tt_XZ,ROI_NoAO_line,'b','LineWidth',1.5);
            [RMS_C_AO,RMS_C_NoAO] = VesselLineCst_v2(ROI_AO_line,ROI_NoAO_line,contrast_cal,Nratio);
            legend(str1,str2,'Location','NorthOutside','NumColumns',2);
            xlim([tt_XZ(1) tt_XZ(end)]);
            ylim([0 1.1*max(cat(1,ROI_AO_line,ROI_NoAO_line))]);
            text(0,max(cat(1,ROI_AO_line,ROI_NoAO_line)),['AO:' num2str(PVD_AO_curr,'%0.2f') ' NoAO: ' num2str(PVD_NoAO_curr,'%0.2f')...
                ' , ' num2str(PVD_AO_curr./PVD_NoAO_curr)]);
            xlabel('t (s)'); ylabel('Signal (A.U.)');
            set(gcf,'Color',[1,1,1]); set(gca,'FontSize',16);
            grid off; set(gca,'LineWidth',1.5);

            drawnow;
            update_d1_size_answer = menu(['d1_size = ' num2str(d1_size_plot) ' for ' ROI_line_file(iii).name], ...
                'd1_size/1.1','d1_size*1.1','Okay');
            if update_d1_size_answer==1
                d1_size_next = round(d1_size_plot/1.1);
                if d1_size_next<1
                    warning('Requested d1_size (%d) is smaller than the minimum reasonable value. Using d1_size = 1.', ...
                        d1_size_next);
                    d1_size_next = 1;
                end
                if d1_size_next==d1_size_plot
                    warning('Minimum d1_size is already used. Choose Okay to continue.');
                else
                    d1_size_plot = d1_size_next;
                end
            elseif update_d1_size_answer==2
                d1_size_next = round(d1_size_plot*1.1);
                if d1_size_next>d1_size_original
                    warning('Requested d1_size (%d) is larger than the initial value (%d). Maximum d1_size %d is used.', ...
                        d1_size_next,d1_size_original,d1_size_original);
                    d1_size_next = d1_size_original;
                end
                if d1_size_next==d1_size_plot
                    warning('Maximum d1_size is already used. Choose Okay to continue.');
                else
                    d1_size_plot = d1_size_next;
                end
            end
        end

        ROI_line_was_updated = d1_size_plot~=d1_size_original;
        if ROI_line_was_updated
            d1_size = d1_size_plot;
            [~,ROI_line_update_name,ROI_line_update_ext] = fileparts(ROI_line_file(iii).name);
            ROI_line_update_file = [pfigsave ROI_line_update_name '_update' ROI_line_update_ext];
            save(ROI_line_update_file,'ROI_AO_line','ROI_NoAO_line','tt_XZ','d1_size','d1_size_original','ROI_line_i1','ROI_line_i2','up_factor');
            if exist('ROI_AO','var')
                save(ROI_line_update_file,'ROI_AO','-append');
                save(ROI_line_update_file,'ROI_AO_i1','ROI_AO_i2','-append');
            end
            if exist('ROI_NoAO','var')
                save(ROI_line_update_file,'ROI_NoAO','-append');
                save(ROI_line_update_file,'ROI_NoAO_i1','ROI_NoAO_i2','-append');
            end
            if exist('ROI_AO_up','var')
                save(ROI_line_update_file,'ROI_AO_up','-append');
            end
            if exist('ROI_NoAO_up','var')
                save(ROI_line_update_file,'ROI_NoAO_up','-append');
            end
            d1_size = d1_size_original;
        end
        XZ_AO_crop_prev{iii} = XZ_AO_crop;
        XZ_NoAO_crop_prev{iii} = XZ_NoAO_crop;
        d1_size_prev(iii) = d1_size_plot;

        if savethis==1
            for figind = 5000 %[500:100:800]
                figure(figind);
                saveas(gcf,[pfigsave 'ROI' num2str(ROIind) '_' num2str(iii) '_LP.fig']);
                saveas(gcf,[pfigsave 'ROI' num2str(ROIind) '_' num2str(iii) '_LP.pdf']);
                close(figure(figind));
            end
        end
    end

    if savethis==1
        saveFigs_mod(pfigsave,0,'.fig',['ROI' num2str(ROIind)]);
        saveFigs_mod(pfigsave,0,'.pdf',['ROI' num2str(ROIind)]);
    end
end

%% Figure format
% Apply consistent formatting to the image panels.
Colormap_create(); colormap(redHotMap);
displayLim = [0 8000]; clim(displayLim);
get(gca,'Position'); set(gca,'Position',ans);
colorbar off;
h = colorbar; h.Position = [0.8 0.1095 0.02 0.2]; % For large FOV data
h.Ticks = displayLim;

%% Line profile format
% Apply consistent formatting to the line-profile panels.
set(gcf,'Position',[100 390 560 260]);
set(gca,'Position',[0.16 0.25 0.75 0.6]);
styleSavedProfileAxes(gca);

%% scalebar location move
% Adjust scale-bar placement for the final layout.
clear; clc;
lines = findobj(gcf, 'Type', 'line');
widths = zeros(1,length(lines));
colors = zeros(length(lines),3);
for iii = 1:length(lines)
    widths(iii) = lines(iii).LineWidth;
    colors(iii,:) = lines(iii).Color;
end

lines(4).YData = lines(4).YData+100;
lines(7).YData = lines(7).YData-100;

lines(3).YData = [lines(3).YData(1) lines(3).YData(1)+(lines(3).YData(2)-lines(3).YData(1))/10];

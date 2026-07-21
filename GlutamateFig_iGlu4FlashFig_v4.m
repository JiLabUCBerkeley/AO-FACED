clear; clc; close all;
pathfolder = 'H:\FACED_AO\20241217_M462_iGlu4\FACED\FOV2data9_VS\';
psave = [pathfolder 'ROI_update\']; cd(psave);
load([psave 'result.mat']);
FS = 1/dT; Nframe = size(traceAO,1); tx = [0:(Nframe-1)]*dT; Tavg = round(0.05*FS);
MoreMotionreg = 0;

Img_AO_reg = dataAO;
Img_NoAO_reg = dataNoAO;

str1 = 'sampleAO';
str2 = 'sysAO';
NFACED = 100;
dX_FACED = 0.8; %um
dy_FACED = FOV_y; %um

if (exist([psave 'pro.mat'],'file')==2)&(ismember('ROIvalid',who('-file', [psave 'pro.mat'])))
    load([psave 'pro.mat']);
else
    croprg1 = [161:800]; croprg2 = [6:137]; % data2
end

zrg = 1:size(Img_AO_reg,3);
sampleuse = mean(Img_AO_reg(croprg1,croprg2,zrg),3); % same plane for functional imaging
sysuse = mean(Img_NoAO_reg(croprg1,croprg2,zrg),3);

aopath = dir([pathfolder '*sampleAO']); aopath = [pathfolder aopath.name '\'];
tmp = dir([aopath '*Pkg']); aopath = [aopath tmp(1).name '\']
load([aopath 'RoiUpdate_sampleAO_seg_F0.mat'],'roimask');
roimaskAO = roimask(croprg1,croprg2,:);
noaopath = dir([pathfolder '*sysAO']); noaopath = [pathfolder noaopath.name '\'];
tmp = dir([noaopath '*Pkg']); noaopath = [noaopath tmp(1).name '\']
load([noaopath 'RoiUpdate_sysAO_seg_F0.mat'],'roimask');
roimaskNoAO = roimask(croprg1,croprg2,:);

if (exist([psave 'pro.mat'],'file')==2)&(ismember('ROIvalid',who('-file', [psave 'pro.mat'])))
    [];
else
    save([psave 'pro.mat'],'croprg1','croprg2','aopath','noaopath','roimaskAO','roimaskNoAO');
end

% options.overwrite = true;
% saveastiff(Img_AO_reg(croprg1,croprg2,zrg),[psave 'Img_AO_reg.tif'],options);
% saveastiff(Img_NoAO_reg(croprg1,croprg2,zrg),[psave 'Img_NoAO_reg.tif'],options);

%% Display the full field of view
% Show mean AO and NoAO images with matched intensity limits.
displayLim = [0 800];

close all;
figure(100); imagesc(sampleuse); colormap gray
pbaspect([dX_FACED*size(sampleuse,2) dy_FACED*size(sampleuse,1) 1]);
colorbar; clim(displayLim); axis off; title(['Mean IP: ' str1])
set(gcf,'Color',[1,1,1]); set(gca,'FontSize',16);
x_dim = dX_FACED*size(sampleuse,2);
y_dim = dy_FACED*size(sampleuse,1);
scalebar_mod(x_dim,y_dim,20,'\mum','w','w',2);
Colormap_create(); colormap(greenHotMap);
get(gca,'Position'); set(gca,'Position',ans);
h = colorbar; h.Position = [0.82 0.1095 0.02 0.2]; % For large FOV data
h.Ticks = [0 displayLim(2)];

figure(200); imagesc(sysuse); colormap gray
pbaspect([dX_FACED*size(sysuse,2) dy_FACED*size(sysuse,1) 1]);
colorbar; clim(displayLim); axis off; title(['Mean IP: ' str2])
set(gcf,'Color',[1,1,1]); set(gca,'FontSize',16);
scalebar_mod(x_dim,y_dim,20,'\mum','w','w',2);
Colormap_create(); colormap(greenHotMap);
get(gca,'Position'); set(gca,'Position',ans);
h = colorbar; h.Position = [0.82 0.1095 0.02 0.2]; % For large FOV data
h.Ticks = [0 displayLim(2)];

% saveFigs_mod(psave,1,'.fig','FOV_AllTrial');
% saveFigs_mod(psave,1,'.pdf','FOV_AllTrial');
% saveFigs_mod(psave,1,'.fig','FOV_oneTrial');
% saveFigs_mod(psave,1,'.pdf','FOV_oneTrial');

%% Configure trial timing
% Define stimulus timing and output settings for trial averaging.
staticMovingT = [1 0.05 1.95]; % iGlu4, since 12/08/2024
TimingCorr = 0;
pfigsave = [psave 'Compare_fig_20251029_trialavg\']; mkdir(pfigsave);

%% Segment and average trials
% Align repeated trials and calculate mean AO/NoAO responses.
[traceAO_Seg,tx_seg] = getdataSeg(traceAO,FS,sum(staticMovingT));
traceAO_trialAvg = squeeze(nanmean(nanmean(traceAO_Seg,2),4));
[traceNoAO_Seg,~] = getdataSeg(traceNoAO,FS,sum(staticMovingT));
traceNoAO_trialAvg = squeeze(nanmean(nanmean(traceNoAO_Seg,2),4));

traceAO_timeAvg = movmean(traceAO_Seg,Tavg);
data_perm = permute(traceAO_timeAvg, [1, 3, 2, 4]);
data_reshaped = reshape(data_perm,size(data_perm,1),size(data_perm,2),[]);
traceAO_trialAvg_timeAvg = nanmean(data_reshaped,3);
traceAO_trialAvg_timeAvg_se = nanstd(data_reshaped,[],3)/sqrt(size(data_reshaped,3));

traceNoAO_timeAvg = movmean(traceNoAO_Seg,Tavg);
data_perm = permute(traceNoAO_timeAvg, [1, 3, 2, 4]);
data_reshaped = reshape(data_perm,size(data_perm,1),size(data_perm,2),[]);
traceNoAO_trialAvg_timeAvg = nanmean(data_reshaped,3);
traceNoAO_trialAvg_timeAvg_se = nanstd(data_reshaped,[],3)/sqrt(size(data_reshaped,3));

%
BaseTimedFF = [0.5 0.5 3];
% BaseTimedFF(1): blank start
% BaseTimedFF(2): blank duration
% BaseTimedFF(3): cycle duration
angleNum = 1; % 1 repeat after averaging
[dF_F_AO_Seg] = getdFF_iGlu(traceAO_Seg,FS,BaseTimedFF,angleNum,size(traceAO_Seg,1));
dF_F_AO_trialAvg = squeeze(nanmean(nanmean(dF_F_AO_Seg,2),4));
[dF_F_NoAO_Seg] = getdFF_iGlu(traceNoAO_Seg,FS,BaseTimedFF,angleNum,size(traceAO_Seg,1));
dF_F_NoAO_trialAvg = squeeze(nanmean(nanmean(dF_F_NoAO_Seg,2),4));

dF_F_AO_timeAvg = movmean(dF_F_AO_Seg,Tavg);
data_perm = permute(dF_F_AO_timeAvg, [1, 3, 2, 4]);
data_reshaped = reshape(data_perm,size(data_perm,1),size(data_perm,2),[]);
dF_F_AO_trialAvg_timeAvg = nanmean(data_reshaped,3);
dF_F_AO_trialAvg_timeAvg_se = nanstd(data_reshaped,[],3)/sqrt(size(data_reshaped,3));

dF_F_NoAO_timeAvg = movmean(dF_F_NoAO_Seg,Tavg);
data_perm = permute(dF_F_NoAO_timeAvg, [1, 3, 2, 4]);
data_reshaped = reshape(data_perm,size(data_perm,1),size(data_perm,2),[]);
dF_F_NoAO_trialAvg_timeAvg = nanmean(data_reshaped,3);
dF_F_NoAO_trialAvg_timeAvg_se = nanstd(data_reshaped,[],3)/sqrt(size(data_reshaped,3));

%
% save([pfigsave 'result_2.mat'],'dataAO','traceAO','dataNoAO','traceNoAO',...
%     'traceAO_Seg','traceNoAO_Seg','traceAO_trialAvg','traceNoAO_trialAvg','dF_F_AO_trialAvg','dF_F_NoAO_trialAvg',...
%     'traceAO_trialAvg_timeAvg','dF_F_AO_trialAvg_timeAvg','traceNoAO_trialAvg_timeAvg','dF_F_NoAO_trialAvg_timeAvg',...
%     'Tavg','tx_seg','traceAO_timeAvg','traceNoAO_timeAvg','dF_F_AO_timeAvg','dF_F_NoAO_timeAvg',...
%     'traceAO_trialAvg_timeAvg_se','traceNoAO_trialAvg_timeAvg_se','dF_F_AO_trialAvg_timeAvg_se','dF_F_NoAO_trialAvg_timeAvg_se','-v7.3');

if TimingCorr
    for iii = 1:size(traceAO,3)
        subplot(2,1,1);
        plot(tx,mean(traceNoAO(:,roiTimeCorrUSE,iii),2)); xlim([0 12]); xticks([0:3:12])
        title(num2str(iii));
        subplot(2,1,2);
        plot(tx_seg,squeeze(nanmean(nanmean(traceNoAO_Seg(:,:,roiTimeCorrUSE,iii),2),3))); xlim([0 3]);
        input('Enter')
    end
end

%% Trace preview
% Preview the response trace for each valid ROI.
[negativeCorr,thre2,yLdff] = SensorProperty('iGlu4'); % 'VADER', 'JEDI2Pe', 'Force1A', 'JEDI2P-nonSoma'
FigorVideo = 4;
displaySz = 10;
displayLim = [-10 500];
ROIvalid = find(squeeze(nansum(nansum(roimaskAO)))>0).';
save([psave 'pro.mat'],'ROIvalid','-append');

ptrace = [psave 'tracePreview\']; mkdir(ptrace);
for ROIind = ROIvalid
    close all;
    figure; set(gcf,'Position',[1 41 1920 963]);
    ylimval = [0 max([max(max(traceAO(:,ROIind,:),[],1),[],3) max(max(traceNoAO(:,ROIind,:),[],1),[],3)])];

    subplot(4,1,1); plot(tx,traceAO(:,ROIind,1));
    xlim([0 12]); xticks([0:3:12]); set(gca,'XTickLabel',''); ylim(ylimval); title(num2str(ROIind));
    subplot(4,1,2); plot(tx,mean(traceAO(:,ROIind,:),3)); xlim([0 12]); xticks([0:3:12]);
    xlim([0 12]); xticks([0:3:12]); set(gca,'XTickLabel',''); ylim([0 ylimval(2)*0.5]);
    subplot(4,1,3); plot(tx,traceNoAO(:,ROIind,1)); xlim([0 12]); xticks([0:3:12]);
    xlim([0 12]); xticks([0:3:12]); set(gca,'XTickLabel',''); ylim(ylimval);
    subplot(4,1,4); plot(tx,mean(traceNoAO(:,ROIind,:),3));
    xlim([0 12]); xticks([0:3:12]); ylim([0 ylimval(2)*0.5]);

    saveFigs_mod(ptrace,1,'.fig',['ROI' num2str(ROIind) '_fig']);
    saveFigs_mod(ptrace,1,'.tif',['ROI' num2str(ROIind) '_fig']);
end

%% Comparison plot
% Generate paired AO/NoAO images, traces, or videos for each ROI.
for ROIind = ROIvalid
    temp = find(roimaskAO(:,:,ROIind)==1);
    [rr,cc] = ind2sub([size(roimaskAO,1),size(roimaskAO,2)],temp);
    zoomd1 = round(nanmean(rr))-displaySz*(dX_FACED/dy_FACED):round(nanmean(rr))+displaySz*(dX_FACED/dy_FACED);
    zoomd1 = intersect(zoomd1,1:size(sampleuse,1));
    zoomd2 = round(nanmean(cc))-displaySz:round(nanmean(cc))+displaySz;
    zoomd2 = intersect(zoomd2,1:size(sampleuse,2));
    if FigorVideo==1
        close all;
        sampleZoom = sampleuse(zoomd1,zoomd2);
        uiopen([psave 'FOV_0002_ROImap.fig'],1);
        drawROI_iGlu4(zoomd1,zoomd2,'w',num2str(ROIind));
        saveas(gcf,[psave 'FOV_0002_ROImap.fig']);
        close all;
    elseif FigorVideo==2
        sampleZoomVolume = Volume_AO_avg(zoomd1,zoomd2,:);
    elseif ismember(FigorVideo,[3,4])
        sampleZoom = sampleuse(zoomd1,zoomd2);
    end

    clear rr cc zoomd1 zoomd2
    temp = find(roimaskNoAO(:,:,ROIind)==1);
    [rr,cc] = ind2sub([size(roimaskNoAO,1),size(roimaskNoAO,2)],temp);
    zoomd1 = round(nanmean(rr))-displaySz*(dX_FACED/dy_FACED):round(nanmean(rr))+displaySz*(dX_FACED/dy_FACED);
    zoomd1 = intersect(zoomd1,1:size(sysuse,1));
    zoomd2 = round(nanmean(cc))-displaySz:round(nanmean(cc))+displaySz;
    zoomd2 = intersect(zoomd2,1:size(sysuse,2));
    if FigorVideo==1
        sysZoom = sysuse(zoomd1,zoomd2);
        uiopen([pathfolder 'FOV_0001_ROImap.fig'],1);
        drawROI_iGlu4(zoomd1,zoomd2,'w',num2str(ROIind));
        saveas(gcf,[pathfolder 'FOV_0001_ROImap.fig']);
        close all;
    elseif FigorVideo==2
        sysZoomVolume = Volume_NoAO_avg(zoomd1,zoomd2,:);
    elseif ismember(FigorVideo,[3,4])
        sysZoom = sysuse(zoomd1,zoomd2);
    end

    if FigorVideo==1
        close all;
        figure(100); imagesc(sampleZoom);
        pbaspect([dX_FACED*size(sampleZoom,2) dy_FACED*size(sampleZoom,1) 1]);
        colorbar; clim(displayLim); axis off; title(['Max IP: ' str1])
        set(gcf,'Color',[1,1,1]); set(gca,'FontSize',16);
        x_dim = dX_FACED*size(sampleZoom,2);
        y_dim = dy_FACED*size(sampleZoom,1);
        scalebar_mod(x_dim,y_dim,5,'\mum','w','w');
        Colormap_create(); colormap(greenHotMap);
        get(gca,'Position'); set(gca,'Position',ans);
        h = colorbar; h.Position = [0.8 0.1095 0.02 0.2]; % For large FOV data
        h.Ticks = displayLim;

        figure(200); imagesc(sysZoom);
        pbaspect([dX_FACED*size(sysZoom,2) dy_FACED*size(sysZoom,1) 1]);
        colorbar; clim(displayLim); axis off; title(['Max IP: ' str2])
        set(gcf,'Color',[1,1,1]); set(gca,'FontSize',16);
        scalebar_mod(x_dim,y_dim,5,'\mum','w','w');
        Colormap_create(); colormap(greenHotMap);
        get(gca,'Position'); set(gca,'Position',ans);
        h = colorbar; h.Position = [0.8 0.1095 0.02 0.2]; % For large FOV data
        h.Ticks = displayLim;

        saveFigs_mod(pfigsave,1,'.fig',['ROI' num2str(ROIind)]);
        saveFigs_mod(pfigsave,1,'.tif',['ROI' num2str(ROIind)]);

    elseif FigorVideo==2
        % v = VideoWriter([pfigsave 'ROI' num2str(ROIind) '.avi'],'Motion JPEG AVI');
        v = VideoWriter([pfigsave 'ROI' num2str(ROIind) '_v3.mp4'],'MPEG-4');
        v.Quality = 100;
        v.FrameRate = 35; %fps
        open(v);
        close all; figure; set(gcf,'unit','pixels','Position',[50 200 1200 500]);
        for iii = 1:length(tx_new)
            %%
            subplot('Position',[0.05 0.55 0.2 0.35]);
            hold off; imagesc(sysZoomVolume(:,:,iii));
            pbaspect([dX_FACED*size(sysZoomVolume,2) dy_FACED*size(sysZoomVolume,1) 1]);
            colorbar; clim(displayLim); axis off;
            text(1,5,'NoAO','Color','w','FontSize',16,'FontWeight','Bold');
            set(gcf,'Color',[1,1,1]); set(gca,'FontSize',16);
            x_dim = dX_FACED*size(sysZoomVolume,2);
            y_dim = dy_FACED*size(sysZoomVolume,1);
            scalebar_mod(x_dim,y_dim,5,'\mum','w','w');
            Colormap_create(); colormap(greenHotMap);
            get(gca,'Position'); set(gca,'Position',ans);
            h = colorbar; h.Position = [0.22 0.42 0.01 0.2]; % For large FOV data
            h.Ticks = displayLim;
            hText = findall(gca, 'Type', 'text');
            delete(hText(1));

            subplot('Position',[0.05 0.15 0.2 0.35]);
            hold off; imagesc(sampleZoomVolume(:,:,iii));
            pbaspect([dX_FACED*size(sampleZoomVolume,2) dy_FACED*size(sampleZoomVolume,1) 1]);
            colorbar; clim(displayLim); axis off;
            text(1,5,'AO','Color','w','FontSize',16,'FontWeight','Bold');
            set(gcf,'Color',[1,1,1]); set(gca,'FontSize',16);
            Colormap_create(); colormap(greenHotMap);
            pause(0.1);
            get(gca,'Position'); set(gca,'Position',ans);
            colorbar off;

            subplot('Position',[0.35 0.55 0.6 0.35]);
            hold off; plot(tx_new(1:iii),traceAO_trialAvg_timeAvg(1:iii,ROIind),'r'); hold on;
            plot(tx_new(1:iii),traceNoAO_avg(1:iii,ROIind),'b');
            YL = [0 1.1*max(traceAO_trialAvg_timeAvg(:,ROIind))];
            ylim(YL);
            title(['t = ' num2str(tx_new(iii),'%0.2f') ' s']);
            xticks([0:5:45]);
            xticklabels([]); xlim([0 48]); ylabel('F (A.U.)');
            legend('AO','NoAO','NumColumns',2,'AutoUpdate','off');
            set(gca,'FontSize',16,'LineWidth',2);
            AddStimuRect_fn(angleNum,staticMovingT,YL,[1,0.647,0],0.2);
            % line([tx_new(iii) tx_new(iii)],[YL],'Color','k');
            LineItems = findobj(gcf, 'Type', 'Line'); set(LineItems, 'LineWidth', 1.5);

            subplot('Position',[0.35 0.15 0.6 0.35]);
            hold off; plot(tx_new(1:iii),dF_F_AO_avg(1:iii,ROIind),'r'); hold on;
            plot(tx_new(1:iii),dF_F_NoAO_avg(1:iii,ROIind),'b');
            ylim(yLdff); xlim([0 48]); ylabel('\DeltaF/F'); xlabel('Time (s)');
            set(gca,'FontSize',16,'LineWidth',2);
            AddStimuRect_fn(angleNum,staticMovingT,yLdff,[1,0.647,0],0.2);
            % line([tx_new(iii) tx_new(iii)],[yLdff],'Color','k');
            LineItems = findobj(gcf, 'Type', 'Line'); set(LineItems, 'LineWidth', 1.5);

            %%
            mov = getframe(gcf);
            writeVideo(v,mov);
        end
        close(v);

    elseif FigorVideo==3
        %%
        displayLim = [0 round(2*max([sysZoom(:);sampleZoom(:)])/100)*100];

        close all; figure; set(gcf,'unit','pixels','Position',[50 200 1200 500]);
        subplot('Position',[0.05 0.55 0.2 0.35]);
        hold off; imagesc_auto(sysZoom);
        pbaspect([dX_FACED*size(sysZoom,2) dy_FACED*size(sysZoom,1) 1]);
        colorbar; clim(displayLim); axis off;
        text(1,5,'NoAO','Color','w','FontSize',16,'FontWeight','Bold');
        set(gcf,'Color',[1,1,1]); set(gca,'FontSize',16);
        x_dim = dX_FACED*size(sysZoom,2);
        y_dim = dy_FACED*size(sysZoom,1);
        scalebar_mod(x_dim,y_dim,5,'\mum','w','w');
        Colormap_create(); colormap(greenHotMap);
        get(gca,'Position'); set(gca,'Position',ans);
        h = colorbar; h.Position = [0.22 0.42 0.01 0.2]; % For large FOV data
        h.Ticks = [0 round(displayLim(2)/100)*100];

        subplot('Position',[0.05 0.15 0.2 0.35]);
        hold off; imagesc(sampleZoom);
        pbaspect([dX_FACED*size(sampleZoom,2) dy_FACED*size(sampleZoom,1) 1]);
        colorbar; clim(displayLim); axis off;
        text(1,5,'AO','Color','w','FontSize',16,'FontWeight','Bold');
        set(gcf,'Color',[1,1,1]); set(gca,'FontSize',16);
        Colormap_create(); colormap(greenHotMap);
        pause(0.1);
        get(gca,'Position'); set(gca,'Position',ans);
        colorbar off;

        temp = cat(1,mean(traceAO(:,ROIind,:),3),mean(traceNoAO(:,ROIind,:),3));
        YL = [0 max(temp(:))*1.05];

        subplot('Position',[0.35 0.55 0.25 0.35]);
        hold on; plot(tx_seg,traceAO_trialAvg(:,ROIind),'r');
        plot(tx_seg,traceNoAO_trialAvg(:,ROIind),'b');
        xlim([0 sum(staticMovingT)]); ylim(YL);
        ylabel('F (A.U.)'); box on;
        xticks([0:1:3]);
        xticklabels([]); xlim([0 sum(staticMovingT)]);
        legend('AO','NoAO','NumColumns',2,'AutoUpdate','off','Location','SouthEast');
        set(gca,'FontSize',16,'LineWidth',2);
        AddStimuRect_fn(1,staticMovingT,YL,[1,0.647,0],0.5);
        LineItems = findobj(gcf, 'Type', 'Line'); set(LineItems, 'LineWidth', 1.5);
        text(0.1,YL(2)*0.9,'Trial averaged','FontSize',16);

        subplot('Position',[0.7 0.55 0.25 0.35]);
        hold on; plot(tx_seg,dF_F_AO_trialAvg(:,ROIind),'r');
        plot(tx_seg,dF_F_NoAO_trialAvg(:,ROIind),'b');
        xlim([0 sum(staticMovingT)]); ylim(yLdff);
        ylabel('\DeltaF/F'); box on;
        xticks([0:1:3]);
        xticklabels([]); xlim([0 sum(staticMovingT)]);
        set(gca,'FontSize',16,'LineWidth',2);
        AddStimuRect_fn(1,staticMovingT,yLdff,[1,0.647,0],0.5);
        LineItems = findobj(gcf, 'Type', 'Line'); set(LineItems, 'LineWidth', 1.5);

        subplot('Position',[0.35 0.15 0.25 0.35]);
        hold on; plot(tx_seg,traceAO_trialAvg_timeAvg(:,ROIind),'r');
        plot(tx_seg,traceNoAO_trialAvg_timeAvg(:,ROIind),'b');
        xlim([0 sum(staticMovingT)]); ylim(YL);
        xlabel('Time (s)'); ylabel('F (A.U.)'); box on;
        xticks([0:1:3]);
        xlim([0 sum(staticMovingT)]);
        set(gca,'FontSize',16,'LineWidth',2);
        AddStimuRect_fn(1,staticMovingT,YL,[1,0.647,0],0.5);
        LineItems = findobj(gcf, 'Type', 'Line'); set(LineItems, 'LineWidth', 1.5);
        text(0.1,YL(2)*0.9,'Trial & time averaged','FontSize',16);

        subplot('Position',[0.7 0.15 0.25 0.35]);
        hold on; plot(tx_seg,dF_F_AO_trialAvg_timeAvg(:,ROIind),'r');
        plot(tx_seg,dF_F_NoAO_trialAvg_timeAvg(:,ROIind),'b');
        xlim([0 sum(staticMovingT)]); ylim(yLdff);
        xlabel('Time (s)'); ylabel('\DeltaF/F'); box on;
        xticks([0:1:3]);
        xlim([0 sum(staticMovingT)]);
        set(gca,'FontSize',16,'LineWidth',2);
        AddStimuRect_fn(1,staticMovingT,yLdff,[1,0.647,0],0.5);
        LineItems = findobj(gcf, 'Type', 'Line'); set(LineItems, 'LineWidth', 1.5);
        % saveFigs_mod(pfigsave,1,'.fig',['ROI' num2str(ROIind) '_fig']);
        % saveFigs_mod(pfigsave,1,'.tif',['ROI' num2str(ROIind) '_fig']);

        % saveFigs_mod(pfigsave,1,'.fig',['ROI' num2str(ROIind) '_fig_adj']);
        % saveFigs_mod(pfigsave,1,'.tif',['ROI' num2str(ROIind) '_fig_adj']);
        % saveFigs_mod(pfigsave,1,'.pdf',['ROI' num2str(ROIind) '_fig_adj']);

    elseif FigorVideo==4
        %%
        displayLim = [0 round(2*max([sysZoom(:);sampleZoom(:)])/100)*100];

        close all; figure; set(gcf,'unit','pixels','Position',[50 200 1200 500]);
        subplot('Position',[0.05 0.55 0.2 0.35]);
        hold off; imagesc_auto(sysZoom);
        pbaspect([dX_FACED*size(sysZoom,2) dy_FACED*size(sysZoom,1) 1]);
        colorbar; clim(displayLim); axis off;
        text(1,5,'NoAO','Color','w','FontSize',16,'FontWeight','Bold');
        set(gcf,'Color',[1,1,1]); set(gca,'FontSize',16);
        x_dim = dX_FACED*size(sysZoom,2);
        y_dim = dy_FACED*size(sysZoom,1);
        scalebar_mod(x_dim,y_dim,5,'\mum','w','w');
        Colormap_create(); colormap(greenHotMap);
        get(gca,'Position'); set(gca,'Position',ans);
        h = colorbar; h.Position = [0.22 0.42 0.01 0.2]; % For large FOV data
        h.Ticks = [0 round(displayLim(2)/100)*100];

        subplot('Position',[0.05 0.15 0.2 0.35]);
        hold off; imagesc(sampleZoom);
        pbaspect([dX_FACED*size(sampleZoom,2) dy_FACED*size(sampleZoom,1) 1]);
        colorbar; clim(displayLim); axis off;
        text(1,5,'AO','Color','w','FontSize',16,'FontWeight','Bold');
        set(gcf,'Color',[1,1,1]); set(gca,'FontSize',16);
        Colormap_create(); colormap(greenHotMap);
        pause(0.1);
        get(gca,'Position'); set(gca,'Position',ans);
        colorbar off;

        temp = cat(1,mean(traceAO(:,ROIind,:),3),mean(traceNoAO(:,ROIind,:),3));
        YL = [0 max(temp(:))*1.05];

        YL = [0 150]; yLdff = [-0.5 1.0];

        subplot('Position',[0.35 0.15 0.25 0.35]); hold on;
        [h,j]=alphaerrorplot(tx_seg,traceAO_trialAvg_timeAvg(:,ROIind),traceAO_trialAvg_timeAvg_se(:,ROIind),0.15,'r');
        set(h,'Marker','none'); set(j,'LineWidth',2,'Color','r');
        [h,j]=alphaerrorplot(tx_seg,traceNoAO_trialAvg_timeAvg(:,ROIind),traceNoAO_trialAvg_timeAvg_se(:,ROIind),0.15,'b');
        set(h,'Marker','none'); set(j,'LineWidth',2,'Color','b');
        xlim([0 sum(staticMovingT)]); ylim(YL);
        xlabel('Time (s)'); ylabel('F (A.U.)'); box on;
        xticks([0:1:3]);
        xlim([0 sum(staticMovingT)]);
        set(gca,'FontSize',16,'LineWidth',2);
        AddStimuRect_fn(1,staticMovingT,YL,[1,0.647,0],0.5);
        LineItems = findobj(gcf, 'Type', 'Line'); set(LineItems, 'LineWidth', 1.5);
        text(0.1,YL(2)*0.9,'Trial & time averaged','FontSize',16);

        subplot('Position',[0.7 0.15 0.25 0.35]); hold on;
        [h,j]=alphaerrorplot(tx_seg,dF_F_AO_trialAvg_timeAvg(:,ROIind),dF_F_AO_trialAvg_timeAvg_se(:,ROIind),0.15,'r');
        set(h,'Marker','none'); set(j,'LineWidth',2,'Color','r');
        [h,j]=alphaerrorplot(tx_seg,dF_F_NoAO_trialAvg_timeAvg(:,ROIind),dF_F_NoAO_trialAvg_timeAvg_se(:,ROIind),0.15,'b');
        set(h,'Marker','none'); set(j,'LineWidth',2,'Color','b');
        xlim([0 sum(staticMovingT)]); ylim(yLdff);
        xlabel('Time (s)'); ylabel('\DeltaF/F'); box on;
        xticks([0:1:3]);
        xlim([0 sum(staticMovingT)]);
        set(gca,'FontSize',16,'LineWidth',2);
        AddStimuRect_fn(1,staticMovingT,yLdff,[1,0.647,0],0.5);
        LineItems = findobj(gcf, 'Type', 'Line'); set(LineItems, 'LineWidth', 1.5);

        % saveFigs_mod(pfigsave,1,'.fig',['ROI' num2str(ROIind) '_fig']);
        % saveFigs_mod(pfigsave,1,'.tif',['ROI' num2str(ROIind) '_fig']);

        % saveFigs_mod(pfigsave,1,'.fig',['ROI' num2str(ROIind) '_fig_adj']);
        % saveFigs_mod(pfigsave,1,'.tif',['ROI' num2str(ROIind) '_fig_adj']);
        % saveFigs_mod(pfigsave,1,'.pdf',['ROI' num2str(ROIind) '_fig_adj']);
    end
end

%% signal comparison
% Compare AO and NoAO response metrics across valid ROIs.
% null - mean is m
% 'left' — Test against the alternative hypothesis that the population mean is less than m.
% h = ttest(x,y) returns a test decision for the null hypothesis that the data in x – y comes from a normal distribution with mean equal to zero and unknown variance, using the paired-sample t-test.
% h = 1, rejects null, accepts alternative

Oneplot = 1;
NormToNoAO = 1;
close all;
stdorSE = 2; %1-std; 2-SE

% ylimVal = [0 300]; % FOV1data1, FOV1data6
% ylimVal = [0 550]; % FOV1data2
% ylimVal = [0 220]; % FOV1data8
% ylimVal = [0 400]; % FOV2data4, FOV2data9
% ylimVal = [0 500]; % FOV2data5, FOV2data2, FOV1data9
% ylimVal = [0 650]; % FOV2data6
% ylimVal = [0 350]; % FOV2data10
% ylimVal = [0 600]; % FOV2data11, FOV2data8, FOV1data7

% ylimVal = [0 450]; % 0313 data data1_VS, data2_VS

if NormToNoAO==1
    ylimVal = [0 2];
    staLine = ylimVal(1)+(ylimVal(2)-ylimVal(1))*0.9; staLine = [staLine staLine];
    staStar = mean(staLine)+(ylimVal(2)-ylimVal(1))*0.03;
    yTickVal = ylimVal(1):0.5:ylimVal(2);
else
    staLine = ylimVal(1)+(ylimVal(2)-ylimVal(1))*0.9; staLine = [staLine staLine];
    staStar = mean(staLine)+(ylimVal(2)-ylimVal(1))*0.03;
    yTickVal = ylimVal(1):100:ylimVal(2);
end

ROIuseInd = ROIvalid;
baselineLim = [0.5 1];
t_baseline = find((tx_seg>=baselineLim(1)) & (tx_seg<baselineLim(2)));
temp = mean(traceAO_trialAvg_timeAvg(t_baseline,ROIuseInd),1);
AA = find((temp<6*mean(temp))==1);
ROIuseInd = ROIuseInd(AA);
save([psave 'pro.mat'],'ROIuseInd','-append');

fplotNoAOUSE = traceNoAO_trialAvg;
fplotAOUSE = traceAO_trialAvg;

plotCNoAO = mean(fplotNoAOUSE(t_baseline,ROIuseInd),1);
plotCAO = mean(fplotAOUSE(t_baseline,ROIuseInd),1);
[plotCNoAO_sti_norm,plotCAO_sti_norm,h,p] = plot_compare(NormToNoAO,plotCNoAO,plotCAO,0,...
    Oneplot,1,stdorSE);
xlim([0.5,2.5]); ylabel('Signal EF');
xticks([1,2]); xticklabels({'No AO','AO'}); yticks(yTickVal);
box on; title(['Pre stimulation, ' num2str(baselineLim(1)-1) ' - 0 s']); ylim(ylimVal)
statistic_plot([1 2],staLine,staStar,p);

% post stimulation, 1-2s
t_sti = find((tx_seg>=1) & (tx_seg<2));
plotCNoAO_sti = mean(fplotNoAOUSE(t_sti,ROIuseInd),1);
plotCAO_sti = mean(fplotAOUSE(t_sti,ROIuseInd),1);
[plotCNoAO_sti_norm,plotCAO_sti_norm,h,p] = plot_compare(NormToNoAO,plotCNoAO_sti,plotCAO_sti,0,...
    Oneplot,2,stdorSE);
xlim([0.5,2.5]); ylabel('Signal EF');
xticks([1,2]); xticklabels({'No AO','AO'}); yticks(yTickVal);
box on; title(['Post stimulation, 0 - 1 s']); ylim(ylimVal)
statistic_plot([1 2],staLine,staStar,p);

% post stimulation, 2-3s
t_sti = find((tx_seg>=2) & (tx_seg<3));
plotCNoAO_sti = mean(fplotNoAOUSE(t_sti,ROIuseInd),1);
plotCAO_sti = mean(fplotAOUSE(t_sti,ROIuseInd),1);
[plotCNoAO_sti_norm,plotCAO_sti_norm,h,p] = plot_compare(NormToNoAO,plotCNoAO_sti,plotCAO_sti,0,...
    Oneplot,3,stdorSE);
xlim([0.5,2.5]); ylabel('Signal EF');
xticks([1,2]); xticklabels({'No AO','AO'}); yticks(yTickVal);
box on; title(['Post stimulation, 1 - 2 s']); ylim(ylimVal)
statistic_plot([1 2],staLine,staStar,p);

% post stimulation, 1-3s
t_sti = find((tx_seg>=1) & (tx_seg<3));
plotCNoAO_sti = mean(fplotNoAOUSE(t_sti,ROIuseInd),1);
plotCAO_sti = mean(fplotAOUSE(t_sti,ROIuseInd),1);
[plotCNoAO_sti_norm,plotCAO_sti_norm,h,p] = plot_compare(NormToNoAO,plotCNoAO_sti,plotCAO_sti,0,...
    Oneplot,4,stdorSE);
xlim([0.5,2.5]); ylabel('Signal EF');
xticks([1,2]); xticklabels({'No AO','AO'}); yticks(yTickVal);
box on; title(['Post stimulation, 0 - 2 s']); ylim(ylimVal)
statistic_plot([1 2],staLine,staStar,p);

% if Oneplot==1
%     saveFigs_mod(pfigsave,1,'.fig',['Summary_F_combine_norm']);
%     saveFigs_mod(pfigsave,1,'.tif',['Summary_F_combine_norm']);
%     saveFigs_mod(pfigsave,1,'.pdf',['Summary_F_combine_norm']);
% else
%     saveFigs_mod(pfigsave,1,'.fig',['Summary_F_norm']);
%     saveFigs_mod(pfigsave,1,'.tif',['Summary_F_norm']);
%     saveFigs_mod(pfigsave,1,'.pdf',['Summary_F_norm']);
% end

%% dF/F change
% Quantify and compare stimulus-evoked dF/F changes.
dffthre = 0.10;
if Oneplot==1
    savestr = ['Summary_dff_combine_dffthre' num2str(dffthre)];
else
    savestr = 'Summary_dff';
end
close all;
% ylimVal = [-0.1 0.35]; % FOV1data1, FOV1data7, FOV2data6
% ylimVal = [-0.1 0.25]; % FOV2data9
% ylimVal = [-0.1 0.5]; % FOV2data11, FOV2data8
% ylimVal = [-0.1 0.45]; % FOV2data10, FOV2data5
% ylimVal = [-0.1 0.2]; % FOV1data2, FOV1data6, FOV1data9
% ylimVal = [-0.1 0.4]; % FOV1data8
% ylimVal = [-0.1 0.6]; % FOV2data4

% ylimVal = [-0.1 0.3]; % 0313 data data1_VS, data2_VS

if NormToNoAO==1
    ylimVal = [0 3];
    staLine = ylimVal(1)+(ylimVal(2)-ylimVal(1))*0.9; staLine = [staLine staLine];
    staStar = mean(staLine)+(ylimVal(2)-ylimVal(1))*0.03;
    yTickVal = ylimVal(1):1:ylimVal(2);
else
    staLine = ylimVal(1)+(ylimVal(2)-ylimVal(1))*0.9; staLine = [staLine staLine];
    staStar = mean(staLine)+(ylimVal(2)-ylimVal(1))*0.03;
    yTickVal = ylimVal(1):0.1:ylimVal(2);
end

dffplotNoAOUSE = dF_F_NoAO_trialAvg;
dffplotAOUSE = dF_F_AO_trialAvg;

% t_sti = find((tx_seg>=1) & (tx_seg<3));
% plotdffNoAO_sti = mean(dffplotNoAOUSE(t_sti,:),1);
% plotdffAO_sti = mean(dffplotAOUSE(t_sti,:),1);
% Indtemp1 = find((plotdffNoAO_sti>0.1)==1);
% Indtemp2 = find((plotdffAO_sti>0.1)==1);
% % t_sti = find((tx_seg>=1) & (tx_seg<3));
% % plotdffAO_sti = mean(dffplotAOUSE(t_sti,:),1);
% % Indtemp2 = find((plotdffAO_sti>0.1)==1);
% ROIuseInd = unique([Indtemp2]);
% savestr = 'Summary_dff_combine_AOactive_thre0.1';

%
plotdffNoAO_sti = mean(dffplotNoAOUSE(t_baseline,ROIuseInd),1);
plotdffAO_sti = mean(dffplotAOUSE(t_baseline,ROIuseInd),1);
[dff_NoAO_norm,dff_AO_norm,h,p] = plot_compare(NormToNoAO,plotdffNoAO_sti,plotdffAO_sti,dffthre,...
    Oneplot,1,stdorSE);
xlim([0.5,2.5]); ylabel('\DeltaF/F EF');
xticks([1,2]); xticklabels({'No AO','AO'}); yticks(yTickVal);
box on; title(['Pre stimulation, ' num2str(baselineLim(1)-1) ' - 0 s']); ylim(ylimVal)
statistic_plot([1 2],staLine,staStar,p);

% post stimulation, 1-2s
t_sti = find((tx_seg>=1) & (tx_seg<2));
plotdffNoAO_sti = mean(dffplotNoAOUSE(t_sti,ROIuseInd),1);
plotdffAO_sti = mean(dffplotAOUSE(t_sti,ROIuseInd),1);
[dff_NoAO_norm,dff_AO_norm,h,p] = plot_compare(NormToNoAO,plotdffNoAO_sti,plotdffAO_sti,dffthre,...
    Oneplot,2,stdorSE);
xlim([0.5,2.5]); ylabel('\DeltaF/F EF');
xticks([1,2]); xticklabels({'No AO','AO'}); yticks(yTickVal);
box on; title(['Post stimulation, 0 - 1 s']); ylim(ylimVal)
statistic_plot([1 2],staLine,staStar,p);

% post stimulation, 2-3s
t_sti = find((tx_seg>=2) & (tx_seg<3));
plotdffNoAO_sti = mean(dffplotNoAOUSE(t_sti,ROIuseInd),1);
plotdffAO_sti = mean(dffplotAOUSE(t_sti,ROIuseInd),1);
[dff_NoAO_norm,dff_AO_norm,h,p] = plot_compare(NormToNoAO,plotdffNoAO_sti,plotdffAO_sti,dffthre,...
    Oneplot,3,stdorSE);
xlim([0.5,2.5]); ylabel('\DeltaF/F EF');
xticks([1,2]); xticklabels({'No AO','AO'}); yticks(yTickVal);
box on; title(['Post stimulation, 1 - 2 s']); ylim(ylimVal)
statistic_plot([1 2],staLine,staStar,p);

% post stimulation, 1-3 s
t_sti = find((tx_seg>=1) & (tx_seg<3));
plotdffNoAO_sti = mean(dffplotNoAOUSE(t_sti,ROIuseInd),1);
plotdffAO_sti = mean(dffplotAOUSE(t_sti,ROIuseInd),1);
[dff_NoAO_norm,dff_AO_norm,h,p] = plot_compare(NormToNoAO,plotdffNoAO_sti,plotdffAO_sti,dffthre,...
    Oneplot,4,stdorSE);
xlim([0.5,2.5]); ylabel('\DeltaF/F EF');
xticks([1,2]); xticklabels({'No AO','AO'}); yticks(yTickVal);
box on; title(['Post stimulation, 0 - 2 s']); ylim(ylimVal)
statistic_plot([1 2],staLine,staStar,p);

% saveFigs_mod(pfigsave,1,'.fig',savestr);
% saveFigs_mod(pfigsave,1,'.tif',savestr);
% saveFigs_mod(pfigsave,1,'.pdf',savestr);

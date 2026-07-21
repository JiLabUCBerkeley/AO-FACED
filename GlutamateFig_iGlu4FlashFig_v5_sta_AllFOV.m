%%
clear; clc; close all;
folderInd = [8:11];
NormToNoAO = 1;
meanormax = 1; %1-mean; 2-max. Use mean during the visual stimulation period for comparison
Displayratio = 0;
trialavgOnly = 1; % 1- using trial-avg data; 2- using trial&time-avg data
pathroot = 'H:\FACED_AO\20241217_M462_iGlu4\FACED\';

Validthre = 0.35 % if the negative value fraction is less than this threshold, then the ROI is valid
dffthre = 0.02 %

pfigsave = 'H:\FACED_AO\20241217_M462_iGlu4\FACED\FOV2data8   9  10  11_VS_All\ROI_update\Compare_fig_202601_trialavg\';
mkdir(pfigsave); cd(pfigsave);

if Displayratio==1
    load([pathroot 'FOV2data' num2str(folderInd) '_VS\ROI_update\result.mat'],'FOV_y');
    load([pathroot 'FOV2data' num2str(folderInd) '_VS\ROI_update\pro.mat']);
    uiopen([pathroot 'FOV2data' num2str(folderInd) '_VS\ROI_update\FOV_AllTrial_0001.fig'],1); hold on;
    IAO = tiffreadVolume([aopath 'F' num2str(0) '_RM_avg_ref1.tif']);
    [~,xy]=readImageJROI_main_mod(aopath,['RoiUpdate_sampleAO_seg_F' num2str(0) '.zip'],IAO); bds = xy.';
    close(figure(31));
    IAO = IAO(croprg1,croprg2);
end

plotFNoAO_base_all = [];
plotFAO_base_all = [];
plotFNoAO_1_2_all = [];
plotFAO_1_2_all = [];
plotFNoAO_2_3_all = [];
plotFAO_2_3_all = [];
plotFNoAO_1_3_all = [];
plotFAO_1_3_all = [];

plotdffNoAO_base_all = [];
plotdffAO_base_all = [];
plotdffNoAO_1_2_all = [];
plotdffAO_1_2_all = [];
plotdffNoAO_2_3_all = [];
plotdffAO_2_3_all = [];
plotdffNoAO_1_3_all = [];
plotdffAO_1_3_all = [];

for iii = folderInd % Different datasets have different acquisition speeds
    pathfolder = ['H:\FACED_AO\20241217_M462_iGlu4\FACED\FOV2data' num2str(iii) '_VS\ROI_update\']
    cd(pathfolder)

    baselineLim = [0.5 1]; % [0.5 1]
    t_activerg = [1 2]; % 0-1 s post-stimulation, which is 1-2 s in the whole data
    VEthre = 3; % N times std

    if trialavgOnly==1
        load([pathfolder 'result_2.mat'],'traceNoAO_trialAvg','traceAO_trialAvg','dF_F_NoAO_trialAvg','dF_F_AO_trialAvg','tx_seg');
        FROIuse = 1:size(traceNoAO_trialAvg,2);
        [ValidInd, VEind] = isiGlu4VE(dF_F_AO_trialAvg,dF_F_NoAO_trialAvg,tx_seg,baselineLim,t_activerg,Validthre,VEthre);
        length(ValidInd)/size(traceNoAO_trialAvg,2)
    elseif trialavgOnly==2
        load([pathfolder 'result_2.mat'],'traceNoAO_trialAvg_timeAvg','traceAO_trialAvg_timeAvg','dF_F_NoAO_trialAvg_timeAvg','dF_F_AO_trialAvg_timeAvg','tx_seg');
        FROIuse = 1:size(traceNoAO_trialAvg_timeAvg,2);
        [ValidInd, VEind] = isiGlu4VE(dF_F_AO_trialAvg_timeAvg,dF_F_NoAO_trialAvg_timeAvg,tx_seg,baselineLim,t_activerg,Validthre,VEthre);
        length(ValidInd)/size(traceNoAO_trialAvg_timeAvg,2)
    end

    %%
    if trialavgOnly==1
        fplotNoAOUSE = traceNoAO_trialAvg;
        fplotAOUSE = traceAO_trialAvg;
        dffplotNoAOUSE = dF_F_NoAO_trialAvg;
        dffplotAOUSE = dF_F_AO_trialAvg;
    elseif trialavgOnly==2
        fplotNoAOUSE = traceNoAO_trialAvg_timeAvg;
        fplotAOUSE = traceAO_trialAvg_timeAvg;
        dffplotNoAOUSE = dF_F_NoAO_trialAvg_timeAvg;
        dffplotAOUSE = dF_F_AO_trialAvg_timeAvg;
    end

    FROIuse = intersect(ValidInd,VEind);
    VEind = intersect(ValidInd,VEind);

    t_baseline = find((tx_seg>=baselineLim(1)) & (tx_seg<baselineLim(2)));
    plotFNoAO_base = mean(fplotNoAOUSE(t_baseline,FROIuse),1);
    plotFAO_base = mean(fplotAOUSE(t_baseline,FROIuse),1);
    if meanormax==1
        plotdffNoAO_base = mean(dffplotNoAOUSE(t_baseline,VEind),1);
        plotdffAO_base = mean(dffplotAOUSE(t_baseline,VEind),1);
    elseif meanormax==2
        plotdffNoAO_base = max(dffplotNoAOUSE(t_baseline,VEind),[],1);
        plotdffAO_base = max(dffplotAOUSE(t_baseline,VEind),[],1);
    end

    t_sti = find((tx_seg>=1) & (tx_seg<2));
    plotFNoAO_1_2 = mean(fplotNoAOUSE(t_sti,FROIuse),1);
    plotFAO_1_2 = mean(fplotAOUSE(t_sti,FROIuse),1);
    if meanormax==1
        plotdffNoAO_1_2 = mean(dffplotNoAOUSE(t_sti,VEind),1);
        plotdffAO_1_2 = mean(dffplotAOUSE(t_sti,VEind),1);
    elseif meanormax==2
        plotdffNoAO_1_2 = max(dffplotNoAOUSE(t_sti,VEind),[],1);
        plotdffAO_1_2 = max(dffplotAOUSE(t_sti,VEind),[],1);
    end

    t_sti = find((tx_seg>=2) & (tx_seg<3));
    plotFNoAO_2_3 = mean(fplotNoAOUSE(t_sti,FROIuse),1);
    plotFAO_2_3 = mean(fplotAOUSE(t_sti,FROIuse),1);
    if meanormax==1
        plotdffNoAO_2_3 = mean(dffplotNoAOUSE(t_sti,VEind),1);
        plotdffAO_2_3 = mean(dffplotAOUSE(t_sti,VEind),1);
    elseif meanormax==2
        plotdffNoAO_2_3 = max(dffplotNoAOUSE(t_sti,VEind),[],1);
        plotdffAO_2_3 = max(dffplotAOUSE(t_sti,VEind),[],1);
    end

    t_sti = find((tx_seg>=1) & (tx_seg<3));
    plotFNoAO_1_3 = mean(fplotNoAOUSE(t_sti,FROIuse),1);
    plotFAO_1_3 = mean(fplotAOUSE(t_sti,FROIuse),1);
    if meanormax==1
        plotdffNoAO_1_3 = mean(dffplotNoAOUSE(t_sti,VEind),1);
        plotdffAO_1_3 = mean(dffplotAOUSE(t_sti,VEind),1);
    elseif meanormax==2
        plotdffNoAO_1_3 = max(dffplotNoAOUSE(t_sti,VEind),[],1);
        plotdffAO_1_3 = max(dffplotAOUSE(t_sti,VEind),[],1);
    end

    plotFNoAO_base_all = [plotFNoAO_base_all plotFNoAO_base];
    plotFAO_base_all = [plotFAO_base_all plotFAO_base];
    plotFNoAO_1_2_all = [plotFNoAO_1_2_all plotFNoAO_1_2];
    plotFAO_1_2_all = [plotFAO_1_2_all plotFAO_1_2];
    plotFNoAO_2_3_all = [plotFNoAO_2_3_all plotFNoAO_2_3];
    plotFAO_2_3_all = [plotFAO_2_3_all plotFAO_2_3];
    plotFNoAO_1_3_all = [plotFNoAO_1_3_all plotFNoAO_1_3];
    plotFAO_1_3_all = [plotFAO_1_3_all plotFAO_1_3];

    plotdffNoAO_base_all = [plotdffNoAO_base_all plotdffNoAO_base];
    plotdffAO_base_all = [plotdffAO_base_all plotdffAO_base];
    plotdffNoAO_1_2_all = [plotdffNoAO_1_2_all plotdffNoAO_1_2];
    plotdffAO_1_2_all = [plotdffAO_1_2_all plotdffAO_1_2];
    plotdffNoAO_2_3_all = [plotdffNoAO_2_3_all plotdffNoAO_2_3];
    plotdffAO_2_3_all = [plotdffAO_2_3_all plotdffAO_2_3];
    plotdffNoAO_1_3_all = [plotdffNoAO_1_3_all plotdffNoAO_1_3];
    plotdffAO_1_3_all = [plotdffAO_1_3_all plotdffAO_1_3];
end

save([pfigsave 'results_all_FOV2data' num2str(folderInd) '_Validthre' num2str(Validthre) '_dffthre' num2str(dffthre) '.mat'],...
    'plotFNoAO_base_all','plotFAO_base_all','plotFNoAO_1_2_all','plotFAO_1_2_all',...
    'plotFNoAO_2_3_all','plotFAO_2_3_all','plotFNoAO_1_3_all','plotFAO_1_3_all',...
    'plotdffNoAO_base_all','plotdffAO_base_all','plotdffNoAO_1_2_all','plotdffAO_1_2_all',...
    'plotdffNoAO_2_3_all','plotdffAO_2_3_all','plotdffNoAO_1_3_all','plotdffAO_1_3_all','VEind','ValidInd',...
    'NormToNoAO','meanormax','Validthre','trialavgOnly','-v7.3');

%%
close all;
Oneplot = 0;
stdorSE = 1; % 1- std; 2- SE

if NormToNoAO==1
    ylimVal = [0 2];
    staLine = ylimVal(1)+(ylimVal(2)-ylimVal(1))*0.9; staLine = [staLine staLine];
    staStar = mean(staLine)+(ylimVal(2)-ylimVal(1))*0.03;
    yTickVal = ylimVal(1):0.5:ylimVal(2);
else
    ylimVal = [0 600]; % FOV2data8-11 combination
    staLine = ylimVal(1)+(ylimVal(2)-ylimVal(1))*0.9; staLine = [staLine staLine];
    staStar = mean(staLine)+(ylimVal(2)-ylimVal(1))*0.03;
    yTickVal = ylimVal(1):100:ylimVal(2);
end

% baseline
[F_NoAO_norm,F_AO_norm,h,p,x_positions] = plot_compare_v3(NormToNoAO,plotFNoAO_base_all,plotFAO_base_all,0,...
    Oneplot,1,stdorSE,ylimVal);
xlim([0.5,1.5]); ylabel('Signal AO/NoAO');
xticks([1]); xticklabels({''}); yticks(yTickVal);
title(['Pre stimulation, ' num2str(baselineLim(1)-1) ' - 0 s']); ylim(ylimVal)
statistic_plot([-0.2 0.2]+1,staLine,staStar,p);

if Displayratio==1
    x_dim = 0.8*size(IAO,2); y_dim = FOV_y*size(IAO,1);
    figure; colormap jet; hold on; xlim([0 size(IAO,2)]); ylim([0 size(IAO,1)]); set(gca,'YDir','reverse'); pbaspect([x_dim y_dim 1]);
    iGluROImap(bds,croprg1,croprg2,F_AO_norm./F_NoAO_norm);
    displayLim = [0 2];
    box on; set(gcf,'Color',[1,1,1]); set(gca,'FontSize',16); title(['Pre stimulation, ' num2str(baselineLim(1)-1) ' - 0 s, AO/NoAO']);
    clim(displayLim); h = colorbar; h.Position = [0.85 0.1095 0.02 0.2]; % For large FOV data
    h.Ticks = displayLim; set(gca,'Position',[0.1164    0.1100    0.6938    0.8150]);
    figure(1); iGluROImap(bds,croprg1,croprg2,F_AO_norm./F_NoAO_norm,'w');
    figure; colormap jet; hold on; xlim([0 size(IAO,2)]); ylim([0 size(IAO,1)]); set(gca,'YDir','reverse'); pbaspect([x_dim y_dim 1]);
    iGluROImap(bds,croprg1,croprg2,F_AO_norm./F_NoAO_norm,'w',1);
    box on; set(gcf,'Color',[1,1,1]); set(gca,'FontSize',16); title(['ROI map']);
    clim(displayLim); h = colorbar; h.Position = [0.85 0.1095 0.02 0.2]; % For large FOV data
    h.Ticks = displayLim; set(gca,'Position',[0.1164    0.1100    0.6938    0.8150]);
end

% post stimulation, 1-2s
[F_NoAO_norm,F_AO_norm,h,p] = plot_compare_v3(NormToNoAO,plotFNoAO_1_2_all,plotFAO_1_2_all,0,...
    Oneplot,2,stdorSE,ylimVal);
xlim([0.5,1.5]); ylabel('Signal AO/NoAO');
xticks([1]); xticklabels({''}); yticks(yTickVal);
title(['Post stimulation, 0 - 1 s']); ylim(ylimVal)
statistic_plot([-0.2 0.2]+1,staLine,staStar,p);

if Displayratio==1
    figure; colormap jet; hold on; xlim([0 size(IAO,2)]); ylim([0 size(IAO,1)]); set(gca,'YDir','reverse'); pbaspect([x_dim y_dim 1]);
    iGluROImap(bds,croprg1,croprg2,F_AO_norm./F_NoAO_norm);
    displayLim = [0 2];
    box on; set(gcf,'Color',[1,1,1]); set(gca,'FontSize',16); title(['Post stimulation, 0 - 1 s, AO/NoAO']);
    clim(displayLim); h = colorbar; h.Position = [0.85 0.1095 0.02 0.2]; % For large FOV data
    h.Ticks = displayLim; set(gca,'Position',[0.1164    0.1100    0.6938    0.8150]);
end

if Oneplot==1
    saveFigs_mod(pfigsave,1,'.fig',['Summary_F_combine_norm']);
    saveFigs_mod(pfigsave,1,'.tif',['Summary_F_combine_norm']);
    saveFigs_mod(pfigsave,1,'.pdf',['Summary_F_combine_norm']);
else
    saveFigs_mod(pfigsave,1,'.fig',['FOV2data' num2str(folderInd) '_Summary_F_norm_202601']);
    saveFigs_mod(pfigsave,1,'.tif',['FOV2data' num2str(folderInd) '_Summary_F_norm_202601']);
end

%%
close all;
if Oneplot==1
    savestr = ['Summary_dff_combine_dffthre' num2str(dffthre)];
else
    savestr = ['FOV2data' num2str(folderInd) '_Summary_Validthre' num2str(Validthre) '_dffthre' num2str(dffthre)];
end

if NormToNoAO==1
    ylimVal = [0 3];
    staLine = ylimVal(1)+(ylimVal(2)-ylimVal(1))*0.9; staLine = [staLine staLine];
    staStar = mean(staLine)+(ylimVal(2)-ylimVal(1))*0.03;
    yTickVal = ylimVal(1):1:ylimVal(2);
else
    ylimVal = [-0.1 0.4]; % FOV2data8-11 combination
    staLine = ylimVal(1)+(ylimVal(2)-ylimVal(1))*0.9; staLine = [staLine staLine];
    staStar = mean(staLine)+(ylimVal(2)-ylimVal(1))*0.03;
    yTickVal = ylimVal(1):0.1:ylimVal(2);
end

% post stimulation, 1-2s, [0-1s] post stimulation
Oneplotind = 2;
[dff_NoAO_g1,dff_AO_g1,ind_g1,dff_NoAO_g2,dff_AO_g2,ind_g2] = plot_compare_iGlu4_v3(plotdffNoAO_1_2_all,plotdffAO_1_2_all,...
    dffthre,Oneplot,Oneplotind,stdorSE,[0 3],[-0.05 0.15],['Post stimulation, 0 - 1 s'],[],[]);
figure(1); text(0.6,2,['Validthre ' num2str(Validthre)]);
text(0.6,1.8,['dffthre ' num2str(dffthre)]);
figure(2); text(0.6,0.1,['Validthre ' num2str(Validthre)]);
text(0.6,0.08,['dffthre ' num2str(dffthre)]);

save([pfigsave 'results_all_FOV2data' num2str(folderInd) '_Validthre' num2str(Validthre) '_dffthre' num2str(dffthre) '.mat'],...
    'ind_g1','ind_g2','-append');

saveFigs_mod(pfigsave,1,'.fig',savestr);
saveFigs_mod(pfigsave,1,'.tif',savestr);

%% Aligning the X axis
clear; clc;
close all;
load('H:\FACED_AO\20241217_M462_iGlu4\FACED\FOV2data8   9  10  11_VS_All\ROI_update\Compare_fig_202601_trialavg\results_all_FOV2data8   9  10  11_Validthre0.35_dffthre0.02.mat');
psave = 'H:\FACED_AO\20241217_M462_iGlu4\FACED\FOV2data8   9  10  11_VS_All\ROI_update\Compare_fig_202601_trialavg\';
pathroot = 'H:\FACED_AO\20241217_M462_iGlu4\FACED\';
Oneplot = 0;
stdorSE = 1; % 1- std; 2- SE
baselineLim = [0.5 1];
Validthre = 0.35; dffthre = 0.02;

if NormToNoAO==1
    ylimVal = [0 2];
    staLine = ylimVal(1)+(ylimVal(2)-ylimVal(1))*0.9; staLine = [staLine staLine];
    staStar = mean(staLine)+(ylimVal(2)-ylimVal(1))*0.03;
    yTickVal = ylimVal(1):0.5:ylimVal(2);
else
    ylimVal = [0 600]; % FOV2data8-11 combination
    staLine = ylimVal(1)+(ylimVal(2)-ylimVal(1))*0.9; staLine = [staLine staLine];
    staStar = mean(staLine)+(ylimVal(2)-ylimVal(1))*0.03;
    yTickVal = ylimVal(1):100:ylimVal(2);
end

% baseline
[F_NoAO_norm,F_AO_norm,h,p,x_positions] = plot_compare_v3(NormToNoAO,plotFNoAO_base_all,plotFAO_base_all,0,...
    Oneplot,1,stdorSE,ylimVal);
xlim([0.5,1.5]); ylabel('Signal AO/NoAO');
xticks([1]); xticklabels({''}); yticks(yTickVal);
title(['Pre stimulation, ' num2str(baselineLim(1)-1) ' - 0 s']); ylim(ylimVal)
statistic_plot([-0.2 0.2]+1,staLine,staStar,p);

% post stimulation, 1-2s, [0-1s] post stimulation
[F_NoAO_norm,F_AO_norm,h,p,x_positions] = plot_compare_v3(NormToNoAO,plotFNoAO_1_2_all,plotFAO_1_2_all,0,...
    Oneplot,2,stdorSE,ylimVal);
xlim([0.5,1.5]); ylabel('Signal AO/NoAO');
xticks([1]); xticklabels({''}); yticks(yTickVal);
title(['Post stimulation, 0 - 1 s']); ylim(ylimVal)
statistic_plot([-0.2 0.2]+1,staLine,staStar,p);
saveFigs_mod(psave,1,'.fig',['data_all_F']);
saveFigs_mod(psave,1,'.tif',['data_all_F']);

%%
close all;
if NormToNoAO==1
    ylimVal = [0 3];
    staLine = ylimVal(1)+(ylimVal(2)-ylimVal(1))*0.9; staLine = [staLine staLine];
    staStar = mean(staLine)+(ylimVal(2)-ylimVal(1))*0.03;
    yTickVal = ylimVal(1):1:ylimVal(2);
else
    ylimVal = [-0.1 0.4]; % FOV2data8-11 combination
    staLine = ylimVal(1)+(ylimVal(2)-ylimVal(1))*0.9; staLine = [staLine staLine];
    staStar = mean(staLine)+(ylimVal(2)-ylimVal(1))*0.03;
    yTickVal = ylimVal(1):0.1:ylimVal(2);
end
% post stimulation, 1-2s
Oneplotind = 2;
[dff_NoAO_g1,dff_AO_g1,ind_g1,dff_NoAO_g2,dff_AO_g2,ind_g2,x_positions_dff_g1,x_positions_dff_g2] = plot_compare_iGlu4_v3(plotdffNoAO_1_2_all,plotdffAO_1_2_all,...
    dffthre,Oneplot,Oneplotind,stdorSE,[0 3],[-0.05 0.15],['Post stimulation, 0 - 1 s'],[],[]);
figure(1); text(0.6,2,['Validthre ' num2str(Validthre)]);
text(0.6,1.8,['dffthre ' num2str(dffthre)]);
figure(2); text(0.6,0.1,['Validthre ' num2str(Validthre)]);
text(0.6,0.08,['dffthre ' num2str(dffthre)]);

saveFigs_mod(psave,1,'.fig',['data_all_dff']);
saveFigs_mod(psave,1,'.tif',['data_all_dff']);
close all;

%%
F_start_ind = 1;
dff_g1_start_ind = 1;
dff_g2_start_ind = 1;
for folderInd = [8:11]
    close all;
    pfigsave = [pathroot 'FOV2data' num2str([8:11]) '_VS_All\ROI_update\Compare_fig_202601_trialavg\']
    load([pfigsave 'results_all_FOV2data' num2str(folderInd) '_Validthre' num2str(Validthre) '_dffthre' num2str(dffthre) '.mat']);

    F_end_ind = F_start_ind+length(plotFNoAO_base_all)-1;
    if NormToNoAO==1
        ylimVal = [0 2];
        staLine = ylimVal(1)+(ylimVal(2)-ylimVal(1))*0.9; staLine = [staLine staLine];
        staStar = mean(staLine)+(ylimVal(2)-ylimVal(1))*0.03;
        yTickVal = ylimVal(1):0.5:ylimVal(2);
    else
        ylimVal = [0 600]; % FOV2data8-11 combination
        staLine = ylimVal(1)+(ylimVal(2)-ylimVal(1))*0.9; staLine = [staLine staLine];
        staStar = mean(staLine)+(ylimVal(2)-ylimVal(1))*0.03;
        yTickVal = ylimVal(1):100:ylimVal(2);
    end
    % baseline
    [F_NoAO_norm,F_AO_norm,h,p,~] = plot_compare_v3(NormToNoAO,plotFNoAO_base_all,plotFAO_base_all,0,...
        Oneplot,1,stdorSE,ylimVal,x_positions(F_start_ind:F_end_ind));
    xlim([0.5,1.5]); ylabel('Signal AO/NoAO');
    xticks([1]); xticklabels({''}); yticks(yTickVal);
    title(['Pre stimulation, ' num2str(baselineLim(1)-1) ' - 0 s']); ylim(ylimVal)
    statistic_plot([-0.2 0.2]+1,staLine,staStar,p);

    % post stimulation, 1-2s
    [F_NoAO_norm,F_AO_norm,h,p,~] = plot_compare_v3(NormToNoAO,plotFNoAO_1_2_all,plotFAO_1_2_all,0,...
        Oneplot,2,stdorSE,ylimVal,x_positions(F_start_ind:F_end_ind));
    xlim([0.5,1.5]); ylabel('Signal AO/NoAO');
    xticks([1]); xticklabels({''}); yticks(yTickVal);
    title(['Post stimulation, 0 - 1 s']); ylim(ylimVal)
    statistic_plot([-0.2 0.2]+1,staLine,staStar,p);
    saveFigs_mod(psave,1,'.fig',['data_' num2str(folderInd) '_F']);
    saveFigs_mod(psave,1,'.tif',['data_' num2str(folderInd) '_F']);

    F_start_ind = F_end_ind + 1;

    close all;
    dff_g1_end_ind = dff_g1_start_ind+length(ind_g1)-1;
    dff_g2_end_ind = dff_g2_start_ind+length(ind_g2)-1;
    if NormToNoAO==1
        ylimVal = [0 3];
        staLine = ylimVal(1)+(ylimVal(2)-ylimVal(1))*0.9; staLine = [staLine staLine];
        staStar = mean(staLine)+(ylimVal(2)-ylimVal(1))*0.03;
        yTickVal = ylimVal(1):1:ylimVal(2);
    else
        ylimVal = [-0.1 0.4]; % FOV2data8-11 combination
        staLine = ylimVal(1)+(ylimVal(2)-ylimVal(1))*0.9; staLine = [staLine staLine];
        staStar = mean(staLine)+(ylimVal(2)-ylimVal(1))*0.03;
        yTickVal = ylimVal(1):0.1:ylimVal(2);
    end
    % post stimulation, 1-2s
    Oneplotind = 2;
    [dff_NoAO_g1,dff_AO_g1,ind_g1,dff_NoAO_g2,dff_AO_g2,ind_g2,~,~] = plot_compare_iGlu4_v3(plotdffNoAO_1_2_all,plotdffAO_1_2_all,...
        dffthre,Oneplot,Oneplotind,stdorSE,[0 3],[-0.05 0.15],['Post stimulation, 0 - 1 s'],...
        x_positions_dff_g1(dff_g1_start_ind:dff_g1_end_ind),structfun(@(x) x(dff_g2_start_ind:dff_g2_end_ind), x_positions_dff_g2, 'UniformOutput', false));
    figure(1); text(0.6,2,['Validthre ' num2str(Validthre)]);
    text(0.6,1.8,['dffthre ' num2str(dffthre)]);
    figure(2); text(0.6,0.1,['Validthre ' num2str(Validthre)]);
    text(0.6,0.08,['dffthre ' num2str(dffthre)]);

    saveFigs_mod(psave,1,'.fig',['data_' num2str(folderInd) '_dff']);
    saveFigs_mod(psave,1,'.tif',['data_' num2str(folderInd) '_dff']);

    dff_g1_start_ind = dff_g1_end_ind+1;
    dff_g2_start_ind = dff_g2_end_ind+1;
end
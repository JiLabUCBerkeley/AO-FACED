%% Improvement quantification
% Compare AO and NoAO line-profile contrast across selected ROIs.
clear; clc; close all;
contrast_cal = 3;

pathfolder = 'C:\Users\JZhu\Box\SPA_JunZhu\JiLab\FACED\AO_FACED\Manuscript_AO_FACED\Structure_Tdtomato_20260327\M507_data1_202606_update\';
roiind = [1 2 3];
% pathfolder = 'C:\Users\JZhu\Box\SPA_JunZhu\JiLab\FACED\AO_FACED\Manuscript_AO_FACED\Structure_Tdtomato_20260327\FOV2data4_202606_update\';
% roiind = [2 11 3];

% cd(pathfolder);

% Initialize contrast, background, and peak-valley results.
SigRatioNobkgAll = [];
ratio = [];
cst_AO = [];
cst_NoAO = [];
bkg_AO_all = [];
bkg_NoAO_all = [];
maskall = cell(1,length(roiind));

cst_AO_Norm = [];
cst_NoAO_Norm = [];
peakValley_AO = [];
peakValley_NoAO = [];
peakValley_ratio = [];
peakValleyMetrics = cell(1,length(roiind));

% Load each ROI profile and calculate the selected contrast metrics.
iicount = 0;
for ROIi = roiind
    iicount = iicount+1;
    load([pathfolder 'ROI_' num2str(ROIi) '.mat']);

    close all;
    uiopen([pathfolder 'ROI_' num2str(ROIi) '_YZ_profile.fig'],1);

    lines = findall(gcf, 'Type', 'line');  % Find all line objects
    figure; plot(lines(1).XData,lines(1).YData); hold on;
    plot(lines(2).XData,lines(2).YData)
    input('Enter')

    Line_AO = lines(2).YData; Line_AO(Line_AO<0) = [];
    Line_NoAO = lines(1).YData; Line_NoAO(Line_NoAO<0) = [];

    ratio = [ratio max(Line_AO)/max(Line_NoAO)];

    % Optionally measure peak-valley contrast from the line profiles.
    respPeakValley = menu('Proceed with peak-valley contrast for this ROI?', ...
        'Yes', ...
        'No');
    if respPeakValley == 1
        if numel(lines(2).XData) == numel(Line_AO)
            LineX = lines(2).XData;
        else
            LineX = 1:numel(Line_AO);
        end
        [pvAO, pvNoAO, pvMetrics] = peak2_valley_contrast(Line_AO, Line_NoAO, LineX, ...
            ['ROI ' num2str(ROIi)]);
        peakValley_AO = [peakValley_AO pvAO];
        peakValley_NoAO = [peakValley_NoAO pvNoAO];
        peakValley_ratio = [peakValley_ratio pvAO/pvNoAO];
        peakValleyMetrics{iicount} = pvMetrics;
    else
        peakValley_AO = [peakValley_AO nan];
        peakValley_NoAO = [peakValley_NoAO nan];
        peakValley_ratio = [peakValley_ratio nan];
        peakValleyMetrics{iicount} = [];
    end

    % if ROIi==2 % only for GFP data2, ROI2
    [];
    % else
    if contrast_cal==1 % old contrast
        cst_AO = [cst_AO (max(Line_AO)-min(Line_AO))/(min(Line_AO))];
        cst_NoAO = [cst_NoAO (max(Line_NoAO)-min(Line_NoAO))/(min(Line_NoAO))];
    elseif contrast_cal==2 % Michelson contrast
        cst_AO = [cst_AO (max(Line_AO)-min(Line_AO))/(min(Line_AO)+max(Line_AO))];
        cst_NoAO = [cst_NoAO (max(Line_NoAO)-min(Line_NoAO))/(min(Line_NoAO)+max(Line_NoAO))];
    elseif contrast_cal==3 % Weber contrast & Normalized difference contrast
        % Select matched background regions in the AO and NoAO images.
        close all;

        uiopen([pathfolder 'ROI_' num2str(ROIi) '_sysAO_YZ.fig'],1);
        fig_a = gcf; ax_a = fig_a.CurrentAxes;
        imgsys = findobj(ax_a,'Type','image').CData;

        uiopen([pathfolder 'ROI_' num2str(ROIi) '_sampleAO_YZ.fig'],1);
        pause(1);
        fig_b = gcf; ax_b = fig_b.CurrentAxes;
        imgAO = findobj(ax_b,'Type','image').CData;

        % ---------- Ask whether to adjust colors first ----------
        resp0 = menu('Adjust colors before ROI selection?', 'Yes', 'No');
        if resp0 == 1
            manualAdjustCLim(ax_a, ax_b);
        end
        % -------------------------------------------------------

        while(1)

            % --- Draw ROI in figure(2) (reference mask) ---
            figure(2);
            roi = drawrectangle();
            mask = createMask(roi);
            B = bwboundaries(mask);

            h  = gobjects(length(B),1);   % boundaries in figure(2)
            h1 = gobjects(length(B),1);   % boundaries in figure(1)

            dx = 0; dy = 0;

            % --- Plot boundaries in figure(2) and figure(1) ---
            for k = 1:length(B)
                boundary = B{k};

                figure(2);
                h(k) = plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 2);

                figure(1);
                h1(k) = plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 2);
            end

            % --- Draggable proxy ROI in figure(1) ---
            figure(1);
            ax1 = gca;

            roi_proxy = drawrectangle(ax1, ...
                'Position', roi.Position, ...
                'Color', 'r', ...
                'LineWidth', 2);

            roi_proxy.InteractionsAllowed = "translate";

            pos0 = roi.Position;
            hasMoved = false;   % track first movement

            addlistener(roi_proxy,'MovingROI', @(src,evt) onMove(evt.CurrentPosition));
            addlistener(roi_proxy,'ROIMoved',  @(src,evt) onMove(evt.CurrentPosition));

            % --- Background calculations (fixed mask from figure 2) ---
            bkg_NoAO = imgsys.*mask;
            bkg_NoAO(bkg_NoAO==0) = nan;
            bkg_NoAO = nanmean(bkg_NoAO(:));
            bkg_NoAO_all = [bkg_NoAO_all, bkg_NoAO];

            bkg_AO = imgAO.*mask;
            bkg_AO(bkg_AO==0) = nan;
            bkg_AO = nanmean(bkg_AO(:));
            bkg_AO_all = [bkg_AO_all, bkg_AO];

            % --- User decision menu ---
            resp = menu('Keep current mask?', ...
                'Yes', ...
                'No', ...
                'Adjust colors (keep ROI)');

            if resp==1
                saveas(gcf,[pathfolder 'ROI_' num2str(ROIi) '_bkg.fig']);
                break;

            elseif resp==2
                delete(roi);
                delete(roi_proxy);
                delete(h(ishandle(h)));
                delete(h1(ishandle(h1)));
                continue;

            elseif resp==3
                manualAdjustCLim(ax_a, ax_b);
                continue;
            end

        end
    end

    %% Calculate background-corrected contrast
    % Store Weber, normalized, and signal-ratio measurements for this ROI.
    maskall{1,iicount} = mask;

    cst_AO = [cst_AO (max(Line_AO)-bkg_AO)/bkg_AO];
    cst_NoAO = [cst_NoAO (max(Line_NoAO)-bkg_NoAO)/bkg_NoAO];

    cst_AO_Norm = [cst_AO_Norm (max(Line_AO)-bkg_AO)/(max(Line_AO)+bkg_AO)];
    cst_NoAO_Norm = [cst_NoAO_Norm (max(Line_NoAO)-bkg_NoAO)/(max(Line_NoAO)+bkg_NoAO)];

    SigRatioNobkgAll = [SigRatioNobkgAll (max(Line_AO)-bkg_AO)/(max(Line_NoAO)-bkg_NoAO)];

    % end
end

%% Save contrast results
% Calculate AO/NoAO ratios and save all ROI measurements.
cst_ratio = cst_AO./cst_NoAO;

if contrast_cal==3
    cst_ratio_norm = cst_AO_Norm./cst_NoAO_Norm;
    save([pathfolder 'ratio_202605.mat'],'contrast_cal','ratio','cst_AO','cst_NoAO','cst_ratio',...
        'cst_AO_Norm','cst_NoAO_Norm','cst_ratio_norm','bkg_NoAO_all','bkg_AO_all','maskall','SigRatioNobkgAll',...
        'peakValley_AO','peakValley_NoAO','peakValley_ratio','peakValleyMetrics');
else
    save([pathfolder 'ratio.mat'],'ratio','cst_AO','cst_NoAO','cst_ratio',...
        'peakValley_AO','peakValley_NoAO','peakValley_ratio','peakValleyMetrics');
end

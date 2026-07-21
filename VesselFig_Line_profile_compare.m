%% Configure datasets and analysis parameters
% Set fixed PVD parameters and input locations for both datasets.
clear; clc; close all;

data2024_values = [0 1];
adaptiveFwhmThreshold_values = 0.1; % using std/mean
slidingWindowSec_values = 0.1;
adaptiveFwhmThreshold = adaptiveFwhmThreshold_values;
slidingWindowSec = slidingWindowSec_values;
pathall_data2024_1 = 'H:\FACED_AO\20241209_M420_Vessel\FACED\';
pathall_data2024_0 = 'K:\AO_FACED2026\20260416_RBead_vessel\FACED\';
allExampleMetricTable = [];

paramTag = ['thr' num2str(adaptiveFwhmThreshold,'%0.3g') '_win' num2str(slidingWindowSec,'%0.3g')];
paramExampleMetricTable = [];

%% Process both datasets
% Load each dataset, calculate ROI metrics, and combine the results.
for data2024 = data2024_values

    if data2024==1
        pathall = pathall_data2024_1; folderinds = [1,4,5,8,9,10];
        folderpattern = '*202505*';
    else
        pathall = pathall_data2024_0; folderinds = [4:11,13:14];
        folderpattern = '*202604*';
    end

    patchoutROI = []; % using plots of those ROIs from previous plots
    contrast_cal = 2; % 2
    Nratio = 0.15; % 0.2
    figuredisplay = 0;
    plotPVD = 0;

    Inten_AO_all = [];
    Inten_NoAO_all = [];
    Inten_RBC_AO_all = [];
    Inten_RBC_NoAO_all = [];
    Inten_AO_cst_all = [];
    Inten_NoAO_cst_all = [];
    Inten_RBC_AO_cst_all = [];
    Inten_RBC_NoAO_cst_all = [];
    Dist_AO_all = [];
    Dist_NoAO_all = [];
    PVDMethod_AO_all = strings(0,1);
    PVDMethod_NoAO_all = strings(0,1);
    namecat_inten_all = {};
    namecat_cst_all = {};
    folderind_inten_all = [];
    folderind_cst_all = [];
    foldername_inten_all = strings(0,1);
    foldername_cst_all = strings(0,1);

    % Process all selected fields of view in the current dataset.
    for folderind = folderinds
        if data2024==1
            folderLabel = ['FOV3data' num2str(folderind)];
            datafolder = [pathall folderLabel '_VS\'];
        else
            folderLabel = ['data' num2str(folderind)];
            datafolder = [pathall folderLabel '_VS\'];
        end
        foldername = dir([datafolder folderpattern]);
        if isempty(foldername)
            warning('No folder matched %s in %s.', folderpattern, datafolder);
            continue;
        end

        pathfolder = [datafolder foldername(1).name '\'];
        figfolder = [datafolder 'Fig\'];
        disp(pathfolder);

        namecat_inten = {};
        namecat_cst = {};
        folderind_inten = [];
        folderind_cst = [];
        foldername_inten = strings(0,1);
        foldername_cst = strings(0,1);
        Inten_AO = [];
        Inten_RBC_AO = [];
        Inten_NoAO = [];
        Inten_RBC_NoAO = [];
        Inten_AO_cst = [];
        Inten_NoAO_cst = [];
        Inten_RBC_AO_cst = [];
        Inten_RBC_NoAO_cst = [];
        Dist_AO = [];
        Dist_NoAO = [];
        PVDMethod_AO = strings(0,1);
        PVDMethod_NoAO = strings(0,1);

        % Find and sort the saved line-profile ROI files.
        roifiles = dir([pathfolder 'ROI_*.mat']);
        roiinds = [];
        line_roiinds = [];
        line_roisubinds = [];
        roifilenames = {};
        for fileind = 1:length(roifiles)
            token_roi = regexp(roifiles(fileind).name,'^ROI_(\d+)\.mat$','tokens','once');
            token_line = regexp(roifiles(fileind).name,'^ROI_(\d+)_(\d+)\.mat$','tokens','once');
            if ~isempty(token_roi)
                roiinds(end+1,1) = str2double(token_roi{1});
            elseif ~isempty(token_line)
                line_roiinds(end+1,1) = str2double(token_line{1});
                line_roisubinds(end+1,1) = str2double(token_line{2});
                roifilenames{end+1,1} = roifiles(fileind).name;
            end
        end

        roiinds = unique(roiinds);
        ROINum = length(roiinds);
        [~, sortind] = sortrows([line_roiinds line_roisubinds], [1 2]);
        line_roiinds = line_roiinds(sortind);
        line_roisubinds = line_roisubinds(sortind);
        roifilenames = roifilenames(sortind);
        ROIFileNum = length(roifilenames);

        % Calculate intensity and PVD metrics for each ROI profile.
        for roiindind = 1:ROIFileNum
            roiind = line_roiinds(roiindind);
            if ismember(roiind, patchoutROI)
                continue;
            end
            filename = roifilenames{roiindind};
            [~, filenamebase, ~] = fileparts(filename);
            updatefilename = [filenamebase '_update.mat'];
            updatefile = [figfolder updatefilename];
            if exist(figfolder,'dir')==7 && exist(updatefile,'file')==2
                loadfile = updatefile;
                isupdatefile = 1;
            else
                loadfile = [pathfolder filename];
                isupdatefile = 0;
            end
            disp(loadfile);
            clear ROI_AO_line ROI_NoAO_line labelN tt_XZ;
            load(loadfile);
            if exist('labelN','var')~=1
                labelN = nan;
            end
            if isupdatefile==1
                labelN = 2;
            end

            [~,~,Inten_AO_curr,Inten_NoAO_curr,Inten_RBC_AO_curr,Inten_RBC_NoAO_curr] = VesselLineCst_v2...
                (ROI_AO_line,ROI_NoAO_line,contrast_cal,Nratio);
            Inten_AO = cat(1,Inten_AO,Inten_AO_curr);
            Inten_NoAO = cat(1,Inten_NoAO,Inten_NoAO_curr);
            Inten_RBC_AO = cat(1,Inten_RBC_AO,Inten_RBC_AO_curr);
            Inten_RBC_NoAO = cat(1,Inten_RBC_NoAO,Inten_RBC_NoAO_curr);
            namecat_inten = cat(1,namecat_inten,{[foldername(1).name filesep filename]});
            folderind_inten = cat(1,folderind_inten,folderind);
            foldername_inten = cat(1,foldername_inten,string(folderLabel));

            % PVD comparison
            if labelN==2
                namecat_cst = cat(1,namecat_cst,{[foldername(1).name filesep filename]});
                folderind_cst = cat(1,folderind_cst,folderind);
                foldername_cst = cat(1,foldername_cst,string(folderLabel));
                Inten_AO_cst = cat(1,Inten_AO_cst,Inten_AO_curr);
                Inten_NoAO_cst = cat(1,Inten_NoAO_cst,Inten_NoAO_curr);
                Inten_RBC_AO_cst = cat(1,Inten_RBC_AO_cst,Inten_RBC_AO_curr);
                Inten_RBC_NoAO_cst = cat(1,Inten_RBC_NoAO_cst,Inten_RBC_NoAO_curr);
                if exist('tt_XZ','var')==1 && numel(tt_XZ)>1
                    dT_curr = mean(diff(tt_XZ));
                else
                    dT_curr = [];
                end
                [PVD_AO_curr,~,PVDMethod_AO_curr] = localPeakValleyDistinguishability_v2( ...
                    ROI_AO_line,Nratio,'adaptive',dT_curr,slidingWindowSec,adaptiveFwhmThreshold);
                if PVDMethod_AO_curr=="sliding"
                    NoAOMethod_curr = 'slidingaverage';
                else
                    NoAOMethod_curr = char(PVDMethod_AO_curr);
                end
                [PVD_NoAO_curr,~,PVDMethod_NoAO_curr] = localPeakValleyDistinguishability_v2( ...
                    ROI_NoAO_line,Nratio,NoAOMethod_curr,dT_curr,slidingWindowSec,adaptiveFwhmThreshold);
                PVD_ratio_curr = PVD_AO_curr/PVD_NoAO_curr;
                Dist_AO = cat(1,Dist_AO,PVD_AO_curr);
                Dist_NoAO = cat(1,Dist_NoAO,PVD_NoAO_curr);
                PVDMethod_AO = cat(1,PVDMethod_AO,PVDMethod_AO_curr);
                PVDMethod_NoAO = cat(1,PVDMethod_NoAO,PVDMethod_NoAO_curr);

                if plotPVD==1
                    x_AO_PVD = (1:numel(ROI_AO_line))';
                    x_NoAO_PVD = (1:numel(ROI_NoAO_line))';
                    xLabel_PVD = 'Position';
                    if exist('tt_XZ','var')==1 && numel(tt_XZ)==numel(ROI_AO_line) && numel(tt_XZ)==numel(ROI_NoAO_line)
                        x_AO_PVD = tt_XZ(:);
                        x_NoAO_PVD = tt_XZ(:);
                        xLabel_PVD = 'Time';
                    end
                    close all;
                    figure; set(gcf,'Position',[100 100 900 520]);
                    plot(x_AO_PVD,ROI_AO_line,'r','LineWidth',1.2); hold on;
                    plot(x_NoAO_PVD,ROI_NoAO_line,'b','LineWidth',1.2);
                    legend('sample AO','sys AO','NumColumns',2,'Location','best');
                    xlabel(xLabel_PVD);
                    ylabel('Intensity');
                    set(gcf,'Position',[100 390 560 260]);
                    set(gca,'Position',[0.16 0.25 0.75 0.6]);
                    styleSavedProfileAxes(gca);
                    title([filenamebase ', PVD ' num2str(PVD_AO_curr,'%0.1f') ', ' num2str(PVD_NoAO_curr,'%0.1f') ', ' num2str(PVD_ratio_curr,'%0.1f') '; '...
                        'PVD method ' char(PVDMethod_AO_curr)],...
                        'Interpreter','none','FontSize',10);
                    saveFigs_mod(pathfolder,1,'.fig',[filenamebase '_PVD_' paramTag]);
                    saveFigs_mod(pathfolder,1,'.tif',[filenamebase '_PVD_' paramTag]);
                end

                if figuredisplay==1
                    close all;
                    figure; set(gcf,'Position',[100 100 720 700]);
                    win = linspace(min([ROI_AO_line;ROI_NoAO_line]),max([ROI_AO_line;ROI_NoAO_line]),50);
                    subplot(2,1,1); plot(ROI_AO_line,'r'); hold on; plot(ROI_NoAO_line,'b');
                    legend('sample AO','sys AO','NumColumns',2);
                    subplot(2,1,2); histogram(ROI_AO_line,win,'FaceColor','r'); hold on;
                    histogram(ROI_NoAO_line,win,'FaceColor','b');
                    for Nratioiii = 0.1:0.1:0.4
                        IntenAO_curr = mean(maxk(ROI_AO_line,round(length(ROI_AO_line)*Nratioiii)));
                        line([IntenAO_curr,IntenAO_curr],[0 100],'Color','r');
                        IntenAO_curr = mean(mink(ROI_AO_line,round(length(ROI_AO_line)*Nratioiii)));
                        line([IntenAO_curr,IntenAO_curr],[0 100],'Color','r','LineStyle','--');
                        IntenNoAO_curr = mean(maxk(ROI_NoAO_line,round(length(ROI_NoAO_line)*Nratioiii)));
                        line([IntenNoAO_curr,IntenNoAO_curr],[0 100],'Color','b');
                        IntenNoAO_curr = mean(mink(ROI_NoAO_line,round(length(ROI_NoAO_line)*Nratioiii)));
                        line([IntenNoAO_curr,IntenNoAO_curr],[0 100],'Color','b','LineStyle','--');
                    end
                    [~, savename] = fileparts(filename);
                    saveFigs_mod(figuresaveloc,1,'.fig',savename);
                    saveFigs_mod(figuresaveloc,1,'.tif',savename);
                end
            end
        end

        % Save metrics for the current field of view.
        save([pathfolder 'summary_Cst_cal_' num2str(contrast_cal) '_PVDupdate.mat'],...
            'contrast_cal','Inten_AO','Inten_NoAO','Inten_RBC_AO','Inten_RBC_NoAO',...
            'Dist_AO','Dist_NoAO','ROINum','ROIFileNum','namecat_inten','namecat_cst');

        Inten_AO_all = cat(1,Inten_AO_all,Inten_AO);
        Inten_NoAO_all = cat(1,Inten_NoAO_all,Inten_NoAO);
        Inten_RBC_AO_all = cat(1,Inten_RBC_AO_all,Inten_RBC_AO);
        Inten_RBC_NoAO_all = cat(1,Inten_RBC_NoAO_all,Inten_RBC_NoAO);
        Inten_AO_cst_all = cat(1,Inten_AO_cst_all,Inten_AO_cst);
        Inten_NoAO_cst_all = cat(1,Inten_NoAO_cst_all,Inten_NoAO_cst);
        Inten_RBC_AO_cst_all = cat(1,Inten_RBC_AO_cst_all,Inten_RBC_AO_cst);
        Inten_RBC_NoAO_cst_all = cat(1,Inten_RBC_NoAO_cst_all,Inten_RBC_NoAO_cst);
        Dist_AO_all = cat(1,Dist_AO_all,Dist_AO);
        Dist_NoAO_all = cat(1,Dist_NoAO_all,Dist_NoAO);
        PVDMethod_AO_all = cat(1,PVDMethod_AO_all,PVDMethod_AO);
        PVDMethod_NoAO_all = cat(1,PVDMethod_NoAO_all,PVDMethod_NoAO);
        namecat_inten_all = cat(1,namecat_inten_all,namecat_inten);
        namecat_cst_all = cat(1,namecat_cst_all,namecat_cst);
        folderind_inten_all = cat(1,folderind_inten_all,folderind_inten);
        folderind_cst_all = cat(1,folderind_cst_all,folderind_cst);
        foldername_inten_all = cat(1,foldername_inten_all,foldername_inten);
        foldername_cst_all = cat(1,foldername_cst_all,foldername_cst);
    end

    pfigS = [pathall 'PopulationTest_std_Nratio_' num2str(Nratio) '\Population20260601_' paramTag '\'];
    if exist(pfigS,'dir')~=7
        mkdir(pfigS);
    end

    save([pfigS 'summary_Cst_cal_' num2str(contrast_cal) '_allFolders_PVDupdate_' paramTag '.mat'],...
        'contrast_cal','Nratio','slidingWindowSec','adaptiveFwhmThreshold','paramTag','folderinds','folderpattern','patchoutROI',...
        'Inten_AO_all','Inten_NoAO_all','Inten_RBC_AO_all','Inten_RBC_NoAO_all',...
        'Inten_AO_cst_all','Inten_NoAO_cst_all','Inten_RBC_AO_cst_all','Inten_RBC_NoAO_cst_all',...
        'Dist_AO_all','Dist_NoAO_all','PVDMethod_AO_all','PVDMethod_NoAO_all',...
        'namecat_inten_all','namecat_cst_all','folderind_inten_all','folderind_cst_all',...
        'foldername_inten_all','foldername_cst_all');

    %% Plot population-level comparisons
    % Remove extreme ratios and compare AO with NoAO across all fields.
    NormToNoAO = 1;
    close all;
    Oneplot = 0;
    stdorSE = 1; % 1- std; 2- SE
    threshold_valid = 10;
    intenValid = ~(Inten_AO_all./Inten_NoAO_all > threshold_valid);
    distValid = ~(Dist_AO_all./Dist_NoAO_all > threshold_valid);
    Inten_AO_valid = Inten_AO_all(intenValid);
    Inten_NoAO_valid = Inten_NoAO_all(intenValid);
    Dist_AO_valid = Dist_AO_all(distValid);
    Dist_NoAO_valid = Dist_NoAO_all(distValid);

    if NormToNoAO==1
        ylimVal = [0 4];
        staLine = ylimVal(1)+(ylimVal(2)-ylimVal(1))*0.9; staLine = [staLine staLine];
        staStar = mean(staLine)+(ylimVal(2)-ylimVal(1))*0.03;
        yTickVal = ylimVal(1):1:ylimVal(2);
    else
        ylimVal = [0 600];
        staLine = ylimVal(1)+(ylimVal(2)-ylimVal(1))*0.9; staLine = [staLine staLine];
        staStar = mean(staLine)+(ylimVal(2)-ylimVal(1))*0.03;
        yTickVal = ylimVal(1):100:ylimVal(2);
    end

    [~,~,~,p,x_positions_Inten] = plot_compare_v3(NormToNoAO,Inten_NoAO_valid,Inten_AO_valid,0,...
        Oneplot,1,stdorSE,ylimVal);
    xlim([0.5,1.5]); ylabel('Plasma signal AO/NoAO');
    xticks([1]); xticklabels({''}); yticks(yTickVal);
    title('Plasma signal ratio'); ylim(ylimVal)
    statistic_plot([-0.2 0.2]+1,staLine,staStar,p);

    [~,~,~,p,x_positions_Dist] = plot_compare_v3(NormToNoAO,Dist_NoAO_valid,Dist_AO_valid,0,...
        Oneplot,1,stdorSE,ylimVal);
    xlim([0.5,1.5]); ylabel('Peak-valley distinguishability AO/NoAO');
    xticks([1]); xticklabels({''}); yticks(yTickVal);
    title('Peak-valley distinguishability ratio'); ylim(ylimVal)
    statistic_plot([-0.2 0.2]+1,staLine,staStar,p);

    saveFigs_mod(pfigS,1,'.fig',['All folder_XaxisMatch_' paramTag]);
    saveFigs_mod(pfigS,1,'.tif',['All folder_XaxisMatch_' paramTag]);
    saveFigs_mod(pfigS,1,'.pdf',['All folder_XaxisMatch_' paramTag]);

    %% Plot individual fields on the population axis
    % Reuse population x positions so all field-of-view plots are comparable.
    for folderind = folderinds
        close all;

        if data2024==1
            folderLabel = ['FOV3data' num2str(folderind)];
            datafolder = [pathall folderLabel '_VS\'];
        else
            folderLabel = ['data' num2str(folderind)];
            datafolder = [pathall folderLabel '_VS\'];
        end
        foldername = dir([datafolder folderpattern]);
        if isempty(foldername)
            warning('No folder matched %s in %s.', folderpattern, datafolder);
            continue;
        end

        pathfolder = [datafolder foldername(1).name '\'];
        disp(pathfolder);

        intenInd = find(folderind_inten_all(intenValid)==folderind);
        distInd = find(folderind_cst_all(distValid)==folderind);
        if isempty(intenInd) || isempty(distInd)
            warning('No current-pass metrics found for folderind %d.', folderind);
            continue;
        end

        [~,~,~,p,~] = plot_compare_v3(NormToNoAO,Inten_NoAO_valid(intenInd),Inten_AO_valid(intenInd),0,...
            Oneplot,1,stdorSE,ylimVal,x_positions_Inten(intenInd));
        xlim([0.5,1.5]); ylabel('Plasma signal AO/NoAO');
        xticks([1]); xticklabels({''}); yticks(yTickVal);
        title('Plasma signal ratio'); ylim(ylimVal)
        statistic_plot([-0.2 0.2]+1,staLine,staStar,p);

        [~,~,~,p,~] = plot_compare_v3(NormToNoAO,Dist_NoAO_valid(distInd),Dist_AO_valid(distInd),0,...
            Oneplot,1,stdorSE,ylimVal,x_positions_Dist(distInd));
        xlim([0.5,1.5]); ylabel('Peak-valley distinguishability AO/NoAO');
        xticks([1]); xticklabels({''}); yticks(yTickVal);
        title('Peak-valley distinguishability ratio'); ylim(ylimVal)
        statistic_plot([-0.2 0.2]+1,staLine,staStar,p);

        saveFigs_mod(pfigS,1,'.fig',['data_' num2str(folderind) '_XaxisMatch']);
        saveFigs_mod(pfigS,1,'.tif',['data_' num2str(folderind) '_XaxisMatch']);
    end

    %% Collect selected examples
    % Extract manuscript examples and store their intensity and PVD metrics.
    if data2024==0
        exampleGroup = ["3"; "4"; "5"; "6"; "7"; "8"];
        exampleData = [14; 14; 4; 9; 11; 11];
        exampleROI = ["ROI_18_2"; "ROI_17_2"; "ROI_7_3"; "ROI_14_1"; "ROI_7_1"; "ROI_7_2"];
    elseif data2024==1
        exampleGroup = ["1"; "2"];
        exampleData = [4; 4];
        exampleROI = ["ROI_51_1";"ROI_17_1"];
    else
        exampleGroup = strings(0,1);
        exampleData = [];
        exampleROI = strings(0,1);
    end

    if ~isempty(exampleROI)
        exampleFile = exampleROI+".mat";
        allNames = string(namecat_cst_all);
        exampleData2024 = repmat(data2024,numel(exampleFile),1);
        exampleThreshold = repmat(adaptiveFwhmThreshold,numel(exampleFile),1);
        exampleWindowSec = repmat(slidingWindowSec,numel(exampleFile),1);
        exampleParamTag = repmat(string(paramTag),numel(exampleFile),1);

        matched = false(numel(exampleFile),1);
        matchedName = strings(numel(exampleFile),1);
        plasmaRatio = nan(numel(exampleFile),1);
        rbcRatio = nan(numel(exampleFile),1);
        pvdRatio = nan(numel(exampleFile),1);
        pvdMethodAO = strings(numel(exampleFile),1);
        pvdMethodNoAO = strings(numel(exampleFile),1);
        plasmaAO = nan(numel(exampleFile),1);
        plasmaNoAO = nan(numel(exampleFile),1);
        rbcAO = nan(numel(exampleFile),1);
        rbcNoAO = nan(numel(exampleFile),1);
        pvdAO = nan(numel(exampleFile),1);
        pvdNoAO = nan(numel(exampleFile),1);

        for exampleInd = 1:numel(exampleFile)
            currInd = find(folderind_cst_all==exampleData(exampleInd) & ...
                endsWith(allNames,exampleFile(exampleInd)),1);
            if isempty(currInd)
                continue;
            end

            matched(exampleInd) = true;
            matchedName(exampleInd) = allNames(currInd);
            plasmaAO(exampleInd) = Inten_AO_cst_all(currInd);
            plasmaNoAO(exampleInd) = Inten_NoAO_cst_all(currInd);
            rbcAO(exampleInd) = Inten_RBC_AO_cst_all(currInd);
            rbcNoAO(exampleInd) = Inten_RBC_NoAO_cst_all(currInd);
            pvdAO(exampleInd) = Dist_AO_all(currInd);
            pvdNoAO(exampleInd) = Dist_NoAO_all(currInd);
            pvdMethodAO(exampleInd) = PVDMethod_AO_all(currInd);
            pvdMethodNoAO(exampleInd) = PVDMethod_NoAO_all(currInd);

            plasmaRatio(exampleInd) = plasmaAO(exampleInd)/plasmaNoAO(exampleInd);
            rbcRatio(exampleInd) = rbcAO(exampleInd)/rbcNoAO(exampleInd);
            pvdRatio(exampleInd) = pvdAO(exampleInd)/pvdNoAO(exampleInd);
        end

        exampleMetricTable = table(exampleData2024,exampleThreshold,exampleWindowSec,exampleParamTag, ...
            exampleGroup,exampleData,exampleROI, ...
            matched,matchedName,plasmaRatio,rbcRatio,pvdRatio,pvdMethodAO,pvdMethodNoAO, ...
            'VariableNames',{'Data2024','Threshold','WindowSec','ParamTag','Group','Data','ROI','Matched','MatchedName', ...
            'PlasmaRatio','RBCRatio','PVDRatio','PVDMethod_AO','PVDMethod_NoAO'});

        if isempty(allExampleMetricTable)
            allExampleMetricTable = exampleMetricTable;
        else
            allExampleMetricTable = cat(1,allExampleMetricTable,exampleMetricTable);
        end
        if isempty(paramExampleMetricTable)
            paramExampleMetricTable = exampleMetricTable;
        else
            paramExampleMetricTable = cat(1,paramExampleMetricTable,exampleMetricTable);
        end
    end

end

%% Save selected-example tables
% Combine the two datasets and export the selected metrics to MAT and CSV.
if ~isempty(paramExampleMetricTable)
    paramExampleMetricTable = paramExampleMetricTable([find(paramExampleMetricTable.Data2024==1); ...
        find(paramExampleMetricTable.Data2024~=1)],:);
    for tablePath = string({pathall_data2024_1,pathall_data2024_0})
        tableSaveFolder = pfigS; %fullfile(char(tablePath),['Population20260515_' paramTag]);
        if exist(tableSaveFolder,'dir')~=7
            mkdir(tableSaveFolder);
        end
        save(fullfile(tableSaveFolder,['SelectedMetricTable_' paramTag '.mat']),'paramExampleMetricTable');
        writetable(paramExampleMetricTable,fullfile(tableSaveFolder,['SelectedMetricTable_' paramTag '.csv']));
    end
end
%
if ~isempty(allExampleMetricTable)
    allExampleMetricTable = allExampleMetricTable([find(allExampleMetricTable.Data2024==1); ...
        find(allExampleMetricTable.Data2024~=1)],:);
    save('SelectedMetricTable_allParams.mat','allExampleMetricTable');
    writetable(allExampleMetricTable,'SelectedMetricTable_allParams.csv');
end

disp('Selected metric comparison examples, data2024==0 and data2024==1:');
disp(allExampleMetricTable);

%% Format the final active figure
% Apply the standard size, axes position, limits, and visual style.
set(gcf,'Position',[100 100 450 420]);
set(gca,'Position',[0.1300 0.1100 0.7750 0.8100]);
xlim([0.7 1.3]);
styleSavedProfileAxes(gca);
% yticks([0:0.5:2.5])

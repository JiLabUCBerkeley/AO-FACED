%% Configure flow analysis
% Set datasets, averaging, shift candidates, and the optional ROI limit.
clear; clc; close all;

data2024_values = [0 1];
NormToNoAO = 0;
pathall_data2024_1 = 'H:\FACED_AO\20241209_M420_Vessel\FACED\';
pathall_data2024_0 = 'K:\AO_FACED2026\20260416_RBead_vessel\FACED\';
numavgs = 50;
shiftamtCandidates = 1:5;
maxRoisToProcess = str2double(getenv('FLOW_QUANT_MAX_ROIS'));
if isnan(maxRoisToProcess) || maxRoisToProcess <= 0
    maxRoisToProcess = inf;
end
processedRois = 0;

allQcSummary = struct('data2024', {}, 'folderind', {}, 'filenamebase', {}, ...
    'kymoName', {}, 'fitR2_median', {}, 'cvVelocity', {}, 'shiftamt', {});

%% Quantify flow in both datasets
% Process accepted AO/NoAO kymographs and collect quality metrics.
for data2024 = data2024_values
    if processedRois >= maxRoisToProcess
        break;
    end

    if data2024 == 1
        pathall = pathall_data2024_1;
        folderinds = [1, 4, 5, 8, 9, 10];
        folderpattern = '*202505*';
        folderPrefix = 'FOV3data';
    else
        pathall = pathall_data2024_0;
        folderinds = [4:11, 13:14];
        folderpattern = '*202604*';
        folderPrefix = 'data';
    end

    qcSummary = struct('data2024', {}, 'folderind', {}, 'filenamebase', {}, ...
        'kymoName', {}, 'fitR2_median', {}, 'cvVelocity', {}, 'shiftamt', {});
    outputDir = fullfile(pathall, sprintf('FlowResults_v4_avg%d_adaptiveShift_data2024_%d_20260528', ...
        numavgs, data2024));
    if exist(outputDir, 'dir') ~= 7
        mkdir(outputDir);
    end

    % Find the saved line-profile ROIs in each field of view.
    for folderind = folderinds
        if processedRois >= maxRoisToProcess
            break;
        end

        datafolder = fullfile(pathall, [folderPrefix num2str(folderind) '_VS']);
        foldername = dir(fullfile(datafolder, folderpattern));
        if isempty(foldername)
            warning('No folder matched %s in %s.', folderpattern, datafolder);
            continue;
        end

        pathfolder = fullfile(datafolder, foldername(1).name);
        figfolder = fullfile(datafolder, 'Fig');
        roifiles = dir(fullfile(pathfolder, 'ROI_*.mat'));
        disp(pathfolder);

        roiinds = [];
        roisubinds = [];
        roifilenames = {};
        for fileind = 1:numel(roifiles)
            token = regexp(roifiles(fileind).name, '^ROI_(\d+)_(\d+)\.mat$', 'tokens', 'once');
            if ~isempty(token)
                roiinds(end + 1, 1) = str2double(token{1});
                roisubinds(end + 1, 1) = str2double(token{2});
                roifilenames{end + 1, 1} = roifiles(fileind).name;
            end
        end
        [~, sortind] = sortrows([roiinds roisubinds], [1 2]);
        roiinds = roiinds(sortind);
        roisubinds = roisubinds(sortind);
        roifilenames = roifilenames(sortind);

        % Load, validate, and quantify each accepted ROI.
        for roiindind = 1:numel(roifilenames)
            if processedRois >= maxRoisToProcess
                break;
            end

            ROIind = roiinds(roiindind);
            plotind = roisubinds(roiindind);
            [~, filenamebase] = fileparts(roifilenames{roiindind});
            roiInfoFile = fullfile(pathfolder, sprintf('ROI_%d.mat', ROIind));
            lineFile = fullfile(pathfolder, roifilenames{roiindind});
            updateFile = fullfile(figfolder, [filenamebase '_update.mat']);

            lineVars = load(lineFile);
            if exist(updateFile, 'file') == 2
                updateVars = load(updateFile);
                updateFields = fieldnames(updateVars);
                for fieldInd = 1:numel(updateFields)
                    lineVars.(updateFields{fieldInd}) = updateVars.(updateFields{fieldInd});
                end
                labelN = 2;
                loadfile = updateFile;
                isUpdateFile = true;
            else
                if isfield(lineVars, 'labelN')
                    labelN = lineVars.labelN;
                else
                    labelN = nan;
                end
                loadfile = lineFile;
                isUpdateFile = false;
            end

            if labelN ~= 2
                continue;
            end
            if ~isfield(lineVars, 'ROI_AO') || ~isfield(lineVars, 'ROI_NoAO')
                warning('Skipping %s: ROI_AO or ROI_NoAO missing.', loadfile);
                continue;
            end
            if ~isfield(lineVars, 'tt_XZ') || ~isfield(lineVars, 'up_factor')
                warning('Skipping %s: tt_XZ or up_factor missing.', loadfile);
                continue;
            end

            % Select the best AO shift and apply it to both conditions.
            roiInfo = load(roiInfoFile, 'dXY_FACED');
            frmRate = 1 / mean(diff(lineVars.tt_XZ)) / lineVars.up_factor;
            [shiftamt, shiftamtScores] = bestShiftamtForAO(lineVars.ROI_AO, frmRate, ...
                roiInfo.dXY_FACED, numavgs, shiftamtCandidates, outputDir, ...
                data2024, folderind, filenamebase);

            kymoData = struct('AO', lineVars.ROI_AO, 'NoAO', lineVars.ROI_NoAO);
            kymoNames = fieldnames(kymoData);
            for kymoInd = 1:numel(kymoNames)
                kymoName = kymoNames{kymoInd};
                imageLines = kymoData.(kymoName);
                params = flowParams(numavgs, shiftamt, frmRate, size(imageLines, 2), ...
                    outputDir, data2024, folderind, filenamebase, kymoName, true);
                resultFile = fullfile(outputDir, [params.savePrefix '_flow_quant.mat']);

                if exist(resultFile, 'file') == 2
                    load(resultFile, 'qc');
                    fprintf('\nUsing saved result: %s\n', resultFile);
                else
                    fprintf('\nQuantifying data2024=%d data%d %s %s, avg%d shift%d...\n', ...
                        data2024, folderind, filenamebase, kymoName, numavgs, shiftamt);
                    result = flow_quant_func(imageLines, frmRate, roiInfo.dXY_FACED, params);
                    qc = flowQc(result);
                    qc.numavgs = numavgs;
                    qc.shiftamt = shiftamt;

                    save(resultFile, 'result', 'params', 'qc', 'shiftamtScores', ...
                        'shiftamtCandidates', 'pathfolder', 'loadfile', 'lineFile', ...
                        'updateFile', 'isUpdateFile', 'roiInfoFile', 'data2024', ...
                        'folderind', 'ROIind', 'plotind', 'filenamebase', 'kymoName');
                end

                qcSummary(end + 1).data2024 = data2024;
                qcSummary(end).folderind = folderind;
                qcSummary(end).filenamebase = filenamebase;
                qcSummary(end).kymoName = kymoName;
                qcSummary(end).fitR2_median = qc.fitR2_median;
                qcSummary(end).cvVelocity = qc.cvVelocity;
                qcSummary(end).shiftamt = shiftamt;
                fprintf('%s %s median R2: %.3f, velocity coefficient of variation: %.3f, shiftamt: %d\n', ...
                    filenamebase, kymoName, qc.fitR2_median, qc.cvVelocity, shiftamt);
                close all;
            end

            processedRois = processedRois + 1;
        end
    end

    %% Save and compare quality-control metrics
    % Save dataset results and plot AO versus NoAO fit quality.
    allQcSummary = cat(2, allQcSummary, qcSummary);
    save(fullfile(outputDir, 'flow_qc_summary_v4.mat'), ...
        'qcSummary', 'NormToNoAO', 'numavgs', 'shiftamtCandidates', ...
        'data2024', 'folderinds', 'folderpattern');

    roiKeys = strings(numel(qcSummary), 1);
    for i = 1:numel(qcSummary)
        roiKeys(i) = sprintf('%d_data%d_%s', qcSummary(i).data2024, ...
            qcSummary(i).folderind, qcSummary(i).filenamebase);
    end
    roiNames = unique(roiKeys, 'stable');
    r2NoAO = []; r2AO = []; cvNoAO = []; cvAO = []; labels = {};
    for i = 1:numel(roiNames)
        parts = regexp(char(roiNames(i)), '^(\d+)_data(\d+)_(.*)$', 'tokens', 'once');
        isCurr = [qcSummary.data2024] == str2double(parts{1}) & ...
            [qcSummary.folderind] == str2double(parts{2}) & ...
            strcmp({qcSummary.filenamebase}, parts{3});
        aoInd = find(isCurr & strcmp({qcSummary.kymoName}, 'AO'), 1);
        noaoInd = find(isCurr & strcmp({qcSummary.kymoName}, 'NoAO'), 1);
        if isempty(aoInd) || isempty(noaoInd)
            continue;
        end
        r2NoAO(end + 1, 1) = qcSummary(noaoInd).fitR2_median;
        r2AO(end + 1, 1) = qcSummary(aoInd).fitR2_median;
        cvNoAO(end + 1, 1) = qcSummary(noaoInd).cvVelocity;
        cvAO(end + 1, 1) = qcSummary(aoInd).cvVelocity;
        labels{end + 1, 1} = sprintf('d%s data%s %s', parts{1}, parts{2}, parts{3});
    end

    if ~isempty(labels)
        Oneplot = 0;
        stdorSE = 1;
        titleText = sprintf('data2024 = %d', data2024);
        tag = regexprep(titleText, '\W+', '_');

        ylimVal = [0 1];
        [~, ~, ~, ~, xPositionsR2] = plot_compare_v3(NormToNoAO, r2NoAO, r2AO, 0, ...
            Oneplot, 1, stdorSE, ylimVal);
        xlim([0.5 2.5]); xticks([1 2]); xticklabels({'NoAO', 'AO'});
        ylabel('Median fit R^2');
        title(['Median fit R^2, ' titleText], 'Interpreter', 'none');
        ylim(ylimVal);
        savefig(gcf, fullfile(outputDir, ['QC_R2_comparison_' tag '.fig']));
        print(gcf, '-r600', '-dtiff', fullfile(outputDir, ['QC_R2_comparison_' tag '.tif']));

        cvMax = max([cvNoAO(:); cvAO(:)], [], 'omitnan');
        if ~isfinite(cvMax) || cvMax <= 0
            cvMax = 1;
        end
        ylimVal = [0 1.1 * cvMax];
        plot_compare_v3(NormToNoAO, cvNoAO, cvAO, 0, ...
            Oneplot, 1, stdorSE, ylimVal);
        xlim([0.5 2.5]); xticks([1 2]); xticklabels({'NoAO', 'AO'});
        ylabel('Velocity coefficient of variation');
        title(['Velocity coefficient of variation, ' titleText], 'Interpreter', 'none');
        ylim(ylimVal);
        savefig(gcf, fullfile(outputDir, ['QC_velocityCoefficientOfVariation_comparison_' tag '.fig']));
        print(gcf, '-r600', '-dtiff', fullfile(outputDir, ['QC_velocityCoefficientOfVariation_comparison_' tag '.tif']));

        save(fullfile(outputDir, ['QC_R2_CV_values_' tag '.mat']), ...
            'r2NoAO', 'r2AO', 'cvNoAO', 'cvAO', 'labels', 'NormToNoAO', 'xPositionsR2');
    end
end

%% Save the combined quality-control summary
% Combine results from both datasets in one MAT file.
if ~isempty(allQcSummary)
    outputDir = fullfile(pathall_data2024_0, sprintf('FlowResults_v4_avg%d_adaptiveShift_allData', numavgs));
    if exist(outputDir, 'dir') ~= 7
        mkdir(outputDir);
    end
    save(fullfile(outputDir, 'flow_qc_summary_v4.mat'), ...
        'allQcSummary', 'NormToNoAO', 'numavgs', 'shiftamtCandidates', ...
        'pathall_data2024_0', 'pathall_data2024_1');
end

%% Local analysis functions
% Choose the shift, build parameters, and summarize fit quality.
function [shiftamt, scores] = bestShiftamtForAO(imageLines, frmRate, pxlSize, ...
    numavgs, shiftamtCandidates, outputDir, data2024, folderind, filenamebase)
scores = nan(size(shiftamtCandidates));
messages = strings(size(shiftamtCandidates));
fprintf('\nSelecting shiftamt for data2024=%d data%d %s AO with fixed numavgs=%d...\n', ...
    data2024, folderind, filenamebase, numavgs);

parfor i = 1:numel(shiftamtCandidates)
    params = flowParams(numavgs, shiftamtCandidates(i), frmRate, size(imageLines, 2), ...
        outputDir, data2024, folderind, filenamebase, 'AO_shift_selection', false);
    try
        result = flow_quant_func(imageLines, frmRate, pxlSize, params);
        scores(i) = median(result.goodness, 'omitnan');
        messages(i) = sprintf('  shiftamt=%d, median goodness=%.4f', ...
            shiftamtCandidates(i), scores(i));
    catch err
        messages(i) = sprintf('  shiftamt=%d failed: %s', shiftamtCandidates(i), err.message);
    end
end
fprintf('%s\n', messages);

scoresForMax = scores;
scoresForMax(~isfinite(scoresForMax)) = -inf;
if all(scoresForMax == -inf)
    error('flow_quant_main_v4:NoValidShiftamt', ...
        'No shiftamt candidate produced finite AO goodness for data2024=%d data%d %s.', ...
        data2024, folderind, filenamebase);
end
[~, bestInd] = max(scoresForMax);
shiftamt = shiftamtCandidates(bestInd);
fprintf('Selected shiftamt=%d for %s AO, median goodness=%.4f\n', ...
    shiftamt, filenamebase, scoresForMax(bestInd));
end

function params = flowParams(numavgs, shiftamt, frmRate, numColumns, ...
    outputDir, data2024, folderind, filenamebase, kymoName, makePlots)
params = struct('numavgs', numavgs, 'shiftamt', shiftamt, 'skipamt', 2, ...
    'windowsize', frmRate / 2, 'numstd', 3, 'maxGaussWidth', 20, ...
    'setNaN', true, 'useMedianFilter', true, 'startColumn', 1, ...
    'endColumn', numColumns, 'makePlots', makePlots, 'outputDir', outputDir, ...
    'showMedianR2Label', true, ...
    'savePrefix', sprintf('data2024_%d_data%d_%s_%s_avg%d_shift%d', ...
    data2024, folderind, filenamebase, kymoName, numavgs, shiftamt));
end

function qc = flowQc(result)
goodness = result.goodness(:);
velocity = result.velocity_mm_s(:);
meanVelocity = mean(velocity, 'omitnan');
qc.fitR2_median = median(goodness, 'omitnan');
if any(isfinite(goodness))
    qc.fitR2_lowFrac = sum(goodness(isfinite(goodness)) < 0.8) / nnz(isfinite(goodness));
else
    qc.fitR2_lowFrac = nan;
end
if isfinite(meanVelocity) && abs(meanVelocity) > eps
    qc.cvVelocity = std(velocity, 'omitnan') / abs(meanVelocity);
else
    qc.cvVelocity = nan;
end
end

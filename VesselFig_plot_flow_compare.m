%% Select flow results
% Choose the dataset and flow-analysis version to compare.
clear; clc; close all;

data2024 = 1;
resultVersion = 4; % The result is from flow_quant_main_v4.m
numavgs = 50;
if data2024 == 1
    pathall = 'H:\FACED_AO\20241209_M420_Vessel\FACED\';
else
    pathall = 'K:\AO_FACED2026\20260416_RBead_vessel\FACED\';
end
switch resultVersion
    case 3
        shiftamt = 1;
        pathfolder = fullfile(pathall, sprintf('FlowResults_v3_avg%d_shift%d_data2024_%d', ...
            numavgs, shiftamt, data2024));
    case 4
        pathfolder = fullfile(pathall, sprintf('FlowResults_v4_avg%d_adaptiveShift_data2024_%d_20260528', ...
            numavgs, data2024));
    otherwise
        error('plot_temp:UnsupportedResultVersion', ...
            'resultVersion must be 3 or 4.');
end
cd(pathfolder);

files = dir(fullfile(pathfolder, '*_flow_quant.mat'));
baseNames = {};
r2AO_median = [];
r2NoAO_median = [];
r2AO_mean = [];
r2NoAO_mean = [];

%% Pair AO and NoAO flow results
% Load fit R2 statistics and retain complete ROI pairs.
for fileInd = 1:numel(files)
    token = regexp(files(fileInd).name, ...
        '^(.*)_(AO|NoAO)(?:_avg\d+_shift\d+)?_flow_quant\.mat$', ...
        'tokens', 'once');
    if isempty(token)
        continue;
    end

    baseName = token{1};
    pairInd = find(strcmp(baseNames, baseName), 1);
    if isempty(pairInd)
        baseNames{end + 1, 1} = baseName;
        r2AO_median(end + 1, 1) = nan;
        r2NoAO_median(end + 1, 1) = nan;
        r2AO_mean(end + 1, 1) = nan;
        r2NoAO_mean(end + 1, 1) = nan;
        pairInd = numel(baseNames);
    end

    resultFile = fullfile(files(fileInd).folder, files(fileInd).name);
    [fitR2Median, fitR2Mean] = flowResultR2Stats(resultFile);
    if strcmp(token{2}, 'AO')
        r2AO_median(pairInd) = fitR2Median;
        r2AO_mean(pairInd) = fitR2Mean;
    else
        r2NoAO_median(pairInd) = fitR2Median;
        r2NoAO_mean(pairInd) = fitR2Mean;
    end
end

validPairs = isfinite(r2AO_median) & isfinite(r2NoAO_median);
r2AO_median = r2AO_median(validPairs);
r2NoAO_median = r2NoAO_median(validPairs);
r2AO_mean = r2AO_mean(validPairs);
r2NoAO_mean = r2NoAO_mean(validPairs);
baseNames = baseNames(validPairs);

%% Compare AO and NoAO fit quality
% Plot paired median and mean R2 statistics with significance tests.
rng('default');
figure; set(gcf, 'Position', [100 100 500 420], 'Color', 'w'); hold on;
r2Diff = r2AO_median - r2NoAO_median;
[~, p] = ttest(r2Diff, 0, 'Tail', 'right');
scatter(1 + (rand(size(r2Diff)) - 0.5) * 0.2, r2Diff, 20, 'r', 'filled', 'MarkerFaceAlpha', 0.6);
errorbar(1.2, mean(r2Diff, 'omitnan'), std(r2Diff, 'omitnan'), 'ro', 'MarkerFaceColor', 'r');
yline(0, '--', 'Color', [0.4 0.4 0.4]);
xlim([0.5 1.5]); xticks(1); xticklabels({'AO - NoAO'});
ylabel('Median fit R^2 difference');
title(['Median R^2, p = ' num2str(p, '%0.3g')]);
set(gca, 'FontSize', 16);
text(0.62, 0.92, ['p = ' num2str(p, '%0.3g')], 'Units', 'normalized');
text(0.62, 0.84, ['N = ' num2str(numel(r2Diff))], 'Units', 'normalized');
text(0.62, 0.76, ['mean = ' num2str(mean(r2Diff, 'omitnan'), '%0.3f')], 'Units', 'normalized');
text(0.62, 0.68, ['std = ' num2str(std(r2Diff, 'omitnan'), '%0.3f')], 'Units', 'normalized');
savefig(gcf, fullfile(pathfolder, 'QC_R2_AO_minus_NoAO_median.fig'));
print(gcf, '-r600', '-dtiff', fullfile(pathfolder, 'QC_R2_AO_minus_NoAO_median.tif'));
% Show the two paired median R^2 distributions without normalization.
rng('default');
xNoAO = 1.1 + (rand(size(r2NoAO_median)) - 0.5) * 0.2;
xAO = 1.9 + (rand(size(r2AO_median)) - 0.5) * 0.2;
[~, p] = ttest(r2NoAO_median, r2AO_median, 'Tail', 'left');
figure; set(gcf, 'Position', [100 100 500 420], 'Color', 'w'); hold on;
scatter(xNoAO, r2NoAO_median, 20, 'b', 'filled', 'MarkerFaceAlpha', 0.6);
scatter(xAO, r2AO_median, 20, 'r', 'filled', 'MarkerFaceAlpha', 0.6);
errorbar(0.9, mean(r2NoAO_median, 'omitnan'), std(r2NoAO_median, 'omitnan'), ...
    'bo', 'MarkerFaceColor', 'b');
errorbar(2.1, mean(r2AO_median, 'omitnan'), std(r2AO_median, 'omitnan'), ...
    'ro', 'MarkerFaceColor', 'r');
xlim([0.5 2.5]); xticks([1 2]); xticklabels({'NoAO', 'AO'});
ylim([0 1]);
ylabel('Median fit R^2');
title(['Median R^2, p = ' num2str(p, '%0.3g')]);
set(gca, 'FontSize', 16);
text(0.64, 0.27, ['N = ' num2str(numel(r2NoAO_median))]);
text(0.64, 0.19, ['mean = ' num2str(mean(r2NoAO_median, 'omitnan'), '%0.3f')]);
text(0.64, 0.11, ['std = ' num2str(std(r2NoAO_median, 'omitnan'), '%0.3f')]);
text(1.64, 0.27, ['N = ' num2str(numel(r2AO_median))]);
text(1.64, 0.19, ['mean = ' num2str(mean(r2AO_median, 'omitnan'), '%0.3f')]);
text(1.64, 0.11, ['std = ' num2str(std(r2AO_median, 'omitnan'), '%0.3f')]);
savefig(gcf, fullfile(pathfolder, 'QC_R2_AO_and_NoAO_median.fig'));
print(gcf, '-r600', '-dtiff', fullfile(pathfolder, 'QC_R2_AO_and_NoAO_median.tif'));

figure; set(gcf, 'Position', [100 100 500 420], 'Color', 'w'); hold on;
r2Diff = r2AO_mean - r2NoAO_mean;
[~, p] = ttest(r2Diff, 0, 'Tail', 'right');
scatter(1 + (rand(size(r2Diff)) - 0.5) * 0.2, r2Diff, 20, 'r', 'filled', 'MarkerFaceAlpha', 0.6);
errorbar(1.2, mean(r2Diff, 'omitnan'), std(r2Diff, 'omitnan'), 'ro', 'MarkerFaceColor', 'r');
yline(0, '--', 'Color', [0.4 0.4 0.4]);
xlim([0.5 1.5]); xticks(1); xticklabels({'AO - NoAO'});
ylabel('Mean fit R^2 difference');
title(['Mean R^2, p = ' num2str(p, '%0.3g')]);
set(gca, 'FontSize', 16);
text(0.62, 0.92, ['p = ' num2str(p, '%0.3g')], 'Units', 'normalized');
text(0.62, 0.84, ['N = ' num2str(numel(r2Diff))], 'Units', 'normalized');
text(0.62, 0.76, ['mean = ' num2str(mean(r2Diff, 'omitnan'), '%0.3f')], 'Units', 'normalized');
text(0.62, 0.68, ['std = ' num2str(std(r2Diff, 'omitnan'), '%0.3f')], 'Units', 'normalized');
savefig(gcf, fullfile(pathfolder, 'QC_R2_AO_minus_NoAO_mean.fig'));
print(gcf, '-r600', '-dtiff', fullfile(pathfolder, 'QC_R2_AO_minus_NoAO_mean.tif'));

%% Collect selected manuscript examples
% Extract R2 values for the predefined example ROIs.
if data2024 == 0
    exampleGroup = ["3"; "4"; "5"; "6"; "7"; "8"];
    exampleData = [14; 14; 4; 9; 11; 11];
    exampleROI = ["ROI_18_2"; "ROI_17_2"; "ROI_7_3"; "ROI_14_1"; "ROI_7_1"; "ROI_7_2"];
else
    exampleGroup = ["1"; "2"];
    exampleData = [4; 4];
    exampleROI = ["ROI_51_1";"ROI_17_1"];
end

matched = false(numel(exampleROI), 1);
matchedName = strings(numel(exampleROI), 1);
exampleR2AO_median = nan(numel(exampleROI), 1);
exampleR2NoAO_median = nan(numel(exampleROI), 1);
exampleR2AO_mean = nan(numel(exampleROI), 1);
exampleR2NoAO_mean = nan(numel(exampleROI), 1);

for exampleInd = 1:numel(exampleROI)
    currInd = find(contains(string(baseNames), "data" + exampleData(exampleInd) + "_") & ...
        contains(string(baseNames), exampleROI(exampleInd)), 1);
    if isempty(currInd)
        continue;
    end
    matched(exampleInd) = true;
    matchedName(exampleInd) = string(baseNames{currInd});
    exampleR2AO_median(exampleInd) = r2AO_median(currInd);
    exampleR2NoAO_median(exampleInd) = r2NoAO_median(currInd);
    exampleR2AO_mean(exampleInd) = r2AO_mean(currInd);
    exampleR2NoAO_mean(exampleInd) = r2NoAO_mean(currInd);
end

exampleR2Table = table(repmat(data2024, numel(exampleROI), 1), exampleGroup, exampleData, ...
    exampleROI, matched, matchedName, exampleR2AO_median, exampleR2NoAO_median, ...
    exampleR2AO_mean, exampleR2NoAO_mean, ...
    'VariableNames', {'Data2024', 'Group', 'Data', 'ROI', 'Matched', 'MatchedName', ...
    'R2AO_Median', 'R2NoAO_Median', 'R2AO_Mean', 'R2NoAO_Mean'});

disp('Selected R2 comparison examples:');
disp(exampleR2Table);

%% AO velocity versus shiftamt
% Summarize AO velocity as a function of the selected shift amount.
close all;
r2Threshold = 0.5;
aoFiles = files(contains({files.name}, '_AO') & contains({files.name}, '_flow_quant.mat'));
aoShiftamt = nan(numel(aoFiles), 1);
aoMedianVelocity = nan(numel(aoFiles), 1);
aoMeanVelocity = nan(numel(aoFiles), 1);
aoMaxVelocity = nan(numel(aoFiles), 1);
aoVelocityAboveR2Threshold = [];

for fileInd = 1:numel(aoFiles)
    load(fullfile(aoFiles(fileInd).folder, aoFiles(fileInd).name), 'result', 'params');
    goodInd = result.goodness(:) > r2Threshold;
    goodVelocity = result.velocity_mm_s(goodInd);
    aoShiftamt(fileInd) = params.shiftamt;
    aoMedianVelocity(fileInd) = median(goodVelocity, 'omitnan');
    aoMeanVelocity(fileInd) = mean(goodVelocity, 'omitnan');
    aoMaxVelocity(fileInd) = max(goodVelocity, [], 'omitnan');
    aoVelocityAboveR2Threshold = [aoVelocityAboveR2Threshold; goodVelocity(:)]; %#ok<AGROW>
end

validAO = isfinite(aoShiftamt) & isfinite(aoMedianVelocity);
shiftVals = unique(aoShiftamt(validAO));
medianVelocityMean = arrayfun(@(x) mean(aoMedianVelocity(validAO & aoShiftamt == x), ...
    'omitnan'), shiftVals);
medianVelocityStd = arrayfun(@(x) std(aoMedianVelocity(validAO & aoShiftamt == x), ...
    'omitnan'), shiftVals);
figure; set(gcf, 'Position', [100 100 500 420], 'Color', 'w'); hold on;
scatter(aoShiftamt(validAO), aoMedianVelocity(validAO), 28, 'r', 'filled', ...
    'MarkerFaceAlpha', 0.6);
errorbar(shiftVals, medianVelocityMean, medianVelocityStd, 'ko', ...
    'MarkerFaceColor', 'k', 'LineWidth', 1.5, 'CapSize', 10);
xlabel('shiftamt');
ylabel(sprintf('Median velocity for R^2 > %0.2f (mm/s)', r2Threshold));
title('AO median velocity vs. shiftamt');
xlim([0 6]); xticks(1:5);
set(gca, 'FontSize', 16);
grid off; box on;
savefig(gcf, fullfile(pathfolder, 'QC_AO_median_velocity_vs_shiftamt.fig'));
print(gcf, '-r600', '-dtiff', fullfile(pathfolder, 'QC_AO_median_velocity_vs_shiftamt.tif'));

validAO = isfinite(aoShiftamt) & isfinite(aoMeanVelocity);
shiftVals = unique(aoShiftamt(validAO));
meanVelocityMean = arrayfun(@(x) mean(aoMeanVelocity(validAO & aoShiftamt == x), ...
    'omitnan'), shiftVals);
meanVelocityStd = arrayfun(@(x) std(aoMeanVelocity(validAO & aoShiftamt == x), ...
    'omitnan'), shiftVals);
figure; set(gcf, 'Position', [100 100 500 420], 'Color', 'w'); hold on;
scatter(aoShiftamt(validAO), aoMeanVelocity(validAO), 28, 'r', 'filled', ...
    'MarkerFaceAlpha', 0.6);
errorbar(shiftVals, meanVelocityMean, meanVelocityStd, 'ko', ...
    'MarkerFaceColor', 'k', 'LineWidth', 1.5, 'CapSize', 10);
xlabel('shiftamt');
ylabel(sprintf('Mean velocity for R^2 > %0.2f (mm/s)', r2Threshold));
title('AO mean velocity vs. shiftamt');
xlim([0 6]); xticks(1:5);
set(gca, 'FontSize', 16);
grid off; box on;
savefig(gcf, fullfile(pathfolder, 'QC_AO_mean_velocity_vs_shiftamt.fig'));
print(gcf, '-r600', '-dtiff', fullfile(pathfolder, 'QC_AO_mean_velocity_vs_shiftamt.tif'));

validAO = isfinite(aoShiftamt) & isfinite(aoMaxVelocity);
shiftVals = unique(aoShiftamt(validAO));
maxVelocityMean = arrayfun(@(x) mean(aoMaxVelocity(validAO & aoShiftamt == x), ...
    'omitnan'), shiftVals);
maxVelocityStd = arrayfun(@(x) std(aoMaxVelocity(validAO & aoShiftamt == x), ...
    'omitnan'), shiftVals);
figure; set(gcf, 'Position', [100 100 500 420], 'Color', 'w'); hold on;
scatter(aoShiftamt(validAO), aoMaxVelocity(validAO), 28, 'r', 'filled', ...
    'MarkerFaceAlpha', 0.6);
errorbar(shiftVals, maxVelocityMean, maxVelocityStd, 'ko', ...
    'MarkerFaceColor', 'k', 'LineWidth', 1.5, 'CapSize', 10);
xlabel('shiftamt');
ylabel(sprintf('Maximum velocity for R^2 > %0.2f (mm/s)', r2Threshold));
title('AO maximum velocity vs. shiftamt');
xlim([0 6]); xticks(1:5);
set(gca, 'FontSize', 16);
grid off; box on;
savefig(gcf, fullfile(pathfolder, 'QC_AO_max_velocity_vs_shiftamt.fig'));
print(gcf, '-r600', '-dtiff', fullfile(pathfolder, 'QC_AO_max_velocity_vs_shiftamt.tif'));

disp([num2str(min(abs(aoVelocityAboveR2Threshold))) ' - ' num2str(max(abs(aoVelocityAboveR2Threshold))) ' mm/s'])

%% Mean flow velocity across all ROIs
% Report AO and NoAO velocity distributions for both datasets.
data2024List = [1 0];
meanFlowSummary = table();

for dataInd = 1:numel(data2024List)
    currData2024 = data2024List(dataInd);
    if currData2024 == 1
        currPathAll = 'H:\FACED_AO\20241209_M420_Vessel\FACED\';
    else
        currPathAll = 'K:\AO_FACED2026\20260416_RBead_vessel\FACED\';
    end

    currPathFolder = fullfile(currPathAll, sprintf('FlowResults_v4_avg%d_adaptiveShift_data2024_%d_20260528', ...
        numavgs, currData2024));
    currFiles = dir(fullfile(currPathFolder, '*_flow_quant.mat'));
    currKymoName = strings(numel(currFiles), 1);
    currMeanVelocity = nan(numel(currFiles), 1);

    for fileInd = 1:numel(currFiles)
        token = regexp(currFiles(fileInd).name, '_(AO|NoAO)(?:_avg\d+_shift\d+)?_flow_quant\.mat$', ...
            'tokens', 'once');
        if isempty(token)
            continue;
        end

        load(fullfile(currFiles(fileInd).folder, currFiles(fileInd).name), 'result');
        currKymoName(fileInd) = token{1};
        currMeanVelocity(fileInd) = mean(result.velocity_mm_s, 'omitnan');
    end

    validVelocity = currKymoName ~= "" & isfinite(currMeanVelocity);
    for kymoName = ["AO", "NoAO"]
        isKymo = validVelocity & currKymoName == kymoName;
        meanFlowSummary = [meanFlowSummary; table(currData2024, kymoName, nnz(isKymo), ...
            mean(currMeanVelocity(isKymo), 'omitnan'), std(currMeanVelocity(isKymo), 'omitnan'), ...
            min(currMeanVelocity(isKymo), [], 'omitnan'), max(currMeanVelocity(isKymo), [], 'omitnan'), ...
            'VariableNames', {'Data2024', 'KymoName', 'NumROIs', ...
            'MeanVelocity_mm_s', 'StdVelocity_mm_s', ...
            'MinMeanVelocity_mm_s', 'MaxMeanVelocity_mm_s'})]; %#ok<AGROW>
    end
end

disp('Mean flow velocity across all ROIs:');
disp(meanFlowSummary);

%% Load fit-quality statistics
% Read median and mean R2 values from each saved flow result.
function [fitR2Median, fitR2Mean] = flowResultR2Stats(resultFile)
S = load(resultFile, 'result', 'qc');
fitR2Median = nan;
fitR2Mean = nan;

if isfield(S, 'result') && isfield(S.result, 'goodness')
    goodness = S.result.goodness(:);
    fitR2Median = median(goodness, 'omitnan');
    fitR2Mean = mean(goodness, 'omitnan');
    return;
end

if isfield(S, 'qc')
    if isfield(S.qc, 'fitR2_median')
        fitR2Median = S.qc.fitR2_median;
    end
    if isfield(S.qc, 'fitR2_mean')
        fitR2Mean = S.qc.fitR2_mean;
    else
        fitR2Mean = fitR2Median;
    end
end
end

%% Format the active figure
% Apply the standard manuscript figure layout.
xlim([0.7 1.3])
set(gcf,'Position',[100 100 450 420]);
set(gca,'Position',[0.13 0.11 0.775 0.81]);
styleSavedProfileAxes(gca);

%% Example figure
% Plot velocity and fit quality over time for one AO/NoAO example.
clear; clc; close all;
pathfolder = 'K:\AO_FACED2026\20260416_RBead_vessel\FACED\FlowResults_v4_avg50_adaptiveShift_data2024_0_20260528\';
aodatastr = [pathfolder 'data2024_0_data11_ROI_7_1_AO_avg50_shift1_flow_quant.mat'];
noaodatastr = [pathfolder 'data2024_0_data11_ROI_7_1_NoAO_avg50_shift1_flow_quant.mat'];
savestr = 'data11_ROI_7_1';
r2ylim = [0 1]; velocity_ylim = [-2 4.5];
% pathfolder = 'H:\FACED_AO\20241209_M420_Vessel\FACED\FlowResults_v4_avg50_adaptiveShift_data2024_1_20260528\';
% aodatastr = [pathfolder 'data2024_1_data4_ROI_17_1_AO_avg50_shift1_flow_quant.mat'];
% noaodatastr = [pathfolder 'data2024_1_data4_ROI_17_1_NoAO_avg50_shift1_flow_quant.mat'];
% savestr = 'data4_ROI_17_1';
% r2ylim = [-0.2 1]; velocity_ylim = [-1.25 1.25];

load([aodatastr],'result');
figure; clf; set(gcf,'Position',[100 390 560 260]);
set(gca,'Position',[0.16 0.25 0.75 0.72]);
styleSavedProfileAxes(gca);
yyaxis left;
plot(result.time_ms/1000,result.velocity_mm_s,'r');
xlim([result.time_ms(1) result.time_ms(end)]/1000); ylim(velocity_ylim);
xlabel('t (s)'); ylabel('Velocity (mm/s)');
yyaxis right;
plot(result.time_ms/1000,result.goodness,'k');
ylabel('Flow-fit R^2'); ylim(r2ylim);
yyaxis left;
set(gcf,'Color',[1,1,1]); set(gca,'FontSize',16);
grid off; set(gca,'LineWidth',1);
set(gcf,'Position',[100 390 560 260]);
set(gca,'Position',[0.16 0.25 0.75 0.6]);
ax = gca; ax.Box = 'off'; ax.YAxis(1).Color = 'k'; ax.YAxis(2).Color = 'k';
title('AO');
text(0.1, 1, sprintf('Flow-fit R^2 median = %0.2f\nFlow-fit R^2 range = %0.2f–%0.2f', ...
    median(result.goodness), min(result.goodness), max(result.goodness)), ...
    'Interpreter', 'tex');

load([noaodatastr],'result');
figure; clf; set(gcf,'Position',[100 390 560 260]);
set(gca,'Position',[0.16 0.25 0.75 0.72]);
styleSavedProfileAxes(gca);
yyaxis left;
plot(result.time_ms/1000,result.velocity_mm_s,'b');
xlim([result.time_ms(1) result.time_ms(end)]/1000); ylim(velocity_ylim);
xlabel('t (s)'); ylabel('Velocity (mm/s)');
yyaxis right;
plot(result.time_ms/1000,result.goodness,'k');
ylabel('Flow-fit R^2'); ylim(r2ylim);
yyaxis left;
set(gcf,'Color',[1,1,1]); set(gca,'FontSize',16);
grid off; set(gca,'LineWidth',1);
set(gcf,'Position',[100 390 560 260]);
set(gca,'Position',[0.16 0.25 0.75 0.6]);
ax = gca; ax.Box = 'off'; ax.YAxis(1).Color = 'k'; ax.YAxis(2).Color = 'k';
title('NoAO');
text(0.1, 1, sprintf('Flow-fit R^2 median = %0.2f\nFlow-fit R^2 range = %0.2f–%0.2f', ...
    median(result.goodness), min(result.goodness), max(result.goodness)), ...
    'Interpreter', 'tex');

% saveFigs_mod(pathfolder,1,'.fig',savestr);
% saveFigs_mod(pathfolder,1,'.pdf',savestr);

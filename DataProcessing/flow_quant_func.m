function result = flow_quant_func(imageLines, frmRate, pxlSize, params)
%FLOW_QUANT_FUNC Quantify line-scan flow from an in-memory kymograph.
%
% result = flow_quant_func(imageLines, frmRate, pxlSize, params)
%
% imageLines : time x space kymograph
% frmRate    : line/frame rate in Hz
% pxlSize    : spatial pixel size in um along the kymograph axis
% params     : optional struct controlling LSPIV settings

if nargin < 4 || isempty(params)
    params = struct();
end

validateattributes(imageLines, {'double', 'single', 'uint8', 'uint16'}, {'2d', 'nonempty'}, ...
    mfilename, 'imageLines', 1);
validateattributes(frmRate, {'double'}, {'scalar', 'positive'}, mfilename, 'frmRate', 2);
validateattributes(pxlSize, {'double'}, {'scalar', 'positive'}, mfilename, 'pxlSize', 3);

defaults = struct();
defaults.numavgs = 50;
defaults.shiftamt = 2;
defaults.skipamt = 2;
defaults.windowsize = frmRate / 2;
defaults.numstd = 3;
defaults.maxGaussWidth = 20;
defaults.setNaN = true;
defaults.useMedianFilter = true;
defaults.startColumn = 1;
defaults.endColumn = size(imageLines, 2);
defaults.makePlots = true;
defaults.outputDir = pwd;
defaults.savePrefix = 'flow_quant';

params = applyDefaults(params, defaults);

imageLines = double(imageLines);
numLines = size(imageLines, 1);
numPixels = size(imageLines, 2);

startColumn = max(1, params.startColumn);
endColumn = min(numPixels, params.endColumn);
if startColumn >= endColumn
    error('flow_quant_func:InvalidColumns', 'startColumn must be smaller than endColumn.');
end

if numLines <= params.shiftamt + params.numavgs
    error('flow_quant_func:KymographTooShort', ...
        'Kymograph has too few lines for shiftamt=%d and numavgs=%d.', ...
        params.shiftamt, params.numavgs);
end

if params.makePlots && ~exist(params.outputDir, 'dir')
    mkdir(params.outputDir);
end

fprintf('DC correction\n');
dcOffset = mean(imageLines, 1);
imageLinesDC = imageLines - dcOffset;

fprintf('LSPIV begin\n');
scene_fft = fft(imageLinesDC(1:end - params.shiftamt, :), [], 2);
test_img = zeros(size(scene_fft));
test_img(:, startColumn:endColumn) = ...
    imageLinesDC(params.shiftamt + 1:end, startColumn:endColumn);
test_fft = fft(test_img, [], 2);

W = 1 ./ sqrt(abs(scene_fft)) ./ sqrt(abs(test_fft));

lspivFFT = scene_fft .* conj(test_fft) .* W;
lspivResult = ifft(lspivFFT, [], 2);
fprintf('LSPIV complete\n');

fprintf('Find the peaks\n');
index_vals = params.skipamt:params.skipamt:(size(lspivResult, 1) - params.numavgs);
maxPixelShift = round(numPixels / 2) - 1;
centerPixel = round(numPixels / 2);
fitStart = max(1, centerPixel - maxPixelShift);
fitEnd = min(numPixels, centerPixel + maxPixelShift);

velocity = nan(size(index_vals));
amps = nan(size(index_vals));
sigmas = nan(size(index_vals));
goodness = nan(size(index_vals));

for index = 1:length(index_vals)
    if mod(index_vals(index), 100) == 0
        fprintf('line: %d\n', index_vals(index));
    end

    rows = index_vals(index):index_vals(index) + params.numavgs;
    avgCorr = real(fftshift(sum(lspivResult(rows, :), 1))) ...
        / max(real(sum(lspivResult(rows, :), 1)));
    if any(~isfinite(avgCorr))
        continue;
    end

    c = zeros(1, numPixels);
    c(fitStart:fitEnd) = avgCorr(fitStart:fitEnd);
    [~, maxIndex] = max(c);

    options = fitoptions('gauss1');
    options.Lower = [0, centerPixel - maxPixelShift, 0, 0];
    options.Upper = [1e9, centerPixel + maxPixelShift, params.maxGaussWidth, 1];
    options.StartPoint = [1, maxIndex, 10, 0.1];

    try
        [q, good] = fit((1:length(avgCorr))', avgCorr', ...
            'a1*exp(-((x-b1)/c1)^2) + d1', options);
        velocity(index) = (q.b1 - numPixels / 2 - 1) / params.shiftamt;
        amps(index) = q.a1;
        sigmas(index) = q.c1;
        goodness(index) = good.rsquare;
    catch
        velocity(index) = nan;
    end
end

medV = median(velocity);
velocity = medV / abs(medV) * velocity;

badvals = findBadVelocityPoints(velocity, params.windowsize, params.skipamt, params.numstd);
goodvals = setdiff(1:length(velocity), badvals);
if isempty(goodvals) && ~isempty(velocity)
    warning('flow_quant_func:AllVelocityPointsFlagged', ...
        ['All velocity points were flagged as bad. Keeping the fitted velocity ', ...
        'trace for conversion instead of replacing every point with NaN.']);
    badvals = [];
    goodvals = 1:length(velocity);
end
meanvel = mean(velocity(goodvals));
stdvel = std(velocity(goodvals));

time_ms = (index_vals + params.numavgs / 2) / frmRate * 1000;
time_bad_ms = time_ms(badvals);

velocity_filtered = velocity;
if params.setNaN
    velocity_filtered(badvals) = nan;
end
if params.useMedianFilter
    velocity_filtered = medfilt1(velocity_filtered, 3);
end

velocity_mm_s = pxlSize / 1000 * frmRate * resampleVelocityTrace(velocity_filtered);
raw_velocity_mm_s = pxlSize / 1000 * frmRate * velocity;

N = length(velocity_mm_s);
T = numLines / frmRate;
Fs = N / T;
f = Fs * (0:N / 2) / N;
Y1 = abs(fft(velocity_mm_s, N)) / N;
P1 = Y1(1:floor(N / 2) + 1);

result = struct();
result.velocity_pixels_per_scan = velocity;
result.velocity_mm_s = velocity_mm_s;
result.raw_velocity_mm_s = raw_velocity_mm_s;
result.time_ms = time_ms;
result.time_bad_ms = time_bad_ms;
result.badvals = badvals;
result.goodvals = goodvals;
result.meanVelocity_pixels_per_scan = meanvel;
result.stdVelocity_pixels_per_scan = stdvel;
result.meanVelocity_mm_s = mean(velocity_mm_s);
result.stdVelocity_mm_s = std(velocity_mm_s);
result.index_vals = index_vals;
result.startColumn = startColumn;
result.endColumn = endColumn;
result.amps = amps;
result.sigmas = sigmas;
result.goodness = goodness;
result.frequency_Hz = f;
result.frequencyAmplitude = P1;

if params.makePlots
    plotFlowQuantification(imageLines, lspivResult, numPixels, result, params);
end
end

function velocityOut = resampleVelocityTrace(velocityIn)
if isempty(velocityIn)
    velocityOut = velocityIn;
    return;
end

try
    velocityOut = resample(velocityIn, 1:length(velocityIn));
catch
    warning('flow_quant_func:ResampleFailed', ...
        'resample failed, likely because the velocity trace contains too many NaN values. Returning the unresampled trace.');
    velocityOut = velocityIn;
end
end

function params = applyDefaults(params, defaults)
fields = fieldnames(defaults);
for i = 1:numel(fields)
    if ~isfield(params, fields{i}) || isempty(params.(fields{i}))
        params.(fields{i}) = defaults.(fields{i});
    end
end
end

function badvals = findBadVelocityPoints(velocity, windowsize, skipamt, numstd)
pixelWindowSize = max(1, round(windowsize / skipamt));
badpixels = zeros(size(velocity));

for index = 1:max(length(velocity) - pixelWindowSize, 1)
    try
        windowVals = velocity(index:index + pixelWindowSize - 1);
    catch
        windowVals = velocity;
    end
    pmean = mean(windowVals);
    pstd = std(windowVals);
    pbadpts = find(windowVals > pmean + pstd * numstd | ...
        windowVals < pmean - pstd * numstd);
    badpixels(index + pbadpts - 1) = badpixels(index + pbadpts - 1) + 1;
end

badvals = find(badpixels > 0);
end

function plotFlowQuantification(imageLines, lspivResult, numPixels, result, params)
summaryFig = figure(100);
clf(summaryFig);
set(summaryFig, 'Color', 'w', 'Units', 'Normalized', ...
    'Position', [0.05 0.08 0.86 0.82]);

layout = tiledlayout(summaryFig, 3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(layout, sprintf('%s flow quantification, mean %.3f mm/s', ...
    params.savePrefix, result.meanVelocity_mm_s), 'Interpreter', 'none');

timeFull_ms = (0:size(imageLines, 1) - 1) / numel(imageLines(:, 1)) * max(result.time_ms);
spacePixels = result.startColumn:result.endColumn;
displacementAxis = linspace(-numPixels / 2, numPixels / 2, numPixels);

nexttile(layout, 1);
imagesc(timeFull_ms, spacePixels, imageLines(:, result.startColumn:result.endColumn)');
axis tight;
colormap(gca, 'gray');
colorbar;
title('Raw kymograph');
xlabel('Time (ms)');
ylabel('Space (pixel)');

nexttile(layout, 2);
imagesc(result.index_vals, displacementAxis, real(fftshift(lspivResult, 2))');
axis tight;
colorbar;
title('LSPIV cross-correlation');
xlabel('Time index');
ylabel({'Displacement'; '(pixels/scan)'});

nexttile(layout, 3);
yyaxis left;
plot(result.index_vals, result.velocity_pixels_per_scan, '.', 'MarkerSize', 8);
hold on;
plot(result.index_vals(result.badvals), result.velocity_pixels_per_scan(result.badvals), ...
    'ro', 'MarkerSize', 5);
yline(result.meanVelocity_pixels_per_scan, '--k');
xlim([result.index_vals(1), result.index_vals(end)]);
if isfinite(result.meanVelocity_pixels_per_scan) && isfinite(result.stdVelocity_pixels_per_scan)
    ylim(result.meanVelocity_pixels_per_scan + result.stdVelocity_pixels_per_scan * [-4, 4]);
end
ylabel({'Displacement'; '(pixels/scan)'});
yyaxis right;
plot(result.index_vals, result.goodness, '-', 'LineWidth', 1);
ylim([-1 1]);
ylabel('Fit R^2');
hold off;
grid on;
if isfield(params, 'showMedianR2Label') && params.showMedianR2Label
    title(sprintf('Fitted pixel displacement, median R^2 = %.3f', ...
        median(result.goodness, 'omitnan')));
else
    title('Fitted pixel displacement');
end
xlabel('Time index');

nexttile(layout, 4);
plot(result.time_ms, result.velocity_mm_s, '-', 'Marker', '.', 'MarkerSize', 5);
hold on;
plot(result.time_bad_ms, result.velocity_mm_s(result.badvals), 'ro', 'MarkerSize', 5);
yline(result.meanVelocity_mm_s, '--k');
hold off;
xlim([result.time_ms(1), result.time_ms(end)]);
grid on;
title('Velocity');
xlabel('Time (ms)');
ylabel('Velocity (mm/s)');

nyquist = numel(result.velocity_mm_s) / (max(result.time_ms) / 1000) / 2;
fLim = min(nyquist, 40);
P1_norm = mat2gray(smoothdata(result.frequencyAmplitude, 'gaussian', 3));

nexttile(layout, 5);
plot(result.frequency_Hz, P1_norm, 'LineWidth', 2);
axis([0 fLim 0 inf]);
grid on;
title('Frequency spectrum');
xlabel('Frequency (Hz)');
ylabel('Amplitude (a.u.)');

nexttile(layout, 6);
semilogy(result.frequency_Hz, P1_norm, 'LineWidth', 2);
axis([0 fLim 1e-3 inf]);
grid on;
title('Frequency spectrum - semilog');
xlabel('Frequency (Hz)');
ylabel('Amplitude (a.u.)');

savefig(summaryFig, fullfile(params.outputDir, [params.savePrefix '_all_results.fig']));
print(summaryFig, '-r600', '-dtiff', fullfile(params.outputDir, [params.savePrefix '_all_results.tif']));
end

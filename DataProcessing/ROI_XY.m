%% XY ROI selection in MATLAB
% Required variables from the parent workspace:
% pathfolder, Img_AO, Img_NoAO_reg, FOV_x, FOV_y, str1, str2
pplotXY = [pathfolder 'ProfileXY_bicubic_202606\'];
if ~exist(pplotXY, 'dir'); mkdir(pplotXY); end
close all;

N_display = 100;
NFACED = 100;
if exist('dX_FACED', 'var') && ~isempty(dX_FACED)
    dX_XY = dX_FACED;
else
    dX_XY = FOV_x/NFACED;
end
zAxisYZ = 1:size(Img_AO,3);
zAxisLabel = 'Slice number';
if exist('dZ_FACED', 'var') && ~isempty(dZ_FACED)
    dZ_YZ = dZ_FACED;
elseif exist('dep_ao', 'var') && ~isempty(dep_ao) && numel(dep_ao) >= 2
    dZ_YZ = abs(mean(diff(dep_ao)));
else
    dZ_YZ = 1;
end

datainterp_ao = [];
datainterp_noao = [];
Line_ao = [];
Line_noao = [];
xx_all = [];
slice_ao_all = [];
slice_noao_all = [];
savedXYAO = {};
savedXYNoAO = {};
savedRoiLabels = [];
[pendingImportAO, pendingImportNoAO, pendingImportLabels, roiIdx, ...
    savedXYAO, savedXYNoAO, savedRoiLabels, ...
    datainterp_ao, datainterp_noao, Line_ao, Line_noao, ...
    xx_all, slice_ao_all, slice_noao_all] = initializeExistingROIs(pplotXY);
if isinf(roiIdx)
    return;
end

sliceAO = max(1, round(size(Img_AO, 3)/2));
sliceNoAO = max(1, round(size(Img_NoAO_reg, 3)/2));
Img_NoAO_MIPXY = max(Img_NoAO_reg, [], 3);
Img_AO_MIPXY = max(Img_AO, [], 3);
LimXY1 = 0;
LimXY2 = max([double(Img_AO(:)); double(Img_NoAO_reg(:)); ...
    double(Img_AO_MIPXY(:)); double(Img_NoAO_MIPXY(:))]);
if LimXY2 <= 0, LimXY2 = 4000; end

f = figure('Color','w');
set(f,'unit','centimeter','Position',[1 1 52 20]);
ax1 = axes('Parent',f,'Position',[0.015 0.5 0.4 0.4]); % NoAO XY
ax2 = axes('Parent',f,'Position',[0.25 0.5 0.4 0.4]); % AO XY
ax4 = axes('Parent',f,'Position',[0.6 0.61 0.12 0.24]); % NoAO YZ
ax8 = axes('Parent',f,'Position',[0.015 0.05 0.4 0.4]); % NoAO XY MIP
ax7 = axes('Parent',f,'Position',[0.25 0.05 0.4 0.4]); % AO XY MIP
ax5 = axes('Parent',f,'Position',[0.6 0.21 0.12 0.24]); % AO YZ
ax3 = axes('Parent',f,'Position',[0.75 0.68 0.20 0.22]); % XY profile
ax6 = axes('Parent',f,'Position',[0.75 0.39 0.20 0.20]); % YZ focus profile
ax9 = axes('Parent',f,'Position',[0.75 0.10 0.20 0.20]); % selected YZ line profile
setappdata(f, 'sliceAO', sliceAO);
setappdata(f, 'sliceNoAO', sliceNoAO);
setappdata(f, 'LimXY2', LimXY2);
setappdata(f, 'setupAction', 'wait');
setappdata(f, 'savedXYAO', savedXYAO);
setappdata(f, 'savedXYNoAO', savedXYNoAO);
setappdata(f, 'savedRoiLabels', savedRoiLabels);

hNoAOMinus = uicontrol(f,'Style','pushbutton','String','-','Units','normalized', ...
    'Position',[0.05 0.93 0.035 0.04]);
hSliceNoAO = uicontrol(f,'Style','edit','String',num2str(sliceNoAO),'Units','normalized', ...
    'Position',[0.09 0.93 0.06 0.04]);
hNoAOPlus = uicontrol(f,'Style','pushbutton','String','+','Units','normalized', ...
    'Position',[0.155 0.93 0.035 0.04]);

hAOMinus = uicontrol(f,'Style','pushbutton','String','-','Units','normalized', ...
    'Position',[0.31 0.93 0.035 0.04]);
hSliceAO = uicontrol(f,'Style','edit','String',num2str(sliceAO),'Units','normalized', ...
    'Position',[0.35 0.93 0.06 0.04]);
hAOPlus = uicontrol(f,'Style','pushbutton','String','+','Units','normalized', ...
    'Position',[0.415 0.93 0.035 0.04]);
set(hNoAOMinus, 'Callback', @(src,evt) changeSlice(f, ax1, ax2, ax7, ax8, hSliceNoAO, hSliceAO, Img_NoAO_reg, Img_AO, Img_NoAO_MIPXY, Img_AO_MIPXY, str2, str1, 'noao', -1));
set(hSliceNoAO, 'Callback', @(src,evt) setSliceFromEdit(f, ax1, ax2, ax7, ax8, hSliceNoAO, hSliceAO, Img_NoAO_reg, Img_AO, Img_NoAO_MIPXY, Img_AO_MIPXY, str2, str1, 'noao'));
set(hNoAOPlus, 'Callback', @(src,evt) changeSlice(f, ax1, ax2, ax7, ax8, hSliceNoAO, hSliceAO, Img_NoAO_reg, Img_AO, Img_NoAO_MIPXY, Img_AO_MIPXY, str2, str1, 'noao', 1));
set(hAOMinus, 'Callback', @(src,evt) changeSlice(f, ax1, ax2, ax7, ax8, hSliceNoAO, hSliceAO, Img_NoAO_reg, Img_AO, Img_NoAO_MIPXY, Img_AO_MIPXY, str2, str1, 'ao', -1));
set(hSliceAO, 'Callback', @(src,evt) setSliceFromEdit(f, ax1, ax2, ax7, ax8, hSliceNoAO, hSliceAO, Img_NoAO_reg, Img_AO, Img_NoAO_MIPXY, Img_AO_MIPXY, str2, str1, 'ao'));
set(hAOPlus, 'Callback', @(src,evt) changeSlice(f, ax1, ax2, ax7, ax8, hSliceNoAO, hSliceAO, Img_NoAO_reg, Img_AO, Img_NoAO_MIPXY, Img_AO_MIPXY, str2, str1, 'ao', 1));
uicontrol(f,'Style','pushbutton','String','Use slices','Units','normalized', ...
    'Position',[0.59 0.93 0.08 0.04], ...
    'Callback',@(src,evt) finishSetup(f, 'use'));
uicontrol(f,'Style','pushbutton','String','Bright +','Units','normalized', ...
    'Position',[0.68 0.93 0.07 0.04], ...
    'Callback',@(src,evt) changeBrightness(f, ax1, ax2, ax7, ax8, hSliceNoAO, hSliceAO, Img_NoAO_reg, Img_AO, Img_NoAO_MIPXY, Img_AO_MIPXY, str2, str1, 1/1.25));
uicontrol(f,'Style','pushbutton','String','Bright -','Units','normalized', ...
    'Position',[0.76 0.93 0.07 0.04], ...
    'Callback',@(src,evt) changeBrightness(f, ax1, ax2, ax7, ax8, hSliceNoAO, hSliceAO, Img_NoAO_reg, Img_AO, Img_NoAO_MIPXY, Img_AO_MIPXY, str2, str1, 1.25));
uicontrol(f,'Style','pushbutton','String','Finish','Units','normalized', ...
    'Position',[0.84 0.93 0.06 0.04], ...
    'Callback',@(src,evt) finishSetup(f, 'finish'));

while true
    skipToNextROI = false;
    roiWasSaved = false;
    if isappdata(f, 'currentROI')
        rmappdata(f, 'currentROI');
    end

    refreshXYImages(f, ax1, ax2, ax7, ax8, hSliceNoAO, hSliceAO, Img_NoAO_reg, Img_AO, Img_NoAO_MIPXY, Img_AO_MIPXY, str2, str1);
    cla(ax3,'reset');
    title(ax3,sprintf('ROI %d XY profile preview', roiIdx));
    cla(ax4,'reset');
    title(ax4,sprintf('ROI %d %s YZ preview', roiIdx, str2));
    cla(ax5,'reset');
    title(ax5,sprintf('ROI %d %s YZ preview', roiIdx, str1));
    cla(ax6,'reset');
    title(ax6,sprintf('ROI %d YZ profile preview', roiIdx));
    cla(ax9,'reset');
    title(ax9,sprintf('ROI %d selected YZ profile', roiIdx));
    setupAction = getappdata(f, 'setupAction');

    okXY = false;
    while ~okXY
        setupAction = getappdata(f, 'setupAction');
        if strcmp(setupAction, 'finish')
            skipToNextROI = true;
            break;
        end
        if ~isempty(pendingImportAO)
            lineAOImport = pendingImportAO{1};
            lineNoAOImport = pendingImportNoAO{1};
            roiIdx = pendingImportLabels(1);
            pendingImportAO(1) = [];
            pendingImportNoAO(1) = [];
            pendingImportLabels(1) = [];

            xao = [lineAOImport(1,1); lineAOImport(end,1)];
            yao = [lineAOImport(1,2); lineAOImport(end,2)];
            xnoao = [lineNoAOImport(1,1); lineNoAOImport(end,1)];
            ynoao = [lineNoAOImport(1,2); lineNoAOImport(end,2)];

            plot(ax7, lineAOImport(:,1), lineAOImport(:,2), 'y-', 'LineWidth', 2);
            plot(ax2, lineAOImport(:,1), lineAOImport(:,2), 'y-', 'LineWidth', 2);
            plot(ax8, lineNoAOImport(:,1), lineNoAOImport(:,2), 'y-', 'LineWidth', 2);
            plot(ax1, lineNoAOImport(:,1), lineNoAOImport(:,2), 'y-', 'LineWidth', 2);
            title(ax7,sprintf('%s XY imported ROI %d', str1, roiIdx));
            title(ax1,sprintf('%s XY imported ROI %d', str2, roiIdx));
            okXY = true;
            break;
        end
        xySourceChoice = menu(sprintf('ROI %d: select sampleAO XY line from?', roiIdx), ...
            'sampleAO MIP', ...
            'sampleAO XY slice', ...
            'Finish');
        if xySourceChoice == 3
            skipToNextROI = true;
            setupAction = 'finish';
            break;
        elseif xySourceChoice == 2
            axSelectXY = ax2;
            selectTitle = sprintf('%s XY slice %d (ROI %d): click two points, Enter', str1, sliceAO, roiIdx);
        else
            axSelectXY = ax7;
            selectTitle = sprintf('%s XY MIP (ROI %d): click two points, Enter', str1, roiIdx);
        end
        axes(axSelectXY);
        title(axSelectXY, selectTitle);
        try
            [xao, yao, ~] = improfile;
        catch ME
            if contains(ME.message, 'Interruption during mouse selection')
                setupAction = getappdata(f, 'setupAction');
                if strcmp(setupAction, 'finish')
                    skipToNextROI = true;
                    break;
                end
                continue;
            end
            rethrow(ME);
        end
        if numel(xao) < 2 || numel(yao) < 2
            q = menu('Invalid sampleAO XY line', 'Retry', 'Finish');
            if q == 2, skipToNextROI = true; setupAction = 'finish'; break; end
            continue;
        end

        xao = [xao(1); xao(end)];
        yao = [yao(1); yao(end)];
        hAOMIPxy = plot(ax7, xao, yao, 'y-', 'LineWidth', 2);
        hAOxy = plot(ax2, xao, yao, 'y-', 'LineWidth', 2);

        xnoao = xao;
        ynoao = yao;
        hNoAOxy = plot(ax1, xnoao, ynoao, 'y-', 'LineWidth', 2);
        hNoAOMIPxy = plot(ax8, xnoao, ynoao, 'y-', 'LineWidth', 2);
        title(ax1,sprintf('%s XY: translation only', str2));

        while true
            r = menu('sysAO XY line', ...
                'Click translated START point', ...
                'Keep XY ROI pair', ...
                'Redo sampleAO XY line', ...
                'Reselect MIP ROI', ...
                'Finish');
            if r == 1
                axes(ax1);
                try
                    [xs, ys] = ginput(1);
                catch ME
                    if contains(ME.message, 'Interruption during mouse selection')
                        setupAction = getappdata(f, 'setupAction');
                        if strcmp(setupAction, 'finish')
                            skipToNextROI = true;
                            break;
                        end
                        continue;
                    end
                    rethrow(ME);
                end
                if isempty(xs) || isempty(ys)
                    continue;
                end
                dx = xs - xao(1); dy = ys - yao(1);
                xnoao = xao + dx; ynoao = yao + dy;
                if isgraphics(hNoAOxy), delete(hNoAOxy); end
                if isgraphics(hNoAOMIPxy), delete(hNoAOMIPxy); end
                hNoAOxy = plot(ax1, xnoao, ynoao, 'y-', 'LineWidth', 2);
                hNoAOMIPxy = plot(ax8, xnoao, ynoao, 'y-', 'LineWidth', 2);
            elseif r == 2
                okXY = true;
                break;
            elseif r == 3
                if isgraphics(hAOMIPxy), delete(hAOMIPxy); end
                if isgraphics(hAOxy), delete(hAOxy); end
                if isgraphics(hNoAOxy), delete(hNoAOxy); end
                if isgraphics(hNoAOMIPxy), delete(hNoAOMIPxy); end
                break;
            elseif r == 4
                skipToNextROI = true;
                if isgraphics(hAOMIPxy), delete(hAOMIPxy); end
                if isgraphics(hAOxy), delete(hAOxy); end
                if isgraphics(hNoAOxy), delete(hNoAOxy); end
                if isgraphics(hNoAOMIPxy), delete(hNoAOMIPxy); end
                break;
            else
                skipToNextROI = true;
                setupAction = 'finish';
                break;
            end
        end

        if skipToNextROI
            break;
        end
    end

    if strcmp(setupAction, 'finish')
        break;
    end
    if skipToNextROI
        continue;
    end
    setappdata(f, 'currentROI', struct('xao', xao, 'yao', yao, 'xnoao', xnoao, 'ynoao', ynoao));

    [x_line_ao, y_line_ao] = getLine(xao.', yao.', N_display);
    [x_line_noao, y_line_noao] = getLine(xnoao.', ynoao.', N_display);

    Lall = sqrt(((abs(x_line_ao(1)-x_line_ao(end))-1)*dX_XY)^2 + ...
        ((abs(y_line_ao(1)-y_line_ao(end))-1)*FOV_y)^2);
    xx = linspace(-Lall/2, Lall/2, N_display);

    YZSlice_AO = SliceXZ(Img_AO, x_line_ao, y_line_ao, N_display, 'bicubic');
    YZSlice_NoAO = SliceXZ(Img_NoAO_reg, x_line_noao, y_line_noao, N_display, 'bicubic');
    zAxisROI = zAxisYZ(1:size(YZSlice_AO,1));
    YZProfile_AO = mean(YZSlice_AO, 2, 'omitnan');
    YZProfile_NoAO = mean(YZSlice_NoAO, 2, 'omitnan');

    isValid = all(isfinite([x_line_ao(:); y_line_ao(:); x_line_noao(:); y_line_noao(:)])) ...
        && numel(unique([x_line_ao(:) y_line_ao(:)], 'rows')) > 1 ...
        && numel(unique([x_line_noao(:) y_line_noao(:)], 'rows')) > 1 ...
        && ~isempty(YZSlice_AO) && ~isempty(YZSlice_NoAO) ...
        && all(isfinite(YZProfile_AO(:))) && all(isfinite(YZProfile_NoAO(:)));
    if ~isValid
        warning('ROI %d is invalid. Reselect XY for ROI %d.', roiIdx, roiIdx);
        retry = menu(sprintf('ROI %d invalid. Reselect?', roiIdx), 'Yes', 'Reselect slices', 'Finish');
        if retry == 1
            continue;
        elseif retry == 2
            continue;
        else
            break;
        end
    end

    LimYZ2 = max([double(YZSlice_AO(:)); double(YZSlice_NoAO(:))]);
    if ~isfinite(LimYZ2) || LimYZ2 <= 0, LimYZ2 = 1; end

    cla(ax4,'reset');
    imagesc(ax4, [min(xx) max(xx)], [zAxisROI(1) zAxisROI(end)], YZSlice_NoAO);
    colormap(ax4, gray); clim(ax4,[0 LimYZ2]);
    set(ax4,'YDir','reverse');
    hold(ax4,'on');
    xlabel(ax4,'Distance (\mum)');
    ylabel(ax4,zAxisLabel);
    title(ax4,sprintf('%s YZ (ROI %d)', str2, roiIdx));
    set(ax4,'LineWidth',1.2);

    cla(ax5,'reset');
    imagesc(ax5, [min(xx) max(xx)], [zAxisROI(1) zAxisROI(end)], YZSlice_AO);
    colormap(ax5, gray); clim(ax5,[0 LimYZ2]);
    set(ax5,'YDir','reverse');
    hold(ax5,'on');
    xlabel(ax5,'Distance (\mum)');
    ylabel(ax5,zAxisLabel);
    title(ax5,sprintf('%s YZ (ROI %d)', str1, roiIdx));
    set(ax5,'LineWidth',1.2);

    yzProfileYMax = max([YZProfile_AO(:); YZProfile_NoAO(:)]);
    if ~isfinite(yzProfileYMax) || yzProfileYMax <= 0, yzProfileYMax = 1; end
    cla(ax6,'reset');
    plot(ax6, zAxisROI, YZProfile_AO, 'r', 'LineWidth', 1.5); hold(ax6,'on');
    plot(ax6, zAxisROI, YZProfile_NoAO, 'b', 'LineWidth', 1.5);
    legend(ax6, str1, str2, 'Location', 'best');
    xlabel(ax6,zAxisLabel);
    ylabel(ax6,'Mean Signal (A.U.)');
    title(ax6,sprintf('YZ Mean Profiles (ROI %d)', roiIdx));
    ylim(ax6,[0 yzProfileYMax]);
    set(ax6,'LineWidth',1.2);
    grid(ax6,'on');
    drawnow;

    slicesAccepted = false;
    while ~slicesAccepted
        sliceAnswer = inputdlg({ ...
            sprintf('%s slice number (1-%d):', str1, size(Img_AO,3)), ...
            sprintf('%s slice number (1-%d):', str2, size(Img_NoAO_reg,3))}, ...
            sprintf('ROI %d: choose XY slices from YZ focus', roiIdx), 1, ...
            {num2str(sliceAO), num2str(sliceNoAO)});
        if isempty(sliceAnswer)
            break;
        end
        requestedSliceAO = str2double(sliceAnswer{1});
        requestedSliceNoAO = str2double(sliceAnswer{2});
        if isfinite(requestedSliceAO)
            sliceAO = min(max(round(requestedSliceAO), 1), size(Img_AO,3));
        end
        if isfinite(requestedSliceNoAO)
            sliceNoAO = min(max(round(requestedSliceNoAO), 1), size(Img_NoAO_reg,3));
        end
        setappdata(f, 'sliceAO', sliceAO);
        setappdata(f, 'sliceNoAO', sliceNoAO);
        refreshXYImages(f, ax1, ax2, ax7, ax8, hSliceNoAO, hSliceAO, Img_NoAO_reg, Img_AO, Img_NoAO_MIPXY, Img_AO_MIPXY, str2, str1);
        cla(ax4,'reset');
        imagesc(ax4, [min(xx) max(xx)], [zAxisROI(1) zAxisROI(end)], YZSlice_NoAO);
        colormap(ax4, gray); clim(ax4,[0 LimYZ2]);
        set(ax4,'YDir','reverse');
        hold(ax4,'on');
        xlabel(ax4,'Distance (\mum)');
        ylabel(ax4,zAxisLabel);
        title(ax4,sprintf('%s YZ (ROI %d)', str2, roiIdx));
        set(ax4,'LineWidth',1.2);

        cla(ax5,'reset');
        imagesc(ax5, [min(xx) max(xx)], [zAxisROI(1) zAxisROI(end)], YZSlice_AO);
        colormap(ax5, gray); clim(ax5,[0 LimYZ2]);
        set(ax5,'YDir','reverse');
        hold(ax5,'on');
        xlabel(ax5,'Distance (\mum)');
        ylabel(ax5,zAxisLabel);
        title(ax5,sprintf('%s YZ (ROI %d)', str1, roiIdx));
        set(ax5,'LineWidth',1.2);
        if sliceNoAO <= numel(zAxisROI)
            plot(ax4, [min(xx) max(xx)], zAxisROI([sliceNoAO sliceNoAO]), 'b--', 'LineWidth', 1.2);
        end
        if sliceAO <= numel(zAxisROI)
            plot(ax5, [min(xx) max(xx)], zAxisROI([sliceAO sliceAO]), 'r--', 'LineWidth', 1.2);
        end

        dataAO = improfile(Img_AO(:,:,sliceAO), x_line_ao, y_line_ao, N_display, 'bicubic');
        dataNoAO = improfile(Img_NoAO_reg(:,:,sliceNoAO), x_line_noao, y_line_noao, N_display, 'bicubic');
        dataAO = dataAO(:).';
        dataNoAO = dataNoAO(:).';
        if ~all(isfinite(dataAO(:))) || ~all(isfinite(dataNoAO(:)))
            warning('ROI %d profile contains invalid values. Reselect slice numbers.', roiIdx);
            continue;
        end

        cla(ax6,'reset');
        plot(ax6, zAxisROI, YZProfile_AO, 'r', 'LineWidth', 1.5); hold(ax6,'on');
        plot(ax6, zAxisROI, YZProfile_NoAO, 'b', 'LineWidth', 1.5);
        if sliceNoAO <= numel(zAxisROI)
            plot(ax6, zAxisROI([sliceNoAO sliceNoAO]), [0 yzProfileYMax], 'b--', 'LineWidth', 1);
        end
        if sliceAO <= numel(zAxisROI)
            plot(ax6, zAxisROI([sliceAO sliceAO]), [0 yzProfileYMax], 'r--', 'LineWidth', 1);
        end
        legend(ax6, str1, str2, 'Location', 'best');
        xlabel(ax6,zAxisLabel);
        ylabel(ax6,'Mean Signal (A.U.)');
        title(ax6,sprintf('YZ Mean Profiles (ROI %d)', roiIdx));
        ylim(ax6,[0 yzProfileYMax]);
        set(ax6,'LineWidth',1.2);
        grid(ax6,'on');

        cla(ax3,'reset');
        plot(ax3, xx, dataAO, 'r', 'LineWidth', 1.5); hold(ax3,'on');
        plot(ax3, xx, dataNoAO, 'b', 'LineWidth', 1.5);
        legend(ax3, str1, str2, 'Location', 'best');
        xlabel(ax3,'Distance (\mum)');
        ylabel(ax3,'Signal (A.U.)');
        title(ax3,sprintf('XY Line Profiles (ROI %d)', roiIdx));
        set(ax3,'LineWidth',1.2);
        profileYMax = max([dataAO(:); dataNoAO(:)]);
        if ~isfinite(profileYMax) || profileYMax <= 0, profileYMax = 1; end
        ylim(ax3,[0 profileYMax]);
        grid(ax3,'on');

        sliceDecision = menu(sprintf('ROI %d: use these slice numbers?', roiIdx), ...
            'Use these slices', ...
            'Reselect XY slice number', ...
            'Reselect XY ROI', ...
            'Finish');
        if sliceDecision == 1
            slicesAccepted = true;
        elseif sliceDecision == 2
            continue;
        elseif sliceDecision == 3
            skipToNextROI = true;
            break;
        else
            setupAction = 'finish';
            break;
        end
    end

    if strcmp(setupAction, 'finish')
        break;
    end
    if skipToNextROI || ~slicesAccepted
        continue;
    end

    YZ_AO = [];
    YZ_NoAO = [];
    YZ_AO_pix = [];
    YZ_NoAO_pix = [];
    YZLine_AO = [];
    YZLine_NoAO = [];
    yy_YZLine = [];
    yy_YZLine_um = [];
    yzLineYMax = [];

    yzLineChoice = menu(sprintf('ROI %d: select YZ line profile?', roiIdx), ...
        'Select YZ line', ...
        'Skip', ...
        'Reselect slices/XY', ...
        'Finish');
    if yzLineChoice == 3
        continue;
    elseif yzLineChoice == 4
        break;
    elseif yzLineChoice == 1
        keepYZLine = false;
        while ~keepYZLine
            okYZLine = false;
            cla(ax9,'reset');
            title(ax9,sprintf('ROI %d selected YZ profile', roiIdx));
            while ~okYZLine
                axes(ax5);
                title(ax5,sprintf('%s YZ (ROI %d): click two points, Enter', str1, roiIdx));
                try
                    [xYZAO, yYZAO, ~] = improfile;
                catch ME
                    if contains(ME.message, 'Interruption during mouse selection')
                        continue;
                    end
                    rethrow(ME);
                end
                if numel(xYZAO) < 2 || numel(yYZAO) < 2
                    q = menu('Invalid sampleAO YZ line', 'Retry', 'Skip YZ line', 'Finish');
                    if q == 2
                        break;
                    elseif q == 3
                        setupAction = 'finish';
                        break;
                    end
                    continue;
                end

                xYZAO = [xYZAO(1); xYZAO(end)];
                yYZAO = [yYZAO(1); yYZAO(end)];
                hAOyz = plot(ax5, xYZAO, yYZAO, 'y-', 'LineWidth', 2);

                xYZNoAO = xYZAO;
                yYZNoAO = yYZAO;
                hNoAOyz = plot(ax4, xYZNoAO, yYZNoAO, 'y-', 'LineWidth', 2);

                while true
                    r = menu('sysAO YZ line', ...
                        'Click translated START point', ...
                        'Keep YZ ROI pair', ...
                        'Redo sampleAO YZ line', ...
                        'Skip YZ line', ...
                        'Finish');
                    if r == 1
                        axes(ax4);
                        try
                            [xs, ys] = ginput(1);
                        catch ME
                            if contains(ME.message, 'Interruption during mouse selection')
                                continue;
                            end
                            rethrow(ME);
                        end
                        if isempty(xs) || isempty(ys)
                            continue;
                        end
                        dx = xs - xYZAO(1); dy = ys - yYZAO(1);
                        xYZNoAO = xYZAO + dx; yYZNoAO = yYZAO + dy;
                        if isgraphics(hNoAOyz), delete(hNoAOyz); end
                        hNoAOyz = plot(ax4, xYZNoAO, yYZNoAO, 'y-', 'LineWidth', 2);
                    elseif r == 2
                        okYZLine = true;
                        break;
                    elseif r == 3
                        if isgraphics(hAOyz), delete(hAOyz); end
                        if isgraphics(hNoAOyz), delete(hNoAOyz); end
                        break;
                    elseif r == 4
                        if isgraphics(hAOyz), delete(hAOyz); end
                        if isgraphics(hNoAOyz), delete(hNoAOyz); end
                        break;
                    else
                        setupAction = 'finish';
                        break;
                    end
                end

                if strcmp(setupAction, 'finish')
                    break;
                end
            end

            if strcmp(setupAction, 'finish')
                break;
            end

            if ~okYZLine
                keepYZLine = true;
                continue;
            end

            colYZAO = interp1(xx, 1:N_display, xYZAO, 'linear', 'extrap');
            rowYZAO = interp1(zAxisROI, 1:numel(zAxisROI), yYZAO, 'linear', 'extrap');
            colYZNoAO = interp1(xx, 1:N_display, xYZNoAO, 'linear', 'extrap');
            rowYZNoAO = interp1(zAxisROI, 1:numel(zAxisROI), yYZNoAO, 'linear', 'extrap');
            N_display_YZLine = max(round(sqrt((colYZAO(1)-colYZAO(2))^2 + (rowYZAO(1)-rowYZAO(2))^2)), 2);

            [colLineAO, rowLineAO] = getLine(colYZAO.', rowYZAO.', N_display_YZLine);
            [colLineNoAO, rowLineNoAO] = getLine(colYZNoAO.', rowYZNoAO.', N_display_YZLine);
            YZLine_AO = SliceXZ(YZSlice_AO, colLineAO, rowLineAO, N_display_YZLine, 'bicubic');
            YZLine_NoAO = SliceXZ(YZSlice_NoAO, colLineNoAO, rowLineNoAO, N_display_YZLine, 'bicubic');
            YZLine_AO = YZLine_AO(:).';
            YZLine_NoAO = YZLine_NoAO(:).';
            yy_YZLine = 1:N_display_YZLine;
            yzLineLength_um = sqrt(((colYZAO(1)-colYZAO(2))*dX_XY)^2 + ((rowYZAO(1)-rowYZAO(2))*dZ_YZ)^2);
            yy_YZLine_um = linspace(-yzLineLength_um/2, yzLineLength_um/2, N_display_YZLine);
            YZ_AO = [xYZAO yYZAO];
            YZ_NoAO = [xYZNoAO yYZNoAO];
            YZ_AO_pix = [colLineAO(:), rowLineAO(:)];
            YZ_NoAO_pix = [colLineNoAO(:), rowLineNoAO(:)];

            cla(ax9,'reset');
            plot(ax9, yy_YZLine_um, YZLine_AO, 'r', 'LineWidth', 1.5); hold(ax9,'on');
            plot(ax9, yy_YZLine_um, YZLine_NoAO, 'b', 'LineWidth', 1.5);
            legend(ax9, str1, str2, 'Location', 'best');
            xlabel(ax9,'Distance (\mum)');
            ylabel(ax9,'Signal (A.U.)');
            title(ax9,sprintf('YZ Line Profiles (ROI %d)', roiIdx));
            yzLineYMax = max([YZLine_AO(:); YZLine_NoAO(:)]);
            if ~isfinite(yzLineYMax) || yzLineYMax <= 0, yzLineYMax = 1; end
            ylim(ax9,[0 yzLineYMax]);
            set(ax9,'LineWidth',1.2);
            grid(ax9,'on');
            drawnow;

            yzReviewChoice = menu(sprintf('ROI %d: keep displayed YZ line profiles?', roiIdx), ...
                'Keep YZ line profiles', ...
                'Reselect YZ line', ...
                'Skip YZ line', ...
                'Finish');
            if yzReviewChoice == 1
                keepYZLine = true;
            elseif yzReviewChoice == 2
                if isgraphics(hAOyz), delete(hAOyz); end
                if isgraphics(hNoAOyz), delete(hNoAOyz); end
                YZ_AO = [];
                YZ_NoAO = [];
                YZ_AO_pix = [];
                YZ_NoAO_pix = [];
                YZLine_AO = [];
                YZLine_NoAO = [];
                yy_YZLine = [];
                yy_YZLine_um = [];
                yzLineYMax = [];
            elseif yzReviewChoice == 3
                if isgraphics(hAOyz), delete(hAOyz); end
                if isgraphics(hNoAOyz), delete(hNoAOyz); end
                YZ_AO = [];
                YZ_NoAO = [];
                YZ_AO_pix = [];
                YZ_NoAO_pix = [];
                YZLine_AO = [];
                YZLine_NoAO = [];
                yy_YZLine = [];
                yy_YZLine_um = [];
                yzLineYMax = [];
                cla(ax9,'reset');
                title(ax9,sprintf('ROI %d selected YZ profile', roiIdx));
                keepYZLine = true;
            else
                setupAction = 'finish';
                break;
            end
        end

        if strcmp(setupAction, 'finish')
            break;
        end
    end

    xyDecision = menu(sprintf('ROI %d: save result?', roiIdx), 'Save', 'No', 'Reselect slices/XY', 'Finish');
    if xyDecision == 3
        continue;
    elseif xyDecision == 4
        break;
    elseif xyDecision == 1
        roiWasSaved = true;
        [savedXYAO, savedXYNoAO, savedRoiLabels] = upsertSavedROI( ...
            savedXYAO, savedXYNoAO, savedRoiLabels, roiIdx, ...
            [x_line_ao(:), y_line_ao(:)], [x_line_noao(:), y_line_noao(:)]);
        setappdata(f, 'savedXYAO', savedXYAO);
        setappdata(f, 'savedXYNoAO', savedXYNoAO);
        setappdata(f, 'savedRoiLabels', savedRoiLabels);

        save(fullfile(pplotXY, sprintf('ROI_%d.mat', roiIdx)), ...
            'roiIdx', ...
            'dataAO','dataNoAO','x_line_ao','y_line_ao','x_line_noao','y_line_noao', ...
            'YZSlice_AO','YZSlice_NoAO','YZProfile_AO','YZProfile_NoAO', ...
            'YZ_AO','YZ_NoAO','YZ_AO_pix','YZ_NoAO_pix', ...
            'YZLine_AO','YZLine_NoAO','yy_YZLine','yy_YZLine_um','yzLineYMax', ...
            'sliceAO','sliceNoAO','xx','zAxisROI','zAxisLabel', ...
            'Lall','N_display','dX_XY','dZ_YZ','FOV_y');

        fSummary = saveCleanSummaryFigure([ax1 ax2 ax8 ax7 ax3 ax4 ax5 ax6 ax9], ...
            fullfile(pplotXY, sprintf('ROI_%d_summary.fig', roiIdx)), ...
            fullfile(pplotXY, sprintf('ROI_%d_summary.tif', roiIdx)));
        if isgraphics(fSummary), close(fSummary); end

        yzXDim = abs(diff(xlim(ax4)));
        yzYDim = abs(diff(ylim(ax4))) * dZ_YZ;
        fSysAOYZ = saveYZImagePanel(ax4, yzXDim, yzYDim, ...
            fullfile(pplotXY, sprintf('ROI_%d_sysAO_YZ.fig', roiIdx)), ...
            fullfile(pplotXY, sprintf('ROI_%d_sysAO_YZ.tif', roiIdx)));
        if isgraphics(fSysAOYZ), close(fSysAOYZ); end

        yzXDim = abs(diff(xlim(ax5)));
        yzYDim = abs(diff(ylim(ax5))) * dZ_YZ;
        fSampleAOYZ = saveYZImagePanel(ax5, yzXDim, yzYDim, ...
            fullfile(pplotXY, sprintf('ROI_%d_sampleAO_YZ.fig', roiIdx)), ...
            fullfile(pplotXY, sprintf('ROI_%d_sampleAO_YZ.tif', roiIdx)));
        if isgraphics(fSampleAOYZ), close(fSampleAOYZ); end

        fProfile = saveXYProfileFigure(xx, dataAO, dataNoAO, profileYMax, str1, str2, roiIdx, ...
            fullfile(pplotXY, sprintf('ROI_%d_profile.fig', roiIdx)), ...
            fullfile(pplotXY, sprintf('ROI_%d_profile.tif', roiIdx)));
        if isgraphics(fProfile), close(fProfile); end

        if ~isempty(YZLine_AO) && ~isempty(YZLine_NoAO)
            fYZProfile = saveYZProfileFigure(yy_YZLine_um, YZLine_AO, YZLine_NoAO, yzLineYMax, str1, str2, roiIdx, ...
                fullfile(pplotXY, sprintf('ROI_%d_YZ_profile.fig', roiIdx)), ...
                fullfile(pplotXY, sprintf('ROI_%d_YZ_profile.tif', roiIdx)));
            if isgraphics(fYZProfile), close(fYZProfile); end
        end

    end

    if roiWasSaved
        roiIdx = roiIdx + 1;
    end
end

function [pendingImportAO, pendingImportNoAO, pendingImportLabels, nextRoiIdx, ...
        savedXYAO, savedXYNoAO, savedRoiLabels, ...
        datainterp_ao, datainterp_noao, Line_ao, Line_noao, ...
        xx_all, slice_ao_all, slice_noao_all] = initializeExistingROIs(pplotXY)
    pendingImportAO = {};
    pendingImportNoAO = {};
    pendingImportLabels = [];
    savedXYAO = {};
    savedXYNoAO = {};
    savedRoiLabels = [];
    datainterp_ao = [];
    datainterp_noao = [];
    Line_ao = [];
    Line_noao = [];
    xx_all = [];
    slice_ao_all = [];
    slice_noao_all = [];
    [roiFiles, roiNums, nextRoiIdx] = findExistingROIFiles(pplotXY);
    if isempty(roiFiles)
        return;
    end

    importChoice = menu(sprintf('Found %d existing ROI files.', numel(roiFiles)), ...
        'Load ROI locations and reprocess', ...
        'Load all saved ROI data only', ...
        'Do not load existing ROIs', ...
        'Finish');
    if importChoice == 1
        [pendingImportAO, pendingImportNoAO, pendingImportLabels] = loadExistingXYROIs(roiFiles, roiNums);
        savedXYAO = pendingImportAO;
        savedXYNoAO = pendingImportNoAO;
        savedRoiLabels = pendingImportLabels;
    elseif importChoice == 2
        [savedXYAO, savedXYNoAO, savedRoiLabels, ...
            datainterp_ao, datainterp_noao, Line_ao, Line_noao, ...
            xx_all, slice_ao_all, slice_noao_all] = loadExistingCompleteROIs(roiFiles, roiNums);
    elseif importChoice == 4
        nextRoiIdx = inf;
    end
end

function [roiFiles, roiNums, nextRoiIdx] = findExistingROIFiles(pplotXY)
    nextRoiIdx = 1;
    roiFiles = dir(fullfile(pplotXY, 'ROI_*.mat'));
    roiNums = [];
    if isempty(roiFiles)
        return;
    end

    roiNumsAll = nan(numel(roiFiles), 1);
    for k = 1:numel(roiFiles)
        tok = regexp(roiFiles(k).name, '^ROI_(\d+)\.mat$', 'tokens', 'once');
        if ~isempty(tok)
            roiNumsAll(k) = str2double(tok{1});
        end
    end

    validFiles = isfinite(roiNumsAll);
    roiFiles = roiFiles(validFiles);
    roiNums = roiNumsAll(validFiles);
    if isempty(roiFiles)
        return;
    end

    [roiNums, sortIdx] = sort(roiNums);
    roiFiles = roiFiles(sortIdx);
    nextRoiIdx = max(roiNums) + 1;
end

function [pendingImportAO, pendingImportNoAO, pendingImportLabels] = loadExistingXYROIs(roiFiles, roiNums)
    pendingImportAO = {};
    pendingImportNoAO = {};
    pendingImportLabels = [];

    for k = 1:numel(roiFiles)
        roiPath = fullfile(roiFiles(k).folder, roiFiles(k).name);
        S = load(roiPath, 'x_line_ao', 'y_line_ao', 'x_line_noao', 'y_line_noao', 'Line_ao', 'Line_noao');
        [lineAO, lineNoAO] = extractXYLinesFromROI(S);
        if isempty(lineAO) || isempty(lineNoAO)
            continue;
        end
        pendingImportAO{end+1} = lineAO;
        pendingImportNoAO{end+1} = lineNoAO;
        pendingImportLabels(end+1) = roiNums(k);
    end
end

function [savedXYAO, savedXYNoAO, savedRoiLabels, ...
        datainterp_ao, datainterp_noao, Line_ao, Line_noao, ...
        xx_all, slice_ao_all, slice_noao_all] = loadExistingCompleteROIs(roiFiles, roiNums)
    savedXYAO = {};
    savedXYNoAO = {};
    savedRoiLabels = [];
    datainterp_ao = [];
    datainterp_noao = [];
    Line_ao = [];
    Line_noao = [];
    xx_all = [];
    slice_ao_all = [];
    slice_noao_all = [];

    for k = 1:numel(roiFiles)
        roiPath = fullfile(roiFiles(k).folder, roiFiles(k).name);
        S = load(roiPath);
        [lineAO, lineNoAO] = extractXYLinesFromROI(S);
        if isempty(lineAO) || isempty(lineNoAO)
            continue;
        end

        savedXYAO{end+1} = lineAO;
        savedXYNoAO{end+1} = lineNoAO;
        savedRoiLabels(end+1) = roiNums(k);

        if isfield(S, 'dataAO') && isfield(S, 'dataNoAO') && ...
                ~isempty(S.dataAO) && ~isempty(S.dataNoAO)
            saveIdx = size(datainterp_ao, 1) + 1;
            datainterp_ao(saveIdx,:) = S.dataAO(:).';
            datainterp_noao(saveIdx,:) = S.dataNoAO(:).';
            Line_ao(saveIdx,:,1) = lineAO(:,1).';
            Line_ao(saveIdx,:,2) = lineAO(:,2).';
            Line_noao(saveIdx,:,1) = lineNoAO(:,1).';
            Line_noao(saveIdx,:,2) = lineNoAO(:,2).';
            if isfield(S, 'xx') && ~isempty(S.xx)
                xx_all(saveIdx,:) = S.xx(:).';
            end
            if isfield(S, 'sliceAO') && ~isempty(S.sliceAO)
                slice_ao_all(saveIdx) = S.sliceAO;
            end
            if isfield(S, 'sliceNoAO') && ~isempty(S.sliceNoAO)
                slice_noao_all(saveIdx) = S.sliceNoAO;
            end
        end
    end
end

function [lineAO, lineNoAO] = extractXYLinesFromROI(S)
    lineAO = [];
    lineNoAO = [];
    if isfield(S, 'x_line_ao') && isfield(S, 'y_line_ao') && ...
            isfield(S, 'x_line_noao') && isfield(S, 'y_line_noao')
        lineAO = [S.x_line_ao(:), S.y_line_ao(:)];
        lineNoAO = [S.x_line_noao(:), S.y_line_noao(:)];
    elseif isfield(S, 'Line_ao') && isfield(S, 'Line_noao')
        lineAO = normalizeXYLine(S.Line_ao);
        lineNoAO = normalizeXYLine(S.Line_noao);
    end

    if size(lineAO,2) ~= 2 || size(lineNoAO,2) ~= 2 || size(lineAO,1) < 2 || size(lineNoAO,1) < 2
        lineAO = [];
        lineNoAO = [];
    end
end

function lineXY = normalizeXYLine(lineData)
    lineXY = [];
    if isempty(lineData)
        return;
    end
    if ndims(lineData) == 3 && size(lineData,3) == 2
        lineXY = [squeeze(lineData(1,:,1)).', squeeze(lineData(1,:,2)).'];
    elseif size(lineData,2) == 2
        lineXY = lineData;
    elseif size(lineData,1) == 2
        lineXY = lineData.';
    end
end

function changeSlice(f, ax1, ax2, ax7, ax8, hSliceNoAO, hSliceAO, Img_NoAO_reg, Img_AO, Img_NoAO_MIPXY, Img_AO_MIPXY, str2, str1, whichStack, delta)
    sliceAO = getappdata(f, 'sliceAO');
    sliceNoAO = getappdata(f, 'sliceNoAO');
    if strcmp(whichStack, 'ao')
        sliceAO = min(max(sliceAO + delta, 1), size(Img_AO,3));
    else
        sliceNoAO = min(max(sliceNoAO + delta, 1), size(Img_NoAO_reg,3));
    end
    setappdata(f, 'sliceAO', sliceAO);
    setappdata(f, 'sliceNoAO', sliceNoAO);
    refreshXYImages(f, ax1, ax2, ax7, ax8, hSliceNoAO, hSliceAO, Img_NoAO_reg, Img_AO, Img_NoAO_MIPXY, Img_AO_MIPXY, str2, str1);
end

function setSliceFromEdit(f, ax1, ax2, ax7, ax8, hSliceNoAO, hSliceAO, Img_NoAO_reg, Img_AO, Img_NoAO_MIPXY, Img_AO_MIPXY, str2, str1, whichStack)
    sliceAO = getappdata(f, 'sliceAO');
    sliceNoAO = getappdata(f, 'sliceNoAO');
    if strcmp(whichStack, 'ao')
        requestedSlice = str2double(get(hSliceAO, 'String'));
        if isfinite(requestedSlice)
            sliceAO = min(max(round(requestedSlice), 1), size(Img_AO,3));
        end
    else
        requestedSlice = str2double(get(hSliceNoAO, 'String'));
        if isfinite(requestedSlice)
            sliceNoAO = min(max(round(requestedSlice), 1), size(Img_NoAO_reg,3));
        end
    end
    setappdata(f, 'sliceAO', sliceAO);
    setappdata(f, 'sliceNoAO', sliceNoAO);
    refreshXYImages(f, ax1, ax2, ax7, ax8, hSliceNoAO, hSliceAO, Img_NoAO_reg, Img_AO, Img_NoAO_MIPXY, Img_AO_MIPXY, str2, str1);
end

function changeBrightness(f, ax1, ax2, ax7, ax8, hSliceNoAO, hSliceAO, Img_NoAO_reg, Img_AO, Img_NoAO_MIPXY, Img_AO_MIPXY, str2, str1, scaleFactor)
    LimXY2 = getappdata(f, 'LimXY2');
    LimXY2 = LimXY2 * scaleFactor;
    setappdata(f, 'LimXY2', LimXY2);
    refreshXYImages(f, ax1, ax2, ax7, ax8, hSliceNoAO, hSliceAO, Img_NoAO_reg, Img_AO, Img_NoAO_MIPXY, Img_AO_MIPXY, str2, str1);
end

function finishSetup(f, action)
    setappdata(f, 'setupAction', action);
    uiresume(f);
end

function fSummary = saveCleanSummaryFigure(sourceAxes, figPath, tifPath)
    fSummary = figure('Color','w','unit','centimeter','Position',[1 1 52 20], 'Visible','off');
    for k = 1:numel(sourceAxes)
        axCopy = copyobj(sourceAxes(k), fSummary);
        set(axCopy, 'Units', 'normalized', 'Position', get(sourceAxes(k), 'Position'));
    end
    drawnow;
    safeSaveFigure(fSummary, figPath, 'fig');
    safeSaveFigure(fSummary, tifPath, 'tif');
end

function fProfile = saveXYProfileFigure(xx, dataAO, dataNoAO, profileYMax, str1, str2, roiIdx, figPath, tifPath)
    fProfile = figure('Color','w','Position',[100 100 560 420], 'Visible','off');
    plot(xx, dataAO, 'r', 'LineWidth', 1.5); hold on;
    plot(xx, dataNoAO, 'b', 'LineWidth', 1.5);
    legend(str1, str2, 'Location', 'best');
    xlabel('Distance (\mum)');
    ylabel('Signal (A.U.)');
    title(sprintf('XY Line Profiles (ROI %d)', roiIdx));
    styleSavedProfileAxes(gca);
    xlim([min(xx) max(xx)]);
    ylim([0 profileYMax]);
    grid off;
    safeSaveFigure(fProfile, figPath, 'fig');
    safeSaveFigure(fProfile, tifPath, 'tif');
end

function fYZProfile = saveYZProfileFigure(yy, dataAO, dataNoAO, yzLineYMax, str1, str2, roiIdx, figPath, tifPath)
    fYZProfile = figure('Color','w','Position',[100 100 560 420], 'Visible','off');
    plot(yy, dataAO, 'r', 'LineWidth', 1.5); hold on;
    plot(yy, dataNoAO, 'b', 'LineWidth', 1.5);
    legend(str1, str2, 'Location', 'best');
    xlabel('Distance (\mum)');
    ylabel('Signal (A.U.)');
    title(sprintf('YZ Line Profiles (ROI %d)', roiIdx));
    styleSavedProfileAxes(gca);
    xlim([min(yy) max(yy)]);
    ylim([0 yzLineYMax]);
    grid off;
    safeSaveFigure(fYZProfile, figPath, 'fig');
    safeSaveFigure(fYZProfile, tifPath, 'tif');
end

function fYZImage = saveYZImagePanel(sourceAx, xDim, yDim, figPath, tifPath)
    fYZImage = figure('Color','w','Position',[100 100 560 420], 'Visible','off');
    sourceXLim = xlim(sourceAx);
    sourceYLim = ylim(sourceAx);
    axCopy = copyobj(sourceAx, fYZImage);
    set(axCopy, 'Units', 'normalized', 'Position', [0.12 0.14 0.80 0.78]);
    set(fYZImage, 'CurrentAxes', axCopy);
    imgObj = findobj(axCopy, 'Type', 'image');
    if ~isempty(imgObj)
        imgData = get(imgObj(1), 'CData');
        [nRows, nCols] = size(imgData);
        lineObjs = findobj(axCopy, 'Type', 'line');
        for k = 1:numel(lineObjs)
            set(lineObjs(k), 'XData', axisToPixel(get(lineObjs(k), 'XData'), sourceXLim, nCols));
            set(lineObjs(k), 'YData', axisToPixel(get(lineObjs(k), 'YData'), sourceYLim, nRows));
        end
        set(imgObj, 'XData', [1 nCols], 'YData', [1 nRows]);
        xlim(axCopy, [1 nCols]);
        ylim(axCopy, [1 nRows]);
        xlabel(axCopy, 'Pixel');
        ylabel(axCopy, 'Pixel');
        if isfinite(xDim) && isfinite(yDim) && xDim > 0 && yDim > 0
            pbaspect(axCopy, [xDim yDim 1]);
        end
    end
    if isfinite(xDim) && isfinite(yDim) && xDim > 0 && yDim > 0
        scalebar_mod2(xDim, yDim, 2, '\mum', 'w', 'w');
    end
    axis(axCopy, 'off');
    drawnow;
    safeSaveFigure(fYZImage, figPath, 'fig');
    safeSaveFigure(fYZImage, tifPath, 'tif');
end

function pixelData = axisToPixel(axisData, axisLim, nPix)
    axisMin = min(axisLim);
    axisMax = max(axisLim);
    if axisMax == axisMin || nPix <= 1
        pixelData = axisData;
        return;
    end
    pixelData = 1 + (axisData - axisMin) ./ (axisMax - axisMin) .* (nPix - 1);
end

function safeSaveFigure(figHandle, savePath, saveType)
    [saveFolder, saveName, saveExt] = fileparts(savePath);
    if isempty(saveFolder)
        saveFolder = pwd;
    end
    tmpPath = [tempname saveExt];
    oldVisible = get(figHandle, 'Visible');
    set(figHandle, 'Visible', 'on');
    drawnow;

    cleanupTmp = onCleanup(@() cleanupFigureSaveTmp(tmpPath, figHandle, oldVisible));
    if strcmpi(saveType, 'fig')
        savefig(figHandle, tmpPath, 'compact');
        testFig = openfig(tmpPath, 'invisible');
        if isgraphics(testFig), close(testFig); end
    else
        saveas(figHandle, tmpPath);
    end

    if exist(tmpPath, 'file') ~= 2
        error('Temporary figure save failed: %s', tmpPath);
    end
    if exist(savePath, 'file') == 2
        delete(savePath);
    end
    copyfile(tmpPath, savePath, 'f');
    if strcmpi(saveType, 'fig')
        testFig = openfig(savePath, 'invisible');
        if isgraphics(testFig), close(testFig); end
    end
end

function cleanupFigureSaveTmp(filePath, figHandle, oldVisible)
    if isgraphics(figHandle)
        set(figHandle, 'Visible', oldVisible);
    end
    deleteIfExists(filePath);
end

function deleteIfExists(filePath)
    if exist(filePath, 'file') == 2
        delete(filePath);
    end
end

function refreshXYImages(f, ax1, ax2, ax7, ax8, hSliceNoAO, hSliceAO, Img_NoAO_reg, Img_AO, Img_NoAO_MIPXY, Img_AO_MIPXY, str2, str1)
    sliceAO = getappdata(f, 'sliceAO');
    sliceNoAO = getappdata(f, 'sliceNoAO');
    LimXY2 = getappdata(f, 'LimXY2');
    LimXY1 = 0;
    oldXLim = xlim(ax1);
    oldYLim = ylim(ax1);
    hasView = all(isfinite([oldXLim oldYLim])) ...
        && abs(diff(oldXLim)) > 0 && abs(diff(oldYLim)) > 0 ...
        && min(oldXLim) >= 0.5 && max(oldXLim) <= size(Img_NoAO_reg,2) + 0.5 ...
        && min(oldYLim) >= 0.5 && max(oldYLim) <= size(Img_NoAO_reg,1) + 0.5;
    set(hSliceAO, 'String', num2str(sliceAO));
    set(hSliceNoAO, 'String', num2str(sliceNoAO));

    cla(ax1,'reset');
    imagesc(ax1, Img_NoAO_reg(:,:,sliceNoAO)); colormap(ax1, gray); clim(ax1,[LimXY1 LimXY2]);
    axis(ax1,'image'); title(ax1,sprintf('%s XY slice %d', str2, sliceNoAO)); hold(ax1,'on');

    cla(ax2,'reset');
    imagesc(ax2, Img_AO(:,:,sliceAO)); colormap(ax2, gray); clim(ax2,[LimXY1 LimXY2]);
    axis(ax2,'image'); title(ax2,sprintf('%s XY slice %d', str1, sliceAO)); hold(ax2,'on');

    cla(ax7,'reset');
    imagesc(ax7, Img_AO_MIPXY); colormap(ax7, gray); clim(ax7,[LimXY1 LimXY2]);
    axis(ax7,'image'); title(ax7,sprintf('%s XY MIP', str1)); hold(ax7,'on');

    cla(ax8,'reset');
    imagesc(ax8, Img_NoAO_MIPXY); colormap(ax8, gray); clim(ax8,[LimXY1 LimXY2]);
    axis(ax8,'image'); title(ax8,sprintf('%s XY MIP', str2)); hold(ax8,'on');

    if isappdata(f, 'savedXYAO')
        savedXYAO = getappdata(f, 'savedXYAO');
        savedXYNoAO = getappdata(f, 'savedXYNoAO');
        savedRoiLabels = getappdata(f, 'savedRoiLabels');
        drawSavedROIs(ax2, ax7, savedXYAO, savedRoiLabels);
        drawSavedROIs(ax1, ax8, savedXYNoAO, savedRoiLabels);
    end

    if isappdata(f, 'currentROI')
        roiLine = getappdata(f, 'currentROI');
        plot(ax7, roiLine.xao, roiLine.yao, 'y-', 'LineWidth', 2);
        plot(ax2, roiLine.xao, roiLine.yao, 'y-', 'LineWidth', 2);
        plot(ax8, roiLine.xnoao, roiLine.ynoao, 'y-', 'LineWidth', 2);
        plot(ax1, roiLine.xnoao, roiLine.ynoao, 'y-', 'LineWidth', 2);
    end

    linkaxes([ax1 ax2 ax7 ax8], 'xy');
    if hasView
        xlim(ax1, oldXLim);
        ylim(ax1, oldYLim);
        xlim(ax2, oldXLim);
        ylim(ax2, oldYLim);
        xlim(ax7, oldXLim);
        ylim(ax7, oldYLim);
        xlim(ax8, oldXLim);
        ylim(ax8, oldYLim);
    end
    drawnow;
end

function drawSavedROIs(axSlice, axMIP, savedLines, savedLabels)
    for k = 1:numel(savedLines)
        lineXY = savedLines{k};
        if isempty(lineXY), continue; end
        plot(axSlice, lineXY(:,1), lineXY(:,2), 'c-', 'LineWidth', 1.3);
        plot(axMIP, lineXY(:,1), lineXY(:,2), 'c-', 'LineWidth', 1.3);
        labelX = mean(lineXY(:,1));
        labelY = mean(lineXY(:,2));
        text(axSlice, labelX, labelY, sprintf('%d', savedLabels(k)), ...
            'Color','c','FontSize',10,'FontWeight','bold', ...
            'HorizontalAlignment','center','VerticalAlignment','middle');
        text(axMIP, labelX, labelY, sprintf('%d', savedLabels(k)), ...
            'Color','c','FontSize',10,'FontWeight','bold', ...
            'HorizontalAlignment','center','VerticalAlignment','middle');
    end
end

function [savedXYAO, savedXYNoAO, savedRoiLabels] = upsertSavedROI(savedXYAO, savedXYNoAO, savedRoiLabels, roiIdx, lineAO, lineNoAO)
    existingIdx = find(savedRoiLabels == roiIdx, 1, 'first');
    if isempty(existingIdx)
        savedXYAO{end+1} = lineAO;
        savedXYNoAO{end+1} = lineNoAO;
        savedRoiLabels(end+1) = roiIdx;
    else
        savedXYAO{existingIdx} = lineAO;
        savedXYNoAO{existingIdx} = lineNoAO;
    end
end

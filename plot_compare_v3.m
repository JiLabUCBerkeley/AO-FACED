function [dff_NoAO_norm,dff_AO_norm,h,p,x_positions] = plot_compare_v3(NormToNoAO,dff_NoAO,dff_AO,dffthre,...
    Oneplot,Oneplotind,stdorSE,ylimVal,x_positions)
if NormToNoAO
    [dff_NoAO_norm,dff_AO_norm] = NormToNoAO_fc(dff_NoAO,dff_AO,dffthre);
else
    dff_AO_norm = dff_AO;
    dff_NoAO_norm = dff_NoAO;
end

if nargin < 9 || isempty(x_positions)
    % Generate default if none given
    if NormToNoAO==1
        rng('default')
        x_positions = 1 + (rand(size(dff_AO_norm)) - 0.5) * 0.2;
    else
        rng('default')
        x_positions.NoAO = 1.1 + (rand(size(dff_NoAO_norm)) - 0.5) * 0.2;
        x_positions.AO   = 1.9 + (rand(size(dff_AO_norm)) - 0.5) * 0.2;
    end
else
    % Use provided x_positions
    [];
end

[h,p] = ttest(dff_NoAO_norm,dff_AO_norm,"Tail","left");
if Oneplot==1
    if Oneplotind==1
        figure; set(gcf,'Position',[100 100 2000 420]);
    end
    subplot(1,4,Oneplotind); hold on;
else
    figure; set(gcf,'Position',[100 100 500 420]); hold on;
end

if NormToNoAO==1
    scatter(x_positions, dff_AO_norm, 20, 'r', 'filled', 'MarkerFaceAlpha', 0.6);

    textyloc = ylimVal(1)+linspace(0.05,0.3,5)*(ylimVal(2)-ylimVal(1));
    text(0.7,textyloc(5),['p = ' num2str(p)]);
    text(0.7,textyloc(4),['N ori = ' num2str(length(dff_AO))]);
    text(0.7,textyloc(3),['N norm = ' num2str(length(dff_AO_norm))]);
    text(0.7,textyloc(2),['mean = ' num2str(mean(dff_AO_norm),'%0.3f')]);
    if stdorSE==1
        errorbar(1.2,mean(dff_AO_norm),std(dff_AO_norm), 'ro', 'MarkerFaceColor','r');
        text(0.7,textyloc(1),['std = ' num2str(std(dff_AO_norm),'%0.3f')])
    elseif stdorSE==2
        errorbar(1.2,mean(dff_AO_norm),standardError(dff_AO_norm), 'ro', 'MarkerFaceColor','r');
        text(0.7,textyloc(1),['SE = ' num2str(standardError(dff_AO_norm),'%0.3f')]);
    end
elseif NormToNoAO==0
    scatter(x_positions.NoAO, dff_NoAO_norm, 20, 'b', 'filled', 'MarkerFaceAlpha', 0.6);

    scatter(x_positions.AO, dff_AO_norm, 20, 'r', 'filled', 'MarkerFaceAlpha', 0.6);

    textyloc = ylimVal(1)+linspace(0.05,0.3,5)*(ylimVal(2)-ylimVal(1));
    text(0.7,textyloc(5),['p = ' num2str(p)]);
    text(0.7,textyloc(4),['N ori = ' num2str(length(dff_NoAO))]);
    text(0.7,textyloc(3),['N norm = ' num2str(length(dff_NoAO_norm))]);
    text(0.7,textyloc(2),['mean = ' num2str(mean(dff_NoAO_norm),'%0.3f')]);

    text(1.5,textyloc(4),['N ori = ' num2str(length(dff_AO))]);
    text(1.5,textyloc(3),['N norm = ' num2str(length(dff_AO_norm))]);
    text(1.5,textyloc(2),['mean = ' num2str(mean(dff_AO_norm),'%0.3f')]);

    if stdorSE==1
        errorbar(0.9,mean(dff_NoAO_norm),std(dff_NoAO_norm), 'bo', 'MarkerFaceColor','b');
        errorbar(2.1,mean(dff_AO_norm),std(dff_AO_norm), 'ro', 'MarkerFaceColor','r');
        text(0.7,textyloc(1),['std = ' num2str(std(dff_NoAO_norm),'%0.3f')])
        text(1.5,textyloc(1),['std = ' num2str(std(dff_AO_norm),'%0.3f')])
    elseif stdorSE==2
        errorbar(0.9,mean(dff_NoAO_norm),standardError(dff_NoAO_norm), 'bo', 'MarkerFaceColor','b');
        errorbar(2.1,mean(dff_AO_norm),standardError(dff_AO_norm), 'ro', 'MarkerFaceColor','r');
        text(0.7,textyloc(1),['SE = ' num2str(standardError(dff_NoAO_norm),'%0.3f')]);
        text(1.5,textyloc(1),['SE = ' num2str(standardError(dff_AO_norm),'%0.3f')]);
    end
end

set(gcf,'Color',[1,1,1]); set(gca,'FontSize',16);
end
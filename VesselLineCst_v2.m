function [Cst_Norm_AO,Cst_Norm_NoAO,Inten_AO,Inten_NoAO,Inten_RBC_AO,Inten_RBC_NoAO] = VesselLineCst_v2...
    (ROI_AO_line,ROI_NoAO_line,contrast_cal,Nratio)

if ismember(contrast_cal,1:5)
    Inten_AO = mean(maxk(ROI_AO_line,round(length(ROI_AO_line)*Nratio)));
    Inten_NoAO = mean(maxk(ROI_NoAO_line,round(length(ROI_NoAO_line)*Nratio)));

    Inten_RBC_AO = mean(mink(ROI_AO_line,round(length(ROI_AO_line)*Nratio)));
    Inten_RBC_NoAO = mean(mink(ROI_NoAO_line,round(length(ROI_NoAO_line)*Nratio)));
end

if contrast_cal==1 % old plots, RMS contrast
    Cst_Norm_AO = sqrt(1/length(ROI_AO_line)*sum((ROI_AO_line-mean(ROI_AO_line)).^2))/Inten_AO(end);
    Cst_Norm_NoAO = sqrt(1/length(ROI_NoAO_line)*sum((ROI_NoAO_line-mean(ROI_NoAO_line)).^2))/Inten_NoAO(end);
elseif contrast_cal==2 % df/f0
    Cst_Norm_AO = (Inten_AO(end)-Inten_RBC_AO(end))/Inten_RBC_AO(end);
    Cst_Norm_NoAO = (Inten_NoAO(end)-Inten_RBC_NoAO(end))/Inten_RBC_NoAO(end);
elseif contrast_cal==3 % df/f0
    Cst_Norm_AO = mean(maxk((ROI_AO_line-Inten_RBC_AO(end))/Inten_RBC_AO(end),round(length(ROI_AO_line)*Nratio)));
    Cst_Norm_NoAO = mean(maxk((ROI_NoAO_line-Inten_RBC_NoAO(end))/Inten_RBC_NoAO(end),round(length(ROI_NoAO_line)*Nratio)));
elseif contrast_cal==4 % RMS contrast
    Cst_Norm_AO = sqrt(1/length(ROI_AO_line)*sum((ROI_AO_line-mean(ROI_AO_line)).^2))/Inten_RBC_AO(end);
    Cst_Norm_NoAO = sqrt(1/length(ROI_NoAO_line)*sum((ROI_NoAO_line-mean(ROI_NoAO_line)).^2))/Inten_RBC_NoAO(end);
elseif contrast_cal==5
    Cst_Norm_AO = (Inten_AO(end)-Inten_RBC_AO(end))/(Inten_AO(end)+Inten_RBC_AO(end));
    Cst_Norm_NoAO = (Inten_NoAO(end)-Inten_RBC_NoAO(end))/(Inten_NoAO(end)+Inten_RBC_NoAO(end));
elseif contrast_cal==6
    bin_edges = linspace(min(cat(1,ROI_AO_line(:),ROI_NoAO_line(:))), ...
        max(cat(1,ROI_AO_line(:),ROI_NoAO_line(:))),50);

    [Inten_AO,Inten_RBC_AO] = localHistogramPeakLevels(ROI_AO_line,bin_edges);
    [Inten_NoAO,Inten_RBC_NoAO] = localHistogramPeakLevels(ROI_NoAO_line,bin_edges);

    Cst_Norm_AO = localDfOverF0(Inten_AO,Inten_RBC_AO);
    Cst_Norm_NoAO = localDfOverF0(Inten_NoAO,Inten_RBC_NoAO);
end

end
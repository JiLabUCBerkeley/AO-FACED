%% Load deformable-mirror measurements
% Load the mirror maps and optional system-aberration reference.
clear; clc; close all;
load(['C:\Users\JZhu\Box\Matlab\Code\FACED\DMMontage_mask.mat']);
pathname = 'K:\AO_FACED2026\20260416_RBead_vessel\20260416_FAM\sampleTest2\sample_AO4\';
cd(pathname);
lambda = 1035*1e-9; % m

subtractSys = 1;
if subtractSys==1
    psys = 'K:\AO_FACED2026\20260416_RBead_vessel\20260416_FAM\sysAO\sys_AO3\';
    sysname = '0.1_0.1_64_64_HWP26_PMT0.8_RF7_5AVG_BMC_mirror_image.tif';
    Vsys = tiffreadVolume([psys sysname]);
end

%% Display the indexed mirror map
% Show piston and wavefront values as images and a surface.
fileall = dir([pathname '*mirror_image*.tif']);
V = tiffreadVolume([pathname fileall.name]);
h = DMimg_display(imrotate(V,180));
title('Mirror image, indexed');
h.Title.String = "Piston (\mum)";
set(gcf,'Color',[1,1,1]); set(gca,'FontSize',16);

[X,Y] = meshgrid(1:size(V,2),1:size(V,1));
figure; set(gca,'Position',[0.1102    0.1100    0.6    0.8034]);
s = surf(X,Y,imrotate(V,180),'FaceAlpha',0.7); h = colorbar; colormap jet;
h.Location = 'eastoutside';
h.Position = [0.85 0.11 0.0381 0.8034];
s.EdgeColor = 'none';
title('Mirror image, indexed');
h.Title.String = "Piston (\mum)";
Viewagle = [140 45]; view(Viewagle); box on;
set(gcf,'Color',[1,1,1]); set(gca,'FontSize',16);

V = V-min(V(:));
h = DMimg_display(imrotate(V*2*1e-6/lambda,180));
title('Mirror image, indexed');
climv = clim; h.Title.String = "Wave";
set(gcf,'Color',[1,1,1]); set(gca,'FontSize',16);

%% Display the unwrapped mirror map
% Load and visualize the unwrapped piston and wavefront when available.
fileall = dir([pathname '*mirror*Unwrap.tif']);
if length(fileall)==1
    V = tiffreadVolume([pathname fileall.name]);
    h = DMimg_display(imrotate(V,180));
    title('Mirror image, unwrapped');
    climv = clim; h.Title.String = "Piston (\mum)";
    set(gcf,'Color',[1,1,1]); set(gca,'FontSize',16);

    figure; set(gca,'Position',[0.1102    0.1100    0.6    0.8034]);
    s = surf(X,Y,imrotate(V,180),'FaceAlpha',0.7); h = colorbar; colormap jet;
    h.Location = 'eastoutside';
    h.Position = [0.85 0.11 0.0381 0.8034];
    s.EdgeColor = 'none';
    title('Mirror image, unwrapped');
    h.Title.String = "Piston (\mum)";
    Viewagle = [140 45]; view(Viewagle); box on;
    set(gcf,'Color',[1,1,1]); set(gca,'FontSize',16);

    V = V-min(V(:));
    h = DMimg_display(imrotate(V*2*1e-6/lambda,180));
    title('Mirror image, unwrapped');
    climv = clim; h.Title.String = "Wave";
    set(gcf,'Color',[1,1,1]); set(gca,'FontSize',16);
end

%% Subtract the system aberration
% Isolate and display the sample-induced mirror correction.
if subtractSys==1
    VsampleOnly = V-Vsys;
    VsampleOnly = VsampleOnly-min(VsampleOnly(:));
    h = DMimg_display(imrotate(VsampleOnly,180));
    title('Mirror image, sample only');
    climv = clim; h.Title.String = "Piston (\mum)";
    set(gcf,'Color',[1,1,1]); set(gca,'FontSize',16);

    h = DMimg_display(imrotate(VsampleOnly*2*1e-6/lambda,180));
    title('Mirror image, sample only');
    climv = clim; h.Title.String = "Wave";
    set(gcf,'Color',[1,1,1]); set(gca,'FontSize',16);
end

%% Compare system, total, and sample maps
% Display the three aligned aberration maps side by side.
close all;
figure; imagesc(fliplr(imrotate(Vsys,180))); colormap jet; axis square
figure; imagesc(fliplr(imrotate(V,180))); colormap jet; axis square
figure; imagesc(fliplr(imrotate(VsampleOnly,180))); colormap jet; axis square

%% calculate RMS
% Calculate the RMS sample wavefront error in wavelengths.
if subtractSys==1
    AA = imrotate(VsampleOnly*2*1e-6/lambda,180);
else
    AA = imrotate(V*2*1e-6/lambda,180);
end
% rms_val = sqrt(mean((AA).^2, 'all', 'omitnan'));
rms_val = sqrt(mean((AA-nanmean(AA(:))).^2, 'all', 'omitnan'));
disp(['rms val = ' num2str(rms_val) ' λ'])

%% Reconstruct the tilt montage
% Load the full montage and individual DM subregion measurements.
fileall = dir([pathname '*tilt montage*.tif']);
V = tiffreadVolume([pathname fileall.name]);
BWall(BWall==0) = nan;
h = DMimg_display(V.*BWall);
title('Tilt montage, ori');
h.Title.String = "Signal (A.U.)";
set(gcf,'Color',[1,1,1]); set(gca,'FontSize',16);

% Normalize each subregion before rebuilding the actuator map.
DMnum = 37;
filesub = dir([pathname '*Subregion*.tif']);
fileind = [];
for iii = 1:length(filesub)
    filesub(iii).name
    fileind = [fileind str2num(filesub(iii).name(end-6:end-4))+1];
    clear temp;
    temp = tiffreadVolume([pathname filesub(iii).name]);
    if iii==1
        sub_all = zeros(size(temp,1),size(temp,2),length(filesub));
    end
    sub_all(:,:,iii) = temp;
end
sub_norm = sub_all./repmat(max(max(sub_all,[],1),[],2),[size(temp,1),size(temp,2),1]);

% Map detected regions to the corresponding DM actuators.
BW = (V>0);
[imgd,N] = bwlabel(BW);
roimask = zeros(size(imgd,1),size(imgd,2),DMnum);
tilt_montage = nan(size(V));
tilt_montage_norm = nan(size(V));

cal_method = 2;
if N<DMnum
    cal_method = 1;
end

if cal_method==1
    load('C:\Users\JZhu\Box\Matlab\Code\FACED\maskcenter_all.mat');
    indcurr = zeros(1,N);
    for iii = 1:N
        roi_current = imgd;
        roi_current(roi_current~=iii) = nan;
        roi_current(roi_current==iii) = 1;
        temp = find(roi_current==1);
        [rr,cc] = ind2sub([size(imgd,1),size(imgd,2)],temp);

        distall = zeros(1,DMnum);
        for disind = 1:DMnum
            [distall(1,disind)] = dist_cal(nanmean(rr),nanmean(cc),maskcenter_all(1,disind),maskcenter_all(2,disind));
        end
        [~,indcurr(iii)] = min(distall);
        roimask(:,:,indcurr(iii)) = roi_current;
        tilt_montage_norm(min(rr):max(rr),min(cc):max(cc)) = sub_norm(:,:,find(fileind==indcurr(iii)));
    end
elseif cal_method==2
    maskcenter_all = zeros(2,DMnum);
    mapind = [37,35,33,28,...
        36,34,29,25,21,...
        32,30,26,22,18,17,...
        31,27,23,19,14,11,10,...
        24,20,15,12,7,4,...
        16,13,8,5,2,...
        9,6,3,1];
    subind = 0;

    for roiind = mapind
        subind = subind+1;
        roi_current = imgd;
        roi_current(roi_current~=roiind) = nan;
        roi_current(roi_current==roiind) = 1;
        roimask(:,:,subind) = roi_current;
        temp = find(roi_current==1);
        [rr,cc] = ind2sub([size(imgd,1),size(imgd,2)],temp);
        maskcenter_all(:,subind) = round([nanmean(rr);nanmean(cc)]);
        tilt_montage(min(rr):max(rr),min(cc):max(cc)) = imrotate(sub_all(:,:,subind),-90);
        tilt_montage_norm(min(rr):max(rr),min(cc):max(cc)) = imrotate(sub_norm(:,:,subind),-90);
        % small image is also rotated to match it
    end
end

h = DMimg_display(tilt_montage.*BWall);
title('Tilt montage, rotated');
h.Title.String = "Norm signal";
set(gcf,'Color',[1,1,1]); set(gca,'FontSize',16);

h = DMimg_display(tilt_montage_norm.*BWall);
title('Tilt montage, rotated, normalized');
h.Title.String = "Norm signal";
set(gcf,'Color',[1,1,1]); set(gca,'FontSize',16);

%
saveFigs_mod(pathname,1,'.fig','Aberration');

%% Display one tip/tilt map
% Show and save a representative normalized actuator response.
h = DMimg_display(imrotate(sub_norm(:,:,22),-90));
title('Tip/tilt map, normalized');
h.Title.String = "Norm signal";
set(gcf,'Color',[1,1,1]); set(gca,'FontSize',16);

saveFigs_mod(pathname,1,'.fig','TipTilt_1');
saveFigs_mod(pathname,1,'.pdf','TipTilt_1');

%% Adjust montage color scale
% Apply the final display range and color-bar layout.
displayLim = [0 2];
get(gca,'Position'); set(gca,'Position',ans);
h = colorbar; h.Position = [0.82 0.1095 0.02 0.2]; % For zoom-in data
clim(displayLim);
h.Ticks = displayLim;


%% display integrated FACED signal trace
% disaply FFT magnitude.
clear; clc; close all;
pixeltime = 1/100/1000; % s, 100 pixel/ms
dwelltime = 50*pixeltime; % s, 0.5 ms each point
FS = 1/dwelltime; % Hz

pathfolder = 'G:\FACED_AO\20240612_M431_Thy1GFP\20240612_FAM\sysAO\sysAO_FACED_1\';
xlsdata = xlsread([pathfolder '0.1_0.1_64_64_HWP26_PMT0.8_RF16.5_3AVG_Optimize angles_FFTdata_Iteration0_Group1_Tiltstep00.csv']);
txtdata = readmatrix([pathfolder '0.1_0.1_64_64_HWP26_PMT0.8_RF16.5_3AVG_FFTdata_Iteration0_Group1_Tiltstep00.txt']);
% 96 before average: integrate clock

% tip/tilt, measurement data summation is the measurement in the txt file
figure; plot(sum(xlsdata,2)); hold on;
plot(txtdata(:,2));

%
close all;
figure; set(gcf,'Position',[100 200 1800 400]);
subplot(1,3,1); plot(txtdata(:,2));
subplot(1,3,2); plot(fftshift(abs(fft(txtdata(:,2))))); set(gca,'Yscale','log');
subplot(1,3,3); plot(fftshift(angle(fft(txtdata(:,2))))/pi*180); % degrees

freq_p = txtdata(:,3); % fft freq
freq_segDM = txtdata(:,7); % DM modulation freq
[interpoint,ifreq,i2] = intersect(freq_p,freq_segDM);

L = size(txtdata,1);
tx = [0:L-1]*dwelltime;
f = FS/L*(0:(L/2)); f = cat(2,f,nan(1,L-length(f)));

figure; set(gca,'Position',[0.13 0.18 0.775 0.7477]);
plot(tx*1000,txtdata(:,2),'k'); xlim([0 90]);
xlabel('Time (ms)'); ylabel(['Integrate FACED' newline 'fluorescence signal (a.u.)']);
set(gcf,'Color',[1,1,1]); set(gca,'FontSize',16);

figure; set(gca,'Position',[0.13 0.18 0.775 0.7477]);
plot(FS/L*(-L/2:L/2-1),fftshift(abs(fft(txtdata(:,2)))))
xlim([480 940]);

figure; set(gca,'Position',[0.13 0.18 0.775 0.7477]);
plot(f,txtdata(:,4),'k'); hold on;
plot(f(ifreq),txtdata(ifreq,4),'r*');
xlim([480 940]); ylim([0 0.5]);
xlabel('Frequency (Hz)'); ylabel(['FFT magnitude (a.u.)']);
set(gcf,'Color',[1,1,1]); set(gca,'FontSize',16);

figure; set(gca,'Position',[0.13 0.18 0.775 0.7477]);
plot(freq_p,txtdata(:,4),'k'); hold on;
plot(freq_p(ifreq),txtdata(ifreq,4),'r*');
xlim([0.24 0.465]);
xticks([freq_segDM(1:5:20)]);  % Set specific x-tick locations
xticklabels({'\omega_1', '\omega_6', '\omega_{11}', '\omega_{16}'});
xlabel('Freq'); ylabel(['FFT magnitude (a.u.)']);
set(gcf,'Color',[1,1,1]); set(gca,'FontSize',16);

% saveFigs_mod(pathfolder,0,'.fig','AM_demo_tip_tilt_202510');
% saveFigs_mod(pathfolder,0,'.pdf','AM_demo_tip_tilt_202510');

%% Display the FACED phase response
% disaply FFT phase.
close all;
xlsdata_phase = xlsread([pathfolder '0.1_0.1_64_64_HWP26_PMT0.8_RF16.5_3AVG_Optimize rel. phases_FFTdata_Iteration0_Group0_Tiltstep00.csv']);
txtdata_phase = readmatrix([pathfolder '0.1_0.1_64_64_HWP26_PMT0.8_RF16.5_3AVG_FFTdata01_Iteration00.txt']);

L = size(txtdata_phase,1);
tx = [0:L-1]*dwelltime;
f = FS/L*(0:(L/2)); f = cat(2,f,nan(1,L-length(f)));

% figure; plot(angle(fftshift(fft(AA)))/pi*180);

% measurement data sum is the measurement in the txt file
figure; plot(mean(xlsdata_phase,2)/max(mean(xlsdata_phase,2))); hold on;
plot(txtdata_phase(:,2)/max(txtdata_phase(:,2)));

freq_p_phase = txtdata_phase(:,3); % fft freq
freq_segDM_phase = txtdata_phase(:,7); % DM modulation freq
[interpoint,ifreq_phase,i2_phase] = intersect(freq_p_phase,freq_segDM_phase);

figure; set(gca,'Position',[0.13 0.18 0.775 0.7477]);
plot(tx*1000,txtdata_phase(:,2),'k'); xlim([0 180]);
xlabel('Time (ms)'); ylabel(['Integrate FACED' newline 'fluorescence signal (a.u.)']);
set(gcf,'Color',[1,1,1]); set(gca,'FontSize',16);

figure; set(gca,'Position',[0.13 0.18 0.775 0.7477]);
plot(f,txtdata_phase(:,5),'k'); hold on;
plot(f(ifreq_phase),txtdata_phase(ifreq_phase,5),'r*');
xlim([480 940]); ylim([-120 120]);
xlabel('Frequency (Hz)'); ylabel(['FFT phase (degree)']);
set(gcf,'Color',[1,1,1]); set(gca,'FontSize',16);

figure; set(gca,'Position',[0.13 0.18 0.775 0.7477]);
plot(freq_p_phase,txtdata_phase(:,5),'k'); hold on;
plot(freq_p_phase(ifreq_phase),txtdata_phase(ifreq_phase,5),'r*');
xlim([0.24 0.465]);
xticks([freq_segDM_phase(1:5:20)]);  % Set specific x-tick locations
xticklabels({'\omega_1', '\omega_6', '\omega_{11}', '\omega_{16}'});
xlabel('Freq'); ylabel(['FFT phase (a.u.)']);
set(gcf,'Color',[1,1,1]); set(gca,'FontSize',16);

% saveFigs_mod(pathfolder,0,'.fig','AM_demo_phase_202510');
% saveFigs_mod(pathfolder,0,'.pdf','AM_demo_phase_202510');

function [ROI_img,ROI_img_up,ROI_img_up_flat,line_profile] = imgflat(Img,XZ_cord,d1_size,d2_size,s,flatmethod,up_factor)
% changes made: flatten the whole image and then crop to get ROI_img

ROI_img = Img(max(1,round(round(mean(s(:)))-d1_size)):min(round(round(mean(s(:)))+d1_size),size(Img,1)),XZ_cord(1,1):XZ_cord(2,1));

if flatmethod==1 % image rotation is hard to explain the unit
    dy = diff(XZ_cord(:,2));
    dx = diff(XZ_cord(:,1));
    angle = atan2d(dy,dx);
    ROI_img_up_flat = imrotate(ROI_img, angle, 'bilinear', 'crop');

    H_center = mean(XZ_cord(:,1))-XZ_cord(1,1)+1;
    line_profile = mean(ROI_img_up_flat(:,round(max(1,H_center-d2_size):min(H_center+d2_size,size(ROI_img_up_flat,2)))),2);

elseif flatmethod==2
    % No upsampling makes the image step look like, due to round to integer pixels
    % img_up = ifft(fft(ROI_img),up_factor*size(ROI_img,1),1); % does not work well with the vessel data

    s = s-round(mean(s(:))-d1_size)+1;
    ROI_img_up = imresize(ROI_img,[up_factor*size(ROI_img,1) size(ROI_img,2)],'bicubic');
    s = s*up_factor;
    sr = -round(s-mean(s(:)));
    ROI_img_up_flat = zeros(size(ROI_img_up));
    for i=1:length(sr(:))
        ROI_img_up_flat(:,i) = circshift(ROI_img_up(:,i),sr(i));
    end

    H_center = mean(XZ_cord(:,1))-XZ_cord(1,1)+1;
    line_profile = mean(ROI_img_up_flat(:,round(max(1,H_center-d2_size):min(H_center+d2_size,size(ROI_img_up_flat,2)))),2);

elseif flatmethod==3
    Imgcrop = Img(:,XZ_cord(1,1):XZ_cord(2,1));
    Img_up = imresize(Imgcrop,[up_factor*size(Imgcrop,1) size(Imgcrop,2)],'bicubic');
    s = s*up_factor;
    sr = -round(s-mean(s(:)));
    Img_up_flat = zeros(size(Img_up));
    for i=1:length(sr(:))
        Img_up_flat(:,i) = circshift(Img_up(:,i),sr(i));
    end

    ROI_img_up = Img_up(max(1,round(round(mean(s(:)))-d1_size*up_factor)):min(round(round(mean(s(:)))+d1_size*up_factor),size(Img_up,1)),:);
    ROI_img_up_flat = Img_up_flat(max(1,round(round(mean(s(:)))-d1_size*up_factor)):min(round(round(mean(s(:)))+d1_size*up_factor),size(Img_up_flat,1)),:);

    H_center = round(size(ROI_img_up_flat,2)/2);
    line_profile = mean(ROI_img_up_flat(:,round(max(1,H_center-d2_size):min(H_center+d2_size,size(ROI_img_up_flat,2)))),2);
end

end
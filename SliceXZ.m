function [XZSlice] = SliceXZ(Img,x_line,y_line,N_display,interpMethod,win_paral,dx_physical,dy_physical)
% Img: Imaging volume or a 2D frame
% x_line, y_line: the coordinates to do slicing
% N_display: number of points to be displayed, resampling in the improfile step
% interpMethod: interp method in the improfile step
% win_paral: range that's parallel to the line, for averaging the data
% dx_physical,dy_physical: physical distance in X and Y dimension in the image

if nargin<6
    usewindow = 0;
else
    usewindow = 1;
end

if usewindow==1
    XZSlice = zeros(size(Img,3),N_display,length(win_paral));
    parallLines = getparallLines(x_line,y_line,win_paral,dx_physical,dy_physical);
    % get lines that are parallel to [x_line, y_line]
else
    XZSlice = zeros(size(Img,3),N_display);
end

for iii = 1:size(Img,3)
    if usewindow==0
        XZSlice(iii,:) = improfile(Img(:,:,iii),x_line,y_line,N_display,interpMethod);
    elseif usewindow==1
        for wind = 1:length(win_paral)
            x_line = parallLines(:,1,wind);
            y_line = parallLines(:,2,wind);
            XZSlice(iii,:,wind) = improfile(Img(:,:,iii),x_line,y_line,N_display,interpMethod);
        end
    end
end
XZSlice = mean(XZSlice,3);

end
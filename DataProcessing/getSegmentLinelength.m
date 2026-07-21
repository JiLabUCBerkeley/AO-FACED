function [total_length] = getSegmentLinelength(ptsi,dX_FACED,FOV_y)
total_length = 0;  % Initialize total length in pixels
X_all = ptsi(:,1);
Y_all = ptsi(:,2);

for i = 1:length(ptsi)-1  % Loop through all consecutive pairs of points
    % Calculate the Euclidean distance between consecutive points in pixels
    dx = X_all(i+1) - X_all(i);
    dy = Y_all(i+1) - Y_all(i);
    segment_length_pixels = sqrt(dx^2 + dy^2);  % Calculate the distance in pixels

    % Convert the segment length from pixels to physical units (mm, cm, etc.)
    segment_length_physical = segment_length_pixels * sqrt(dX_FACED^2 + FOV_y^2);

    % Add the physical length to the total length
    total_length = total_length + segment_length_physical;
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code generates disparity map using energy minimization with 
% smoothing for defined window size.
% 
% Submitted by: Ashwin Goyal (UID - 115526297)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function disparity_map_energy_minimization(I_left, I_right, window_size)

% Add boundary to apply SSD on the edge pixels
I_left_mod = zeros(window_size-1+size(I_left,1),window_size-1+size(I_left,2));
I_left_mod(ceil(window_size/2):end-floor(window_size/2),ceil(window_size/2):end-floor(window_size/2)) = I_left;
I_right_mod = zeros(window_size-1+size(I_right,1),window_size-1+size(I_right,2));
I_right_mod(ceil(window_size/2):end-floor(window_size/2),ceil(window_size/2):end-floor(window_size/2)) = I_right;

% Compute disparity using SSD
D = zeros(size(I_left));
SSD = zeros([size(I_left_mod),size(I_left,2)]);
for i=ceil(window_size/2):floor(window_size/2)+size(I_left,1)
    for j=ceil(window_size/2):floor(window_size/2)+size(I_left,2)
        % Create left window
        left_window = I_left_mod(i-floor(window_size/2):i+floor(window_size/2),j-floor(window_size/2):j+floor(window_size/2));
        for k=ceil(window_size/2):floor(window_size/2)+size(I_left,2)
            % Create right window
            right_window = I_right_mod(i-floor(window_size/2):i+floor(window_size/2),k-floor(window_size/2):k+floor(window_size/2));
            % Store SSD values
            SSD(i,j,k-floor(window_size/2)) = sum(sum((left_window - right_window).^2));
        end
        % Compute disparity for minimum SSD
        D(i-floor(window_size/2),j-floor(window_size/2)) = abs(find(SSD(i,j,:)==min(SSD(i,j,:)),1) - j);
    end
end

% Apply error smoothing
SSD_new = zeros(size(SSD));
for i=ceil(window_size/2):size(SSD,1)-floor(window_size/2)
    for j=ceil(window_size/2):size(SSD,2)-floor(window_size/2)
        for k=1:size(SSD,3)
            % Create error window
            error_window = SSD(i-floor(window_size/2):i+floor(window_size/2),j-floor(window_size/2):j+floor(window_size/2),k);
            % Compute new SSD as an average
            SSD_new(i,j,k) = sum(sum(error_window));
        end
        % Compute new disparity for new minima
        D(i-floor(window_size/2),j-floor(window_size/2)) = abs(find(SSD_new(i,j,:)==min(SSD_new(i,j,:)),1) - j);
    end
end

% Change range to 0 to 255
D = uint8(round(D*255/max(max(D))));

% Generate a histogram of the disparities
H = zeros(1,256);
for i=1:size(D,1)
    for j=1:size(D,2)
        H(D(i,j)+1) = H(D(i,j)+1) + 1;
    end
end

% Normalize the histogram to a value of 255
H = H*255/numel(D);

% Generate the CDF of the histogram
new_H = zeros(1,256);
for i=1:256
    new_H(i) = round(sum(H(1:i)));
end

% Generate the new disparities using CDF
new_D = D;
for i=1:size(new_D,1)
    for j=1:size(new_D,2)
        new_D(i,j) = new_H(D(i,j)+1);
    end
end

% Display the image with colormap
imshow(new_D)
colormap(gca,jet)
colorbar

% Save the image
saveas(gca,['../output/SSD_energy_minima/Disparity Map using Energy Minimization for window_size = ' num2str(window_size) '.jpg'])

end
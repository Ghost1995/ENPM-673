%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code implements the second question of the mid-term.
% 
% Submitted by: Ashwin Goyal (UID - 115526297)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read the example image
I = imread('../input/low-contrast-ex.png');

% Generate a histogram of the intensities in the image
H = zeros(1,256);
for i=1:size(I,1)
    for j=1:size(I,2)
        H(I(i,j)+1) = H(I(i,j)+1) + 1;
    end
end

% Normalize the histogram to a value of 255
H = H*255/numel(I);

% Generate the CDF of the histogram
new_H = zeros(1,256);
for i=1:256
    new_H(i) = round(sum(H(1:i)));
end

% Generate the new image using CDF
new_I = I;
for i=1:size(new_I,1)
    for j=1:size(new_I,2)
        new_I(i,j) = new_H(I(i,j)+1);
    end
end

% Save the image
imwrite(new_I,'../output/improved-contrast.png');
imshow(new_I)
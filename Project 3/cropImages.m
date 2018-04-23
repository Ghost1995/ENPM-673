%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code extracts regions of interest in an image
% 
% Input:
%    imageFolder --> Location of the images to be cropped
%   outputFolder --> Location where the cropped images need to be saved
% 
% Submitted by: Ashwin Goyal (UID - 115526297)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cropImages(imageFolder, outputFolder)

    % Read image names
    imgFiles = dir([imageFolder '*.jpg']);
    
    % Crop images using roipoly
    for i = 1:length(imgFiles)
        I = imread([imageFolder imgFiles(i).name]);
        I_y = roipoly(I);
        if ~isempty(I_y)
            imwrite(I_y,[outputFolder 'Y_' imgFiles(i).name]);
        end
        I_r = roipoly(I);
        if ~isempty(I_r)
            imwrite(I_r,[outputFolder 'R_' imgFiles(i).name]);
        end
        I_g = roipoly(I);
        if ~isempty(I_g)
            imwrite(I_g,[outputFolder 'G_' imgFiles(i).name]);
        end
    end

end
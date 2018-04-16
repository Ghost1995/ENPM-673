%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code extracts regions of interest in an image
% 
% Input:
%   imageFolder --> Location of the images to get ROI
% 
% Submitted by: Ashwin Goyal (UID - 115526297)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cropImages

    imageFolder = '.\P1_submission\ColorSeg\Images\TrainingSet\Frames\';
    imgFiles = dir([imageFolder '*.jpg']);
    for i = 1:length(imgFiles)
        I = imread([imageFolder imgFiles(i).name]);
        I_y = roipoly(I);
        if ~isempty(I_y)
            imwrite(I_y,['.\P1_submission\ColorSeg\Images\TrainingSet\CroppedBuoys\Y_' imgFiles(i).name]);
        end
        I_r = roipoly(I);
        if ~isempty(I_r)
            imwrite(I_r,['.\P1_submission\ColorSeg\Images\TrainingSet\CroppedBuoys\R_' imgFiles(i).name]);
        end
        I_g = roipoly(I);
        if ~isempty(I_g)
            imwrite(I_g,['.\P1_submission\ColorSeg\Images\TrainingSet\CroppedBuoys\G_' imgFiles(i).name]);
        end
    end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code converts a video into images and store the images in two sets,
% namely Training Set and Testing Set.
% 
% Input:
%    vidFilename --> Filename of the video to be converted into images
%   outputFolder --> Location where the extracted images need to be saved
% 
% Submitted by: Ashwin Goyal (UID - 115526297)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function video2images(vidFilename, outputFolder)

    % Initialize a video object
    vidObj = VideoReader(vidFilename);
    
    % Extract frames
    i = 1;
    while hasFrame(vidObj)
        vidFrame = readFrame(vidObj);
        if rem(i,10) == 1
            filename = [outputFolder{1} num2str(i,'%03d') '.jpg'];
        else
            filename = [outputFolder{2} num2str(i,'%03d') '.jpg'];
        end
        imwrite(vidFrame,filename);
        i = i+1;
    end

end
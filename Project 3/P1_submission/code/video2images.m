%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code converts a video into images at a given frame rate
% 
% Input:
%    vidFilename --> Filename of the video to be converted into images
%      frameRate --> Frame rate to be used to extract images
%   outputFolder --> Location where the extracted images need to be saved
% 
% Submitted by: Ashwin Goyal (UID - 115526297)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function video2images(vidFilename, frameRate, outputFolder)

    % Initialize a video object
    vidObj = VideoReader(vidFilename);
    if ~isempty(frameRate)
        vidObj.FrameRate = frameRate;
    end
    
    % Extract frames
    i = 1;
    while hasFrame(vidObj)
        vidFrame = readFrame(vidObj);
        filename = [outputFolder num2str(i,'%03d') '.jpg'];
        imwrite(vidFrame,filename);
        i = i+1;
    end

end
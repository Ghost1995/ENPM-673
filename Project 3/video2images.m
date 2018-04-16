%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code converts a video into images
% 
% Input:
%   vidFilename --> Filename of the video to be converted into images
% 
% Submitted by: Ashwin Goyal (UID - 115526297)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function video2images

    vidFilename = 'detectbuoy.avi';
    vidObj = VideoReader(vidFilename);
    i = 1;
    while hasFrame(vidObj)
        vidFrame = readFrame(vidObj);
        if i < 10
            filename = ['00' num2str(i) '.jpg'];
        elseif i < 100
            filename = ['0' num2str(i) '.jpg'];
        else
            filename = [num2str(i) '.jpg'];
        end
        if (rem(floor(i/10),2) == 0)
            if rem(i,10) == 0
                filename = ['.\P1_submission\ColorSeg\Images\TestSet\Frames\' filename];
            else
                filename = ['.\P1_submission\ColorSeg\Images\TrainingSet\Frames\' filename];
            end
        else
            if rem(i,10) == 0
                filename = ['.\P1_submission\ColorSeg\Images\TrainingSet\Frames\' filename];
            else
                filename = ['.\P1_submission\ColorSeg\Images\TestSet\Frames\' filename];
            end
        end
        imwrite(vidFrame,filename);
        i = i+1;
    end

end
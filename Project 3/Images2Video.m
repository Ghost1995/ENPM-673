vidObj = VideoReader('detectbuoy.avi');
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


% outputVideo = VideoWriter('ImageData.avi');
% outputVideo.FrameRate = 30;
% open(outputVideo)
% imgFolder = '.\P2_Submission\VisualOdometry\input\stereo\centre\';
% modelFolder = '.\P2_Submission\VisualOdometry\input\model\';
% % Get intrinsic parameters
% [~, ~, ~, ~, ~, LUT] = ReadCameraModel(imgFolder,modelFolder);
% % Get image file names
% imgFiles = dir([imgFolder '*.png']);
% for i = 1:length(imgFiles)
%     % Read the image frame
%     I = histeq(rgb2gray(UndistortImage(demosaic(imread([imgFolder imgFiles(i).name]), 'gbrg'), LUT)));
%     writeVideo(outputVideo,I)
% end
% close(outputVideo)
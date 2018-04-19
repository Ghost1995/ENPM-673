%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code runs all the commands used in this Project.
% 
% Submitted by: Ashwin Goyal (UID - 115526297)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define the folder of training set
trainFolder = '..\input\Images\TrainingSet\Frames\';
% Read all training image names
trainFiles = dir([trainFolder '*.jpg']);

% Define the folder of testing set
testFolder = '..\input\Images\TestSet\Frames\';
% Read all training image names
testFiles = dir([testFolder '*.jpg']);

% Define the folder of cropped buoys
cropFolder = '..\input\Images\TrainingSet\CroppedBuoys\';

% % Extract frames from the given video
% video2images('..\input\detectbuoy.avi',{trainFolder,testFolder})

% % Crop the images to get the training set
% cropImages(trainFolder,cropFolder)

% % Compute average histogram as well as the color distribution
% averageHistogram(trainFolder,cropFolder,'HSV')

% Create video of segmented images using 1-D gaussian
% vidObj = VideoWriter('..\output\segment1D.mp4','MPEG-4');
% vidObj.Quality = 100;
% open(vidObj)
count = 1;
for i = 1:length(testFiles)+length(trainFiles)
    if i == 1
        I = segment1D('RGB',[trainFolder trainFiles(1).name],false);
    elseif rem(i,10) == 1
        count = count + 1;
        I = segment1D('RGB',[trainFolder trainFiles(count).name],false);
    else
        I = segment1D('RGB',[testFolder testFiles(i-count).name],false);
    end
%     for j = 1:6
%         writeVideo(vidObj,I)
%     end
end
% close(vidObj)
        
% I = segment1D('RGB',[testFolder testFiles(35).name],false);



% Define the folder of training set
trainFolder = '..\..\Images\TrainingSet\Frames\';
% Read all training image names
trainFiles = dir([trainFolder '*.jpg']);

% Define the folder of testing set
testFolder = '..\..\Images\TestSet\Frames\';
% Read all training image names
testFiles = dir([testFolder '*.jpg']);

% Segment testing data
for i = [8,15,18,20,24,25,27,28,30,34,36,37,38,39,86,89,91,93,101,106,107,111,112,113,128,129,130,131,132,133,134,141,154,155,156,157,158,160,168,173]%1:length(testFiles)
    segment1D_2('RGB',[testFolder testFiles(i).name],false);
end

% Segment training data
for i = 4%1:length(trainFiles)
    segment1D_2('RGB',[trainFolder trainFiles(i).name],false);
end

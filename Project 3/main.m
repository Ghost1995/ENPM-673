%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code runs all the commands used in this Project.
% 
% Submitted by: Ashwin Goyal (UID - 115526297)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Add all required paths
addpath('.\ColorSeg\Scripts\Part0','.\ColorSeg\Scripts\Part1','.\ColorSeg\Scripts\Part2')

% Define the folder of training set
trainFolder = '.\ColorSeg\Images\TrainingSet\Frames\';
% Read all training image names
trainFiles = dir([trainFolder '*.jpg']);

% Define the folder of testing set
testFolder = '.\ColorSeg\Images\TestSet\Frames\';
% Read all training image names
testFiles = dir([testFolder '*.jpg']);

% Define the folder of cropped buoys
cropFolder = '.\ColorSeg\Images\TrainingSet\CroppedBuoys\';

% Extract frames from the given video
% video2images('detectbuoy.avi',{trainFolder,testFolder})

% Crop the images to get the training set
% cropImages(trainFolder,cropFolder)

% Compute average histogram as well as the color distribution
averageHistogram(trainFolder,cropFolder,'RGB')

% Get color distributions
greenDist = []; redDist = []; yellowDist = [];
load('.\ColorSeg\Output\Part0\colorDistributions_RGB.mat','greenDist','redDist','yellowDist')
gmObj_green = gmdistribution(mean(greenDist(:,2)),var(greenDist(:,2)));
gmObj_red = gmdistribution(mean(redDist(:,1)),var(redDist(:,1)));
gmObj_yellow = gmdistribution(mean(mean(yellowDist(:,1:2),2)),var(mean(yellowDist(:,1:2),2)));
gmObj = {gmObj_green; gmObj_red; gmObj_yellow};
% Plot the three gaussians
figure('units','normalized','outerposition',[0 0 1 1])
plot(0:255,gauss(gmObj{1},(0:255)'))
title('1-D Gaussian to Detect Green Buoy')
xlabel('Intensity')
ylabel('Probability')
saveas(gcf,'.\ColorSeg\Output\Part0\G_gauss1D.jpg')
plot(0:255,gauss(gmObj{2},(0:255)'))
title('1-D Gaussian to Detect Red Buoy')
xlabel('Intensity')
ylabel('Probability')
saveas(gcf,'.\ColorSeg\Output\Part0\R_gauss1D.jpg')
plot(0:255,gauss(gmObj{3},(0:255)'))
title('1-D Gaussian to Detect Yellow Buoy')
xlabel('Intensity')
ylabel('Probability')
saveas(gcf,'.\ColorSeg\Output\Part0\Y_gauss1D.jpg')
% Create video of segmented images using 1-D gaussian
vidObj = VideoWriter('.\ColorSeg\Output\Part0\segment1D.mp4','MPEG-4');
vidObj.Quality = 100;
open(vidObj)
count = 0;
for i = 1:length(testFiles)+length(trainFiles)
    if rem(i,10) == 1
        count = count + 1;
        I = segment1D(gmObj,[trainFolder trainFiles(count).name]);
    else
        I = segment1D(gmObj,[testFolder testFiles(i-count).name]);
    end
    for j = 1:6
        writeVideo(vidObj,I)
    end
end
close(vidObj)

% Create data using 3 1-D gaussians
data = cat(3,linspace(10,30)',linspace(30,50)',linspace(50,70)');
mu = [mean(data(:,:,1));mean(data(:,:,2));mean(data(:,:,3))];
sigma = cat(3,var(data(:,:,1)),var(data(:,:,2)),var(data(:,:,3)));
X = [data(:,:,1); data(:,:,2); data(:,:,3)];
X = sort(X);
figure('units','normalized','outerposition',[0 0 1 1])
gmObj = gmdistribution(mu,sigma);
plot(X,gauss(gmObj,X))
hold on
% Use EM to retrieve the three gaussians used
[gmObj_1D3N,isConverged] = EM(X,3);
if isConverged
    plot(X,gauss(gmObj_1D3N,X))
    xlabel('Data Points')
    ylabel('Probability')
    title('Probability Distribution')
    legend('Actual PDF','Derived PDF')
    saveas(gcf,'.\ColorSeg\Output\Part1\EM1D3N.jpg')
end
hold off

% Plot the data generated using 3 1-D gaussians again
plot(X,gauss(gmObj,X))
hold on
% Use EM to retrieve four gaussians instead of three
[gmObj_1D4N,isConverged] = EM(X,4);
if isConverged
    plot(X,gauss(gmObj_1D4N,X))
    xlabel('Data Points')
    ylabel('Probability')
    title('Probability Distribution')
    legend('Actual PDF','Derived PDF')
    saveas(gcf,'.\ColorSeg\Output\Part1\EM1D4N.jpg')
end
hold off

% Generate 1-D Gaussian for each buoy
% colorModels('RGB','..\output\ColorModels\',1,5,1);
% Generate 2-D Gaussian for each buoy
% colorModels('RGB','..\output\ColorModels\',1,5,2);

% Generate Gaussian model to be used to model each buoy
[gmObj_green,isConverged] = EM(greenDist(:,1),2);
while ~isConverged
    [gmObj_green,isConverged] = EM(greenDist(:,1),2);
end
figure('units','normalized','outerposition',[0 0 1 1])
plot(0:255,gauss(gmObj_green,(0:255)'))
xlabel('Data Points')
ylabel('Probability')
title('Probability Distribution used for Green Buoy')
saveas(gcf,'.\ColorSeg\Output\Part2\EM_G.jpg')
[gmObj_red,isConverged] = EM(redDist(:,1),1);
while ~isConverged
    [gmObj_red,isConverged] = EM(redDist(:,1),1);
end
plot(0:255,gauss(gmObj_red,(0:255)'))
xlabel('Data Points')
ylabel('Probability')
title('Probability Distribution used for Red Buoy')
saveas(gcf,'.\ColorSeg\Output\Part2\EM_R.jpg')
[gmObj_yellow,isConverged] = EM(yellowDist(:,1),1);
while ~isConverged
    [gmObj_yellow,isConverged] = EM(yellowDist(:,1),1);
end
plot(0:255,gauss(gmObj_yellow,(0:255)'))
xlabel('Data Points')
ylabel('Probability')
title('Probability Distribution used for Yellow Buoy')
saveas(gcf,'.\ColorSeg\Output\Part2\EM_Y.jpg')
gmObjs = {gmObj_green; gmObj_red; gmObj_yellow};
% Create video of segmented images using gaussians generated from EM
vidObj = VideoWriter('.\ColorSeg\Output\Part3\Video\segmentEM.mp4','MPEG-4');
vidObj.Quality = 100;
open(vidObj)
count = 0;
for i = 1:length(testFiles)+length(trainFiles)
    if rem(i,10) == 1
        count = count + 1;
        I = detectBuoy(gmObjs,[trainFolder trainFiles(count).name]);
    else
        I = detectBuoy(gmObjs,[testFolder testFiles(i-count).name]);
    end
    for j = 1:6
        writeVideo(vidObj,I)
    end
end
close(vidObj)

function N = gauss(gmObj, X)
% This function computes N(x|mu,sigma) for N-D Gaussian

    mean = gmObj.mu;
    sigma = gmObj.Sigma;
    mixtureCoeff = gmObj.ComponentProportion;
    N = 0;
    for i = 1:length(mixtureCoeff)
        N = N + mixtureCoeff(i)*(1/(2*pi)^(size(X,2)/2))*(1/sqrt(det(sigma(:,:,i))))*exp(sum(-0.5*((X - mean(i,:))/sigma(:,:,i)).*(X - mean(i,:)),2));
    end

end
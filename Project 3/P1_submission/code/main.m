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

% % Create video of segmented images using 1-D gaussian
% vidObj = VideoWriter('..\output\segment1D.mp4','MPEG-4');
% vidObj.Quality = 100;
% open(vidObj)
% count = 1;
% for i = 1:length(testFiles)+length(trainFiles)
%     if i == 1
%         I = segment1D('RGB',[trainFolder trainFiles(1).name],true,false);
%     elseif rem(i,10) == 1
%         count = count + 1;
%         I = segment1D('RGB',[trainFolder trainFiles(count).name],false,false);
%     else
%         I = segment1D('RGB',[testFolder testFiles(i-count).name],false,false);
%     end
%     for j = 1:6
%         writeVideo(vidObj,I)
%     end
% end
% close(vidObj)

% Create data using 3 1-D gaussians
data = cat(3,linspace(10,30)',linspace(30,50)',linspace(50,70)');
mu = [mean(data(:,:,1));mean(data(:,:,2));mean(data(:,:,3))];
sigma = [std(data(:,:,1));std(data(:,:,2));std(data(:,:,3))];
X = [data(:,:,1); data(:,:,2); data(:,:,3)];
X = sort(X);
figure('units','normalized','outerposition',[0 0 1 1])
Y = (normpdf(X,mu(1),sqrt(sigma(1))) + normpdf(X,mu(2),sqrt(sigma(2))) + normpdf(X,mu(3),sqrt(sigma(3))))/3;
plot(X,Y)
hold on
% Use EM to retrieve the three gaussians used
[new_mu,new_sigma,new_mixtureCoeff] = EM(3,X);
for i = 1:300
    for j = 1:3
        postProb(i,j) = new_mixtureCoeff(j)*gauss(X(i,:),new_mu(j,:),new_sigma(j,:,:));
    end
end
plot(X,sum(postProb,2))
hold off
xlabel('Data Points')
ylabel('Probability')
title('Probability Distribution')
legend('Actual PDF','Derived PDF')
saveas(gcf,'..\output\EM1D3N.jpg')

% Plot the data generated using 3 1-D gaussians again
figure('units','normalized','outerposition',[0 0 1 1])
plot(X,Y)
hold on
% Use EM to retrieve four gaussians instead of three
[new_mu,new_sigma,new_mixtureCoeff] = EM(4,X);
for i = 1:300
    for j = 1:4
        postProb(i,j) = new_mixtureCoeff(j)*gauss(X(i,:),new_mu(j,:),new_sigma(j,:,:));
    end
end
plot(X,sum(postProb,2))
hold off
xlabel('Data Points')
ylabel('Probability')
title('Probability Distribution')
legend('Actual PDF','Derived PDF')
saveas(gcf,'..\output\EM1D4N.jpg')

% Check if the algorithm works for 3-D gaussian
data = cat(3,(1 + rand(100,3))*10,(1 + rand(100,3))*20,(1 + rand(100,3))*30);
mu = [mean(data(:,:,1));mean(data(:,:,2));mean(data(:,:,3))];
sigma = cat(3,cov(data(:,:,1)),cov(data(:,:,2)),cov(data(:,:,3)));
X = [data(:,:,1); data(:,:,2); data(:,:,3)];
X = sort(X);
mixtureCoeff = ones(3,1)/3;
for i = 1:300
    for j = 1:3
        postProb(i,j) = mixtureCoeff(j)*gauss(X(i,:),mu(j,:),sigma(j,:,:));
    end
end
figure('units','normalized','outerposition',[0 0 1 1])
scatter3(X(:,1),X(:,2),X(:,3),40,sum(postProb,2),'filled')
hold on
[new_mu,new_sigma,new_mixtureCoeff] = EM(3,X);
for i = 1:300
    for j = 1:3
        postProb(i,j) = new_mixtureCoeff(j)*gauss(X(i,:),new_mu(j,:),new_sigma(j,:,:));
    end
end
scatter3(X(:,1),X(:,2),X(:,3),40,sum(postProb,2),'filled')
hold off
xlabel('X - Data Points')
ylabel('Y - Data Points')
zlabel('Z - Data Points')
cb = colorbar;
cb.Label.String = 'Probability';
title('Probability Distribution')
legend('Actual PDF','Derived PDF')
saveas(gcf,'..\output\EM3D3N.jpg')




function N = gauss(x, mu, sigma)
% This function computes N(x|mu,sigma)

    sigma = reshape(sigma,[size(sigma,2) size(sigma,3)]);
    N = (1/(2*pi)^(size(x,2)/2))*(1/sqrt(det(sigma)))*exp(-0.5*((x - mu)/sigma)*(x - mu)');
    if isnan(N)
        N = 0;
    end
    
end
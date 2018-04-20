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
% averageHistogram(trainFolder,cropFolder,'RGB')

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

% % Create data using 3 1-D gaussians
% data = cat(3,linspace(10,30)',linspace(30,50)',linspace(50,70)');
% mu = [mean(data(:,:,1));mean(data(:,:,2));mean(data(:,:,3))];
% sigma = cat(3,var(data(:,:,1)),var(data(:,:,2)),var(data(:,:,3)));
% X = [data(:,:,1); data(:,:,2); data(:,:,3)];
% X = sort(X);
% figure('units','normalized','outerposition',[0 0 1 1])
% gmObj = gmdistribution(mu,sigma);
% Y = pdf(gmObj,X);
% plot(X,Y)
% hold on
% % Use EM to retrieve the three gaussians used
% gmObj_1D3N = EM(X,3);
% Y_1D3N = pdf(gmObj_1D3N,X);
% plot(X,Y_1D3N)
% hold off
% xlabel('Data Points')
% ylabel('Probability')
% title('Probability Distribution')
% legend('Actual PDF','Derived PDF')
% saveas(gcf,'..\output\EM1D3N.jpg')
% 
% % Plot the data generated using 3 1-D gaussians again
% figure('units','normalized','outerposition',[0 0 1 1])
% plot(X,Y)
% hold on
% % Use EM to retrieve four gaussians instead of three
% gmObj_1D4N = EM(X,4);
% Y_1D4N = pdf(gmObj_1D4N,X);
% plot(X,Y_1D4N)
% hold off
% xlabel('Data Points')
% ylabel('Probability')
% title('Probability Distribution')
% legend('Actual PDF','Derived PDF')
% saveas(gcf,'..\output\EM1D4N.jpg')


% title('Contour lines of pdf');
% 
% 
% Generate 1000 random variates from the GMM.
% 
% rng('default'); % For reproducibility
% X = random(gm,1000);
% Plot the variates with the pdf contours.
% 
% hold on
% scatter(X(:,1),X(:,2),10,'.') % Scatter plot with points of size 10
% title('Contour lines of pdf and Simulated Data');
% 




% for N = 5
% figure('units','normalized','outerposition',[0 0 1 1])
% bar(0:255,redHist(:,1),'r');
% hold on
% redDist_red = double(sort(redDist(:,1)));
% redObj_red = EM(redDist_red,N);
% redY_red = pdf(redObj_red,redDist_red);
% plot(redDist_red,redY_red)
% xlabel('Intensity')
% ylabel('Frequency')
% title(['Red Color Band of Green Buoy represented by ' num2str(N) ' Gaussians'])
% saveas(gcf,['..\output\redGauss_red_' num2str(N) '.jpg']);
% 
% figure('units','normalized','outerposition',[0 0 1 1])
% bar(0:255,redHist(:,2),'g');
% hold on
% redDist_green = double(sort(redDist(:,2)));
% redObj_green = EM(redDist_green,N);
% redY_green = pdf(redObj_green,redDist_green);
% plot(redDist_green,redY_green)
% xlabel('Intensity')
% ylabel('Frequency')
% title(['Green Color Band of Green Buoy represented by ' num2str(N) ' Gaussians'])
% saveas(gcf,['..\output\redGauss_green_' num2str(N) '.jpg']);
% 
% figure('units','normalized','outerposition',[0 0 1 1])
% bar(0:255,redHist(:,3),'b');
% hold on
% redDist_blue = double(sort(redDist(:,3)));
% redObj_blue = EM(redDist_blue,N);
% redY_blue = pdf(redObj_blue,redDist_blue);
% plot(redDist_blue,redY_blue,'r')
% xlabel('Intensity')
% ylabel('Frequency')
% title(['Blue Color Band of Green Buoy represented by ' num2str(N) ' Gaussians'])
% saveas(gcf,['..\output\redGauss_blue_' num2str(N) '.jpg']);
% end





% Check if the algorithm works for 3-D gaussian
data = cat(3,(1 + rand(100,2))*125,(1 + rand(100,2))*255,(1 + rand(100,2))*200);
mu = [mean(data(:,:,1));mean(data(:,:,2));mean(data(:,:,3))];
sigma = cat(3,cov(data(:,:,1)),cov(data(:,:,2)),cov(data(:,:,3)));
X = [data(:,:,1); data(:,:,2); data(:,:,3)];
X = sort(X);
gmObj = gmdistribution(mu,sigma);
gmPDF = @(x,y)pdf(gmObj,[x y]);
ezcontour(gmPDF,[0 500],[0 500]);
Xtemp = random(gmObj,300);
hold on
scatter(X(:,1),X(:,2),10,'.')

hold off

% Y = pdf(gmObj,X);
% figure('units','normalized','outerposition',[0 0 1 1])
% plot3(X(:,1),X(:,2),Y)
% hold on
% % subplot(1,2,3)
% % plot(X(:,3),Y)
% % hold on
% % 
% % fimplicit3(@(x,y,z)pdf(gmObj,[x, y, z]),[0 10 0 20 0 30])
% % scatter3(X(:,1),X(:,2),X(:,3),40,sum(postProb,2),'filled')
gmObj_2D3N = EM(X,3);
gmPDF = @(x,y)pdf(gmObj_2D3N,[x y]);
ezcontour(gmPDF,[0 500],[0 500]);

% Y_2D3N = pdf(gmObj_2D3N,X);
% plot3(X(:,1),X(:,2),Y_2D3N)
% % hold on
% % subplot(1,2,2)
% % plot3(X(:,2),X(:,3),Y_3D3N)
% % hold on
% % subplot(2,2,3)
% % plot(X(:,3),Y_3D3N)
% % hold on
gmm = fitgmdist(X,3,'Options',statset('MaxIter',1500));
gmPDF = @(x,y)pdf(gmm,[x y]);
fcontour(gmPDF,[0 255 0 255]);
% Ygmm = pdf(gmm,X);
% plot3(X(:,1),X(:,2),Ygmm)
% % subplot(1,2,2)
% % plot3(X(:,2),X(:,3),Ygmm)
% % subplot(2,2,3)
% % plot(X(:,3),Ygmm)
% plot(X,Ygmm)
% scatter3(X(:,1),X(:,2),X(:,3),40,sum(postProb,2),'filled')
% hold off
% xlabel('X - Data Points')
% ylabel('Y - Data Points')
% zlabel('Z - Data Points')
% cb = colorbar;
% cb.Label.String = 'Probability';
% title('Probability Distribution')
% legend('Actual PDF','Derived PDF')
% % saveas(gcf,'..\output\EM3D3N.jpg')
% gmm = fitgmdist(X,3);




function N = gauss(x, mu, sigma)
% This function computes N(x|mu,sigma)

    sigma = reshape(sigma,[size(sigma,2) size(sigma,3)]);
    N = (1/(2*pi)^(size(x,2)/2))*(1/sqrt(det(sigma)))*exp(-0.5*((x - mu)/sigma)*(x - mu)');
    if isnan(N)
        N = 0;
    end
    
end
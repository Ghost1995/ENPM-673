%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code computes average histogram for individual colors
% 
% Input:
%   colorSpace --> Color space to be used
%        frame --> Location of the images of the buoy
%   plot_gauss --> States whether to plot the gaussian used or not
% 
% Submitted by: Ashwin Goyal (UID - 115526297)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function segment1D_2(colorSpace, frame, plot_gauss)

    % Get color distributions
    try
        load(['..\..\..\output\colorDistributions_' colorSpace '.mat'],'greenDist','redDist','yellowDist')
    catch
%         averageHistogram('..\input\Images\TrainingSet\Frames\','..\input\Images\TrainingSet\CroppedBuoys\','RGB')
        averageHistogram('RGB')
        load(['..\..\..\output\colorDistributions_' colorSpace '.mat'],'greenDist','redDist','yellowDist')
    end
    
    % Generate 1-D gaussian for green buoy
    [greenMean,greenSigma] = normfit(greenDist(:,2));
    % Generate 1-D gaussian for red buoy
    [redMean,redSigma] = normfit(redDist(:,1));
    % Generate 1-D gaussian for yellow buoy
    [yellowMean,yellowSigma] = normfit(mean(yellowDist(:,1:2),2));
    
    % Plot the three gaussians if asked
    if plot_gauss
        plot(0:255,normpdf(0:255,greenMean,greenSigma))
        title('1-D Gaussian to Detect Green Buoy')
        xlabel('Intensity')
        ylabel('Probability')
        saveas(gcf,'../../Output/Part0/G_gauss1D.jpg')
        
        plot(0:255,normpdf(0:255,redMean,redSigma))
        title('1-D Gaussian to Detect Red Buoy')
        xlabel('Intensity')
        ylabel('Probability')
        saveas(gcf,'../../Output/Part0/R_gauss1D.jpg')
        
        plot(0:255,normpdf(0:255,yellowMean,yellowSigma))
        title('1-D Gaussian to Detect Yellow Buoy')
        xlabel('Intensity')
        ylabel('Probability')
        saveas(gcf,'../../Output/Part0/Y_gauss1D.jpg')
    end
    
    % Read the image
    I = imread(frame);
    imshow(I) %%%%%%%%%%%%%% to be removed
    I_green = double(I(:,:,2));
    I_red = double(I(:,:,1));
    I_yellow = mean(double(I(:,:,1:2)),3);
    
    % Compute gaussian probabilities
    greenProb = zeros(size(I_green));
    redProb = zeros(size(I_red));
    yellowProb = zeros(size(I_yellow));
    for i = 1:size(I,1)
        for j = 1:size(I,2)
            greenProb(i,j) = gauss(I_green(i,j),greenMean,greenSigma);
            redProb(i,j) = gauss(I_red(i,j),redMean,redSigma);
            yellowProb(i,j) = gauss(I_yellow(i,j),yellowMean,yellowSigma);
        end
    end
    
    % Identify green buoy
    greenBuoy = greenProb > std2(greenProb);
    greenBuoy = bwareafilt(bwmorph(imfill(bwmorph(bwmorph(greenBuoy,'thicken',10),'close',5),'holes'),'thin',7),[450 625]);
    greenProperty = regionprops(greenBuoy);
    greenArea = [];
    for i = 1:length(greenProperty)
        if (greenProperty(i).Centroid(2) > 150)&&(greenProperty(i).Centroid(2) < size(greenProb,1)-100)
            greenArea = [greenArea; greenProperty(i).Area];
        end
    end
    if ~isempty(greenArea)
        [greenArea,greenInd] = sort(greenArea,'descend');
        greenIndex = 1;
        greenExist = true;
    else
        greenIndex = 0;
        greenExist = false;
    end
    
    % Identify red buoy
    redBuoy = redProb > std2(redProb);
    redBuoy = bwareafilt(bwmorph(imfill(bwmorph(bwmorph(redBuoy,'thicken',10),'close',5),'holes'),'thin',8),[400 5000]);
    redProperty = regionprops(redBuoy);
    redArea = [];
    for i = 1:length(redProperty)
        if (redProperty(i).Centroid(2) > 150)&&(redProperty(i).Centroid(2) < size(redProb,1)-100)
            redArea = [redArea; redProperty(i).Area];
        end
    end
    if ~isempty(redArea)
        [redArea,redInd] = sort(redArea,'descend');
        redIndex = 1;
        redExist = true;
    else
        redIndex = 0;
        redExist = false;
    end
    
    % Identify yellow buoy
    yellowBuoy = yellowProb > std2(yellowProb);
    yellowBuoy = bwareafilt(bwmorph(imfill(bwmorph(bwmorph(yellowBuoy,'thicken',10),'close',5),'holes'),'thin',8),[500 4000]);
    yellowProperty = regionprops(yellowBuoy);
    yellowArea = [];
    for i = 1:length(yellowProperty)
        if (yellowProperty(i).Centroid(2) > 150)&&(yellowProperty(i).Centroid(2) < size(yellowProb,1)-100)
            yellowArea = [yellowArea; yellowProperty(i).Area];
        end
    end
    if ~isempty(yellowArea)
        [yellowArea,yellowInd] = sort(yellowArea,'descend');
        yellowIndex = 1;
        yellowExist = true;
    else
        yellowIndex = 0;
        yellowExist = false;
    end
    
    % Create a interdependency grid
    if greenExist && redExist
        greenGrid = false(length(greenArea),length(redArea));
        for i = 1:length(greenInd)
            for j = 1:length(redInd)
                if (redProperty(redInd(j)).Centroid(1) < greenProperty(greenInd(i)).Centroid(1))&&...
                        (abs(redProperty(redInd(j)).Area - greenProperty(greenInd(i)).Area) < 1000)&&...
                        (norm(redProperty(redInd(j)).Centroid - greenProperty(greenInd(i)).Centroid) > 25)
                    greenGrid(i,j) = true;
                end
            end
        end
        temp = find(any(greenGrid,2),1);
        if ~isempty(temp)
            greenIndex = temp;
            redIndex = find(greenGrid(greenIndex,:),1);
        else
            greenExist = false;
        end
    end
    if redExist && yellowExist
        redGrid = false(length(redArea),length(yellowArea));
        for i = 1:length(redInd)
            for j = 1:length(yellowInd)
                if (yellowProperty(yellowInd(j)).Centroid(1) < redProperty(redInd(i)).Centroid(1))&&...
                        (abs(yellowProperty(yellowInd(j)).Area - redProperty(redInd(i)).Area) < 1000)&&...
                        (norm(yellowProperty(yellowInd(j)).Centroid - redProperty(redInd(i)).Centroid) > 25)
                    redGrid(i,j) = true;
                end
            end
        end
        temp = find(any(redGrid,2),1);
        if ~isempty(temp)
            if greenExist
                if redIndex == temp
                    yellowIndex = find(redGrid(redIndex,:),1);
                else
                    grid_3D = false(length(greenArea),length(redArea),length(yellowArea));
                    for i = 1:length(greenInd)
                        for j = 1:length(redInd)
                            if (redProperty(redInd(j)).Centroid(1) < greenProperty(greenInd(i)).Centroid(1))&&...
                                    (abs(redProperty(redInd(j)).Area - greenProperty(greenInd(i)).Area) < 1000)&&...
                                    (norm(redProperty(redInd(j)).Centroid - greenProperty(greenInd(i)).Centroid) > 25)
                                for k = 1:length(yellowInd)
                                    if (yellowProperty(yellowInd(k)).Centroid(1) < redProperty(redInd(j)).Centroid(1))&&...
                                            (abs(yellowProperty(yellowInd(k)).Area - redProperty(redInd(j)).Area) < 1000)&&...
                                            (norm(yellowProperty(yellowInd(k)).Centroid - redProperty(redInd(j)).Centroid) > 25)
                                        grid_3D(i,j,k) = true;
                                    end
                                end
                            end
                        end
                    end
                    if isempty(find(any(grid_3D,3),1))
                        greenExist = false;
                        redExist = false;
                        yellowExist = false;
                    else
                        temp = find(any(any(grid_3D,3),2),1);
                        greenIndex = temp;
                        redIndex = find(any(grid_3D(greenIndex,:,:),2),1);
                        yellowIndex = find(grid_3D(greenIndex,redIndex,:),1);
                    end
                end
            else
                redIndex = temp;
                yellowIndex = find(redGrid(redIndex,:),1);
            end
        else
            redExist = false;
        end
    end
    
    % Improve detection
    if redExist
        if length(redArea) > redIndex
            if (redArea(redIndex) - redArea(redIndex+1) > 250)&&(redArea(redIndex) - redArea(redIndex+1) < 1000)
                redNextIndex = redIndex;
                redIndex = redIndex + 1;
            else
                redNextIndex = redIndex + 1;
            end
        end
    end
    if yellowExist
        if length(yellowArea) > yellowIndex
            if (yellowArea(yellowIndex) - yellowArea(yellowIndex+1) > 250)&&(yellowArea(yellowIndex) - yellowArea(yellowIndex+1) < 1000)
                yellowNextIndex = yellowIndex;
                yellowIndex = yellowIndex + 1;
            else
                yellowNextIndex = yellowIndex + 1;
            end
        end
    end
    
    % Verify detected Buoy
    if redExist && yellowExist && ~redGrid(redIndex,yellowIndex)
        if (exist('redNextIndex','var') ~= 0) && (exist('yellowNextIndex','var') ~= 0)
            if redGrid(redNextIndex,yellowIndex) && ~redGrid(redIndex,yellowNextIndex)
                redIndex = redNextIndex;
            elseif ~redGrid(redNextIndex,yellowIndex) && redGrid(redIndex,yellowNextIndex)
                yellowIndex = yellowNextIndex;
            elseif ~redGrid(redNextIndex,yellowIndex) && ~redGrid(redIndex,yellowNextIndex)
                redExist = false;
                yellowExist = false;
            else
                if abs(redArea(redIndex) - yellowArea(yellowNextIndex)) < abs(redArea(redNextIndex) - yellowArea(yellowIndex))
                    yellowIndex = yellowNextIndex;
                else
                    redIndex = redNextIndex;
                end
            end
        elseif (exist('redNextIndex','var') == 0) && (exist('yellowNextIndex','var') ~= 0)
            if length(redArea) > redIndex
                redExist = false;
            else
                yellowIndex = yellowNextIndex;
            end
        elseif (exist('redNextIndex','var') ~= 0) && (exist('yellowNextIndex','var') == 0)
            if length(yellowArea) > yellowIndex
                yellowExist = false;
            else
                redIndex = redNextIndex;
            end
        else
            redExist = false;
            yellowExist = false;
        end
    end
    if greenExist && redExist && ~greenGrid(greenIndex,redIndex)
        clear redNextIndex
        if redExist
            if length(redArea) > redIndex
                if (redArea(redIndex) - redArea(redIndex+1) > 250)&&(redArea(redIndex) - redArea(redIndex+1) < 1000)
                    redNextIndex = redIndex;
                    redIndex = redIndex + 1;
                else
                    redNextIndex = redIndex + 1;
                end
            end
        end
        if exist('redNextIndex','var') ~= 0
            if greenGrid(greenIndex,redNextIndex)
                redIndex = redNextIndex;
            else
                if any(find(greenGrid(:,redIndex)) ~= greenIndex)
                    greenIndex = find(find(greenGrid(:,redIndex)) ~= greenIndex,1);
                elseif any(find(greenGrid(:,redNextIndex)) ~= greenIndex)
                    greenIndex = find(find(greenGrid(:,redNextIndex),1) ~= greenIndex,1);
                    redIndex = redNextIndex;
                else
                    greenExist = false;
                    redExist = false;
                end
            end
        elseif any(find(greenGrid(:,redIndex)) ~= greenIndex)
            greenIndex = find(greenGrid(:,redIndex),1);
        else
            greenExist = false;
            redExist = false;
        end
    end
    
    % Plot green buoy
    if greenExist
        greenConnected = bwconncomp(greenBuoy);
        greenBuoy = zeros(size(greenBuoy));
        greenBuoy(greenConnected.PixelIdxList{greenInd(greenIndex)}) = 1;
        greenBoundary = bwboundaries(greenBuoy);
        greenBoundary = reshape(flip(greenBoundary{1}'),1,numel(greenBoundary{1}));
        I = insertShape(I,'Polygon',greenBoundary,'LineWidth',2,'Color','g');
    end
    
    % Plot red buoy
    if redExist
        redConnected = bwconncomp(redBuoy);
        redBuoy = zeros(size(redBuoy));
        redBuoy(redConnected.PixelIdxList{redInd(redIndex)}) = 1;
        redBoundary = bwboundaries(redBuoy);
        redBoundary = reshape(flip(redBoundary{1}'),1,numel(redBoundary{1}));
        I = insertShape(I,'Polygon',redBoundary,'LineWidth',2,'Color','r');
    end
    
    % Plot yellow buoy
    if yellowExist
        yellowConnected = bwconncomp(yellowBuoy);
        yellowBuoy = zeros(size(yellowBuoy));
        yellowBuoy(yellowConnected.PixelIdxList{yellowInd(yellowIndex)}) = 1;
        yellowBoundary = bwboundaries(yellowBuoy);
        yellowBoundary = reshape(flip(yellowBoundary{1}'),1,numel(yellowBoundary{1}));
        I = insertShape(I,'Polygon',yellowBoundary,'LineWidth',2,'Color','y');
    end
    
    % Plot the image before saving it
    imshow(I)
%     imwrite(I,['../../Output/Part0/seg_' Frame(end-6:end)]);
end

function N = gauss(x, mu, sigma)
% This function computes N(x|mu,sigma)

    N = (1/(2*pi)^(length(x)/2))*(1/sqrt(det(sigma)))*exp(-0.5*((x - mu)/sigma)*(x - mu)');
    if isnan(N)
        N = 0;
    end
    
end
%     % Create 1-D gaussians
%     value = (0.5:255.5)';
%     [muG,sigmaG] = normfit(value,0.05,zeros(size(value)),greenHist);
%     [muR,sigmaR] = normfit(value,0.05,zeros(size(value)),redHist);
%     [muY,sigmaY] = normfit(value,0.05,zeros(size(value)),yellowHist);
%     totalHist = (sum(greenHist)*greenHist + sum(redHist)*redHist + sum(yellowHist)*yellowHist)/sum(greenHist + redHist + yellowHist);
%     [mu,sigma] = normfit(value,0.05,zeros(size(value)),totalHist);
%     
%     % Save gaussians being used
%     normG = normpdf(value,muG,sigmaG);
%     normR = normpdf(value,muR,sigmaR);
%     normY = normpdf(value,muY,sigmaY);
%     norm = normpdf(value,mu,sigma);
%     plot(value,normG);
%     saveas(gcf,'../../Output/Part0/G_gauss1D.jpg')
%     plot(value,normR);
%     saveas(gcf,'../../Output/Part0/R_gauss1D.jpg')
%     plot(value,normY);
%     saveas(gcf,'../../Output/Part0/Y_gauss1D.jpg')
%     plot(value,norm);
%     saveas(gcf,'../../Output/Part0/gauss1D.jpg')
    
    
% end
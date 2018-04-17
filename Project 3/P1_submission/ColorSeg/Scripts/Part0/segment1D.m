%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code computes average histogram for individual colors
% 
% Input:
%   Frame --> Location of the images of the buoy
% 
% Submitted by: Ashwin Goyal (UID - 115526297)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function segment1D(Frame)

%     % Define the folder of cropped buoys
%     imageFolder = '..\..\Images\TrainingSet\Frames\';
%     
%     % Read image names
%     imgFiles = dir([imageFolder '*.jpg']);
    
    % Compute average color histogram for all images
    [greenHist,redHist,yellowHist] = averageHistogram('RGB');
    
    % Generate 1-D gaussian for green buoy
    [greenMean,greenSigma] = normfit(greenHist(:,2));
    plot(0:255,normpdf(0:255,greenMean,greenSigma))
    title('1-D Gaussian to Detect Green Buoy')
    xlabel('Intensity')
    ylabel('Probability')
    saveas(gcf,'../../Output/Part0/G_gauss1D.jpg')
    
    % Generate 1-D gaussian for red buoy
    [redMean,redSigma] = normfit(redHist(:,1));
    plot(0:255,normpdf(0:255,redMean,redSigma))
    title('1-D Gaussian to Detect Red Buoy')
    xlabel('Intensity')
    ylabel('Probability')
    saveas(gcf,'../../Output/Part0/R_gauss1D.jpg')
    
    % Generate 1-D gaussian for yellow buoy
    [yellowMean,yellowSigma] = normfit(mean(yellowHist(:,1:2),2));
    plot(0:255,normpdf(0:255,yellowMean,yellowSigma))
    title('1-D Gaussian to Detect Yellow Buoy')
    xlabel('Intensity')
    ylabel('Probability')
    saveas(gcf,'../../Output/Part0/Y_gauss1D.jpg')
    
    % Read the image
    I = imread(Frame);
    I_double = double(I);
    
    % Compute gaussian probabilities
    greenProb = zeros(size(I_double,1),size(I_double,2));
    redProb = zeros(size(I_double,1),size(I_double,2));
    yellowProb = zeros(size(I_double,1),size(I_double,2));
    for i = 1:size(I_double,1)
        for j = 1:size(I_double,2)
            greenProb(i,j) = gauss(I_double(i,j,2),greenMean,greenSigma);
            redProb(i,j) = gauss(I_double(i,j,1),redMean,redSigma);
            yellowProb(i,j) = gauss(mean(I_double(i,j,1:2)),yellowMean,yellowSigma);
        end
    end
    
    % Identify green buoy
    greenBuoy = greenProb > std2(greenProb);
    greenBuoy = bwareafilt(bwmorph(imfill(bwmorph(bwmorph(greenBuoy,'thicken',10),'close',5),'holes'),'thin',7),[475 625]);
    greenProperty = regionprops(greenBuoy);
    maxArea = 0;
    greenInd = [];
    for i = 1:length(greenProperty)
        if maxArea < greenProperty(i).Area
            if (all(greenProperty(i).Centroid > 30))&&(all(greenProperty(i).Centroid < flip(size(greenProb))-30))
                maxArea = greenProperty(i).Area;
                greenInd = i;
            end
        end
    end
    
    % Identify red buoy
    redBuoy = redProb > std2(redProb);
    redBuoy = bwareafilt(bwmorph(imfill(bwmorph(bwmorph(redBuoy,'thicken',10),'close',5),'holes'),'thin',8),[450 5000]);
    redProperty = regionprops(redBuoy);
    redMaxArea = 0;
    redNextMax = 0;
    redInd = [];
    redNextInd = [];
    for i = 1:length(redProperty)
        if redMaxArea < redProperty(i).Area
            if (all(redProperty(i).Centroid > 30))&&(all(redProperty(i).Centroid < flip(size(redProb))-30))
                if ~isempty(redInd)
                    redNextMax = redMaxArea;
                    redNextInd = redInd;
                end
                redMaxArea = redProperty(i).Area;
                redInd = i;
            end
        elseif redNextMax < redProperty(i).Area
            if (all(redProperty(i).Centroid > 30))&&(all(redProperty(i).Centroid < flip(size(redProb))-30))
                redNextMax = redProperty(i).Area;
                redNextInd = i;
            end
        end
    end
    
    % Identify yellow buoy
    yellowBuoy = yellowProb > std2(yellowProb);
    yellowBuoy = bwareafilt(bwmorph(imfill(bwmorph(bwmorph(yellowBuoy,'thicken',10),'close',5),'holes'),'thin',8),[500 4000]);
    yellowProperty = regionprops(yellowBuoy);
    yellowMaxArea = 0;
    yellowNextMax = 0;
    yellowInd = [];
    yellowNextInd = [];
    for i = 1:length(yellowProperty)
        if yellowMaxArea < yellowProperty(i).Area
            if (all(yellowProperty(i).Centroid > 30))&&(all(yellowProperty(i).Centroid < flip(size(yellowProb))-30))
                if ~isempty(yellowInd)
                    yellowNextMax = yellowMaxArea;
                    yellowNextInd = yellowInd;
                end
                yellowMaxArea = yellowProperty(i).Area;
                yellowInd = i;
            end
        elseif yellowNextMax < yellowProperty(i).Area
            if (all(yellowProperty(i).Centroid > 30))&&(all(yellowProperty(i).Centroid < flip(size(yellowProb))-30))
                yellowNextMax = yellowProperty(i).Area;
                yellowNextInd = i;
            end
        end
    end
    
    % Check for overlap
    if (~isempty(redInd))&&(~isempty(yellowInd))
        if (redMaxArea - redNextMax > 500)&&(redMaxArea - redNextMax < 1000)&&(~isempty(redNextInd))
            tempInd = redNextInd;
            redNextInd = redInd;
            redInd = tempInd;
        end
        if (yellowMaxArea - yellowNextMax > 500)&&(yellowMaxArea - yellowNextMax < 1000)&&(~isempty(yellowNextInd))
            tempInd = yellowNextInd;
            yellowNextInd = yellowInd;
            yellowInd = tempInd;
        end
        if norm(redProperty(redInd).Centroid - yellowProperty(yellowInd).Centroid) < 20
            if (~isempty(redNextInd))&&(~isempty(yellowNextInd))
                if (yellowProperty(yellowInd).Centroid(1) < redProperty(redNextInd).Centroid(1))&&(yellowProperty(yellowNextInd).Centroid(1) > redProperty(redInd).Centroid(1))
                    redInd = redNextInd;
                elseif (yellowProperty(yellowInd).Centroid(1) > redProperty(redNextInd).Centroid(1))&&(yellowProperty(yellowNextInd).Centroid(1) < redProperty(redInd).Centroid(1))
                    yellowInd = yellowNextInd;
                elseif (yellowProperty(yellowInd).Centroid(1) > redProperty(redNextInd).Centroid(1))&&(yellowProperty(yellowNextInd).Centroid(1) > redProperty(redInd).Centroid(1))
                    redInd = [];
                    yellowInd = [];
                else
                    if abs(redProperty(redInd).Area - yellowProperty(yellowNextInd).Area) < abs(redProperty(redNextInd).Area - yellowProperty(yellowInd).Area)
                        yellowInd = yellowNextInd;
                    else
                        redInd = redNextInd;
                    end
                end
            elseif (isempty(redNextInd))&&(length(redProperty)>1)
                redInd = [];
            elseif (isempty(yellowNextInd))&&(length(yellowProperty)>1)
                yellowInd = [];
            elseif (isempty(redNextInd))&&(~isempty(yellowNextInd))
                yellowInd = yellowNextInd;
            elseif (~isempty(redNextInd))&&(isempty(yellowNextInd))
                redInd = redNextInd;
            end
        end
        % Improve the detection
        if (~isempty(redInd))&&(~isempty(yellowInd))&&(abs(redProperty(redInd).Area - yellowProperty(yellowInd).Area) > 1000)
            initRedInd = redInd;
            initYellowInd = yellowInd;
            if ~isempty(redNextInd)
                if abs(redProperty(redNextInd).Area - yellowProperty(yellowInd).Area) > 1000
                    if ~isempty(yellowNextInd)
                        if abs(redProperty(redNextInd).Area - yellowProperty(yellowNextInd).Area) > 1000
                            if length(redProperty) > 2
                                redLimit = false;
                            else
                                redLimit = true;
                            end
                            if length(yellowProperty) > 2
                                yellowLimit = false;
                            else
                                yellowLimit = true;
                            end
                            if ~(redLimit && yellowLimit)
                                redInd = redNextInd;
                                yellowInd = yellowNextInd;
                                while ~(redLimit && yellowLimit)
                                    if ~redLimit
                                        redMaxArea = redNextMax;
                                        redNextMax = 0;
                                        redInd = redNextInd;
                                        redNextInd = [];
                                        for i = 1:length(redProperty)
                                            if (redMaxArea > redProperty(i).Area)&&(redNextMax < redProperty(i).Area)
                                                if (all(redProperty(i).Centroid > 30))&&(all(redProperty(i).Centroid < flip(size(redProb))-30))
                                                    redNextMax = redProperty(i).Area;
                                                    redNextInd = i;
                                                end
                                            end
                                        end
                                        if isempty(redNextInd)
                                            redLimit = true;
                                        else
                                            redInd = redNextInd;
                                            if abs(redProperty(redInd).Area - yellowProperty(yellowInd).Area) < 1000
                                                redLimit = true;
                                                yellowLimit = true;
                                                continue;
                                            end
                                        end
                                    end
                                    if ~yellowLimit
                                        yellowMaxArea = yellowNextMax;
                                        yellowNextMax = 0;
                                        yellowInd = yellowNextInd;
                                        yellowNextInd = [];
                                        for i = 1:length(yellowProperty)
                                            if (yellowMaxArea > yellowProperty(i).Area)&&(yellowNextMax < yellowProperty(i).Area)
                                                if (all(yellowProperty(i).Centroid > 30))&&(all(yellowProperty(i).Centroid < flip(size(yellowProb))-30))
                                                    yellowNextMax = yellowProperty(i).Area;
                                                    yellowNextInd = i;
                                                end
                                            end
                                        end
                                        if isempty(yellowNextInd)
                                            yellowLimit = true;
                                        else
                                            yellowInd = yellowNextInd;
                                            if abs(redProperty(redInd).Area - yellowProperty(yellowInd).Area) < 1000
                                                redLimit = true;
                                                yellowLimit = true;
                                                continue;
                                            end
                                        end
                                    end
                                end
                            else
                                redInd = [];
                                yellowInd = [];
                            end
                        else
                            redInd = redNextInd;
                            yellowInd = yellowNextInd;
                        end
                    elseif length(redProperty) > 2
                        while abs(redProperty(redNextInd).Area - yellowProperty(yellowInd).Area) > 1000
                            redMaxArea = redNextMax;
                            redNextMax = 0;
                            redNextInd = [];
                            for i = 1:length(redProperty)
                                if (redMaxArea > redProperty(i).Area)&&(redNextMax < redProperty(i).Area)
                                    if (all(redProperty(i).Centroid > 30))&&(all(redProperty(i).Centroid < flip(size(redProb))-30))
                                        redNextMax = redProperty(i).Area;
                                        redNextInd = i;
                                    end
                                end
                            end
                            if isempty(redNextInd)
                                yellowInd = [];
                                break;
                            end
                        end
                        redInd = redNextInd;
                    else
                        redInd = [];
                        yellowInd = [];
                    end
                else
                    redInd = redNextInd;
                end
            elseif ~isempty(yellowNextInd)
                if abs(redProperty(redInd).Area - yellowProperty(yellowNextInd).Area) > 1000
                    if length(yellowProperty) > 2
                        while abs(redProperty(redInd).Area - yellowProperty(yellowNextInd).Area) > 1000
                            yellowMaxArea = yellowNextMax;
                            yellowNextMax = 0;
                            yellowNextInd = [];
                            for i = 1:length(yellowProperty)
                                if (yellowMaxArea > yellowProperty(i).Area)&&(yellowNextMax < yellowProperty(i).Area)
                                    if (all(yellowProperty(i).Centroid > 30))&&(all(yellowProperty(i).Centroid < flip(size(yellowProb))-30))
                                        yellowNextMax = yellowProperty(i).Area;
                                        yellowNextInd = i;
                                    end
                                end
                            end
                            if isempty(yellowNextInd)
                                redInd = [];
                                break;
                            end
                        end
                        yellowInd = yellowNextInd;
                    else
                        redInd = [];
                        yellowInd = [];
                    end
                else
                    yellowInd = yellowNextInd;
                end
            else
                redInd = [];
                yellowInd = [];
            end
            % Check for overlap once more
            if (~isempty(redInd))&&(~isempty(yellowInd))&&(norm(redProperty(redInd).Centroid - yellowProperty(yellowInd).Centroid) < 20)
                if (redInd ~= initRedInd)&&(yellowInd ~= initYellowInd)
                    if (yellowProperty(yellowInd).Centroid(1) < redProperty(initRedInd).Centroid(1))&&(yellowProperty(initYellowInd).Centroid(1) > redProperty(redInd).Centroid(1))
                        redInd = initRedInd;
                    elseif (yellowProperty(yellowInd).Centroid(1) > redProperty(initRedInd).Centroid(1))&&(yellowProperty(initYellowInd).Centroid(1) < redProperty(redInd).Centroid(1))
                        yellowInd = initYellowInd;
                    elseif (yellowProperty(yellowInd).Centroid(1) > redProperty(initRedInd).Centroid(1))&&(yellowProperty(initYellowInd).Centroid(1) > redProperty(redInd).Centroid(1))
                        redInd = [];
                        yellowInd = [];
                    else
                        if abs(redProperty(redInd).Area - yellowProperty(initYellowInd).Area) < abs(redProperty(initRedInd).Area - yellowProperty(yellowInd).Area)
                            yellowInd = initYellowInd;
                        else
                            redInd = initRedInd;
                        end
                    end
                elseif redInd ~= initRedInd
                    redInd = initRedInd;
                elseif yellowInd ~= initYellowInd
                    yellowInd = initYellowInd;
                end
            end
        end
    elseif ~isempty(redInd)
        if (redMaxArea - redNextMax > 325)&&(redMaxArea - redNextMax < 1000)&&(~isempty(redNextInd))
            redInd = redNextInd;
        end
    elseif ~isempty(yellowInd)
        if (yellowMaxArea - yellowNextMax > 325)&&(yellowMaxArea - yellowNextMax < 1000)&&(~isempty(yellowNextInd))
            yellowInd = yellowNextInd;
        end
    end
    
    % Verify green buoys existence
    if (~isempty(yellowInd))&&(~isempty(greenInd))&&(abs(greenProperty(greenInd).Area - yellowProperty(yellowInd).Area) > 1000)
        greenInd = [];
    elseif (~isempty(redInd))&&(~isempty(greenInd))&&(abs(redProperty(redInd).Area - greenProperty(greenInd).Area) > 1000)
        greenInd = [];
    end
    
    % Plot green buoy as there are no disparities
    if ~isempty(greenInd)
        greenConnected = bwconncomp(greenBuoy);
        greenBuoy = zeros(size(greenBuoy));
        greenBuoy(greenConnected.PixelIdxList{greenInd}) = 1;
        greenBoundary = bwboundaries(greenBuoy);
        greenBoundary = reshape(flip(greenBoundary{1}'),1,numel(greenBoundary{1}));
        I = insertShape(I,'Polygon',greenBoundary,'LineWidth',2,'Color','g');
    end
    
    % Plot red buoy
    if ~isempty(redInd)
        redConnected = bwconncomp(redBuoy);
        redBuoy = zeros(size(redBuoy));
        redBuoy(redConnected.PixelIdxList{redInd}) = 1;
        redBoundary = bwboundaries(redBuoy);
        redBoundary = reshape(flip(redBoundary{1}'),1,numel(redBoundary{1}));
        I = insertShape(I,'Polygon',redBoundary,'LineWidth',2,'Color','r');
    end
    
    % Plot yellow buoy
    if ~isempty(yellowInd)
        yellowConnected = bwconncomp(yellowBuoy);
        yellowBuoy = zeros(size(yellowBuoy));
        yellowBuoy(yellowConnected.PixelIdxList{yellowInd}) = 1;
        yellowBoundary = bwboundaries(yellowBuoy);
        yellowBoundary = reshape(flip(yellowBoundary{1}'),1,numel(yellowBoundary{1}));
        I = insertShape(I,'Polygon',yellowBoundary,'LineWidth',2,'Color','y');
    end
    
    % Plot the image before saving it
    imshow(I)
    imwrite(I,['../../Output/Part0/seg_' Frame(end-6:end)]);
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
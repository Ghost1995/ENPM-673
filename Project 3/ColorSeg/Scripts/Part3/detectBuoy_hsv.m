%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code computes segmented images using given Gaussian models for HSV
% 
% Input:
%   gmObjs --> Mean and standard deviation for the three buoys
%    frame --> Location of the images of the buoy
% 
% Output:
%   I --> Segmented image
% 
% Submitted by: Ashwin Goyal (UID - 115526297)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function I = detectBuoy_hsv(gmObjs, frame)

    greenGauss = gmObjs{1};
    redGauss = gmObjs{2};
    yellowGauss = gmObjs{3};
    
    % Read the image
    I = imread(frame);
    I_hsv = rgb2hsv(I);
    imgData = double(reshape(I_hsv,numel(I_hsv)/3,3));

    % Compute gaussian probabilities
    greenProb = reshape(gauss(greenGauss,imgData(:,3)),size(I_hsv,1),size(I_hsv,2));
    redProb = reshape(gauss(redGauss,imgData(:,1)),size(I_hsv,1),size(I_hsv,2));
    yellowProb = reshape(gauss(yellowGauss,imgData(:,1)),size(I_hsv,1),size(I_hsv,2));
    
    % Identify green buoy
    greenBuoy = greenProb > 3*std2(greenProb);
    greenBuoy = bwareafilt(bwmorph(imfill(bwmorph(greenBuoy,'close'),'holes'),'thin'),[175 650]);
    greenProperty = regionprops(greenBuoy);
    greenArea = [];
    greenInd = [];
    for i = 1:length(greenProperty)
        if (greenProperty(i).Centroid(2) > 200)&&(greenProperty(i).Centroid(2) < size(greenProb,1)-100)
            greenArea = [greenArea; greenProperty(i).Area];
            greenInd = [greenInd; i];
        end
    end
    if ~isempty(greenArea)
        [greenArea,sequence] = sort(greenArea,'descend');
        greenInd = greenInd(sequence);
        greenIndex = 1;
        greenExist = true;
    else
        greenIndex = 0;
        greenExist = false;
    end
    
    % Identify red buoy
    redBuoy = redProb > std2(redProb);
    redBuoy = bwareafilt(imfill(redBuoy,'holes'),[275 6750]);
    redProperty = regionprops(redBuoy);
    redArea = [];
    redInd = [];
    for i = 1:length(redProperty)
        if (redProperty(i).Centroid(2) > 150)&&(redProperty(i).Centroid(2) < size(redProb,1)-100)
            redArea = [redArea; redProperty(i).Area];
            redInd = [redInd; i];
        end
    end
    if ~isempty(redArea)
        [redArea,sequence] = sort(redArea,'descend');
        redInd = redInd(sequence);
        redIndex = 1;
        redExist = true;
    else
        redIndex = 0;
        redExist = false;
    end
    
    % Identify yellow buoy
    yellowBuoy = yellowProb > std2(yellowProb);
    yellowBuoy = bwareafilt(imfill(yellowBuoy,'holes'),[600 4250]);
    yellowProperty = regionprops(yellowBuoy);
    yellowArea = [];
    yellowInd = [];
    for i = 1:length(yellowProperty)
        if (yellowProperty(i).Centroid(2) > 150)&&(yellowProperty(i).Centroid(2) < size(yellowProb,1)-100)
            yellowArea = [yellowArea; yellowProperty(i).Area];
            yellowInd = [yellowInd; i];
        end
    end
    if ~isempty(yellowArea)
        [yellowArea,sequence] = sort(yellowArea,'descend');
        yellowInd = yellowInd(sequence);
        yellowIndex = 1;
        yellowExist = true;
    else
        yellowIndex = 0;
        yellowExist = false;
    end
    
    % Create an interdependency grid
    if redExist && yellowExist
        if greenExist
            grid_3D = false(length(greenArea),length(redArea),length(yellowArea));
            for i = 1:length(greenInd)
                for j = 1:length(redInd)
                    if (redProperty(redInd(j)).Centroid(1) < greenProperty(greenInd(i)).Centroid(1))&&...
                            (norm(redProperty(redInd(j)).Centroid - greenProperty(greenInd(i)).Centroid) > 25)
                        rg_dist = norm(redProperty(redInd(j)).Centroid - greenProperty(greenInd(i)).Centroid);
                        for k = 1:length(yellowInd)
                            if (yellowProperty(yellowInd(k)).Centroid(1) < redProperty(redInd(j)).Centroid(1))&&...
                                    (norm(yellowProperty(yellowInd(k)).Centroid - redProperty(redInd(j)).Centroid) > 25)
                                yr_dist = norm(yellowProperty(yellowInd(k)).Centroid - redProperty(redInd(j)).Centroid);
                                if abs(rg_dist - yr_dist) < 25
                                    grid_3D(i,j,k) = true;
                                end
                            end
                        end
                    end
                end
            end
            if ~isempty(find(grid_3D,1))
                greenIndex = find(any(any(grid_3D,3),2),1);
                redIndex = find(any(grid_3D(greenIndex,:,:),3),1);
                yellowIndex = find(grid_3D(greenIndex,redIndex,:),1);
            else
                greenExist = false;
            end
        end
        if ~greenExist
            grid2D = false(length(redArea),length(yellowArea));
            for i = 1:length(redInd)
                for j = 1:length(yellowInd)
                    if (yellowProperty(yellowInd(j)).Centroid(1) < redProperty(redInd(i)).Centroid(1))&&...
                            (norm(yellowProperty(yellowInd(j)).Centroid - redProperty(redInd(i)).Centroid) > 25)
                        grid2D(i,j) = true;
                    end
                end
            end
            if ~isempty(find(grid2D,1))
                redIndex = find(any(grid2D,2),1);
                yellowIndex = find(grid2D(redIndex,:),1);
            else
                redExist = false;
            end
        end
    end
    
    % Plot green buoy
    if greenExist
        greenConnected = bwconncomp(greenBuoy);
        greenBuoy = zeros(size(greenBuoy));
        greenBuoy(greenConnected.PixelIdxList{greenInd(greenIndex)}) = 1;
        greenBoundary = bwboundaries(greenBuoy);
        greenBoundary = reshape(flip(greenBoundary{1}'),1,numel(greenBoundary{1}));
        I = insertShape(I,'Polygon',greenBoundary,'LineWidth',3,'Color','g');
    else
        greenBuoy = zeros(size(greenBuoy));
    end
    
    % Plot red buoy
    if redExist
        redConnected = bwconncomp(redBuoy);
        redBuoy = zeros(size(redBuoy));
        redBuoy(redConnected.PixelIdxList{redInd(redIndex)}) = 1;
        redBoundary = bwboundaries(redBuoy);
        redBoundary = reshape(flip(redBoundary{1}'),1,numel(redBoundary{1}));
        I = insertShape(I,'Polygon',redBoundary,'LineWidth',3,'Color','r');
    else
        redBuoy = zeros(size(redBuoy));
    end
    
    % Plot yellow buoy
    if yellowExist
        yellowConnected = bwconncomp(yellowBuoy);
        yellowBuoy = zeros(size(yellowBuoy));
        yellowBuoy(yellowConnected.PixelIdxList{yellowInd(yellowIndex)}) = 1;
        yellowBoundary = bwboundaries(yellowBuoy);
        yellowBoundary = reshape(flip(yellowBoundary{1}'),1,numel(yellowBoundary{1}));
        I = insertShape(I,'Polygon',yellowBoundary,'LineWidth',3,'Color','y');
    else
        yellowBuoy = zeros(size(yellowBuoy));
    end
    
    % Save the image
    imwrite(greenBuoy | redBuoy | yellowBuoy,['.\ColorSeg\Output\Part3\HSV\binary_' frame(end-6:end)]);
    imwrite(I,['.\ColorSeg\Output\Part3\HSV\Frames\out_' frame(end-6:end)]);

end

function N = gauss(gmObj, X)

    mean = gmObj.mu;
    sigma = gmObj.Sigma;
    mixtureCoeff = gmObj.ComponentProportion;
    N = 0;
    for i = 1:length(mixtureCoeff)
        N = N + mixtureCoeff(i)*(1/(2*pi)^(size(X,2)/2))*(1/sqrt(det(sigma(:,:,i))))*exp(sum(-0.5*((X - mean(i,:))/sigma(:,:,i)).*(X - mean(i,:)),2));
    end

end
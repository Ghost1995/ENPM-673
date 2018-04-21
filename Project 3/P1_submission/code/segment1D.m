%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code computes average histogram for individual colors
% 
% Input:
%         var --> Mean and standard deviation for the three buoys
%       frame --> Location of the images of the buoy
%   saveFrame --> States whether to save the segmented frames or not
% 
% Output:
%   I --> Segmented image
% 
% Submitted by: Ashwin Goyal (UID - 115526297)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function I = segment1D(var, frame, saveFrame)

    % Read the image
    I = imread(frame);
    I_green = double(I(:,:,2));
    I_red = double(I(:,:,1));
    I_yellow = mean(double(I(:,:,1:2)),3);
    
    % Compute gaussian probabilities
    greenProb = zeros(size(I_green));
    redProb = zeros(size(I_red));
    yellowProb = zeros(size(I_yellow));
    for i = 1:size(I,1)
        for j = 1:size(I,2)
            greenProb(i,j) = (1/sqrt(2*pi*var(1,2)^2))*exp(-(1/(2*var(1,2)^2))*(I_green(i,j) - var(1,1))^2);
            redProb(i,j) = (1/sqrt(2*pi*var(2,2)^2))*exp(-(1/(2*var(2,2)^2))*(I_red(i,j) - var(2,1))^2);
            yellowProb(i,j) = (1/sqrt(2*pi*var(3,2)^2))*exp(-(1/(2*var(3,2)^2))*(I_yellow(i,j) - var(3,1))^2);
        end
    end
    
    % Identify green buoy
    greenBuoy = greenProb > std2(greenProb);
    greenBuoy = bwareafilt(bwmorph(imfill(bwmorph(bwmorph(greenBuoy,'thicken',10),'close'),'holes'),'thin',8),[150 700]);
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
    redBuoy = bwareafilt(imfill(bwmorph(bwmorph(redBuoy,'clean',5),'close'),'holes'),[250 5500]);
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
    yellowBuoy = yellowProb > 2*std2(yellowProb);
    yellowBuoy = bwareafilt(imfill(bwmorph(bwmorph(yellowBuoy,'clean',5),'close'),'holes'),[400 4000]);
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
                            (abs(redProperty(redInd(j)).Area - greenProperty(greenInd(i)).Area) < 500)&&...
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
    end
    
    % Plot red buoy
    if redExist
        redConnected = bwconncomp(redBuoy);
        redBuoy = zeros(size(redBuoy));
        redBuoy(redConnected.PixelIdxList{redInd(redIndex)}) = 1;
        redBoundary = bwboundaries(redBuoy);
        redBoundary = reshape(flip(redBoundary{1}'),1,numel(redBoundary{1}));
        I = insertShape(I,'Polygon',redBoundary,'LineWidth',3,'Color','r');
    end
    
    % Plot yellow buoy
    if yellowExist
        yellowConnected = bwconncomp(yellowBuoy);
        yellowBuoy = zeros(size(yellowBuoy));
        yellowBuoy(yellowConnected.PixelIdxList{yellowInd(yellowIndex)}) = 1;
        yellowBoundary = bwboundaries(yellowBuoy);
        yellowBoundary = reshape(flip(yellowBoundary{1}'),1,numel(yellowBoundary{1}));
        I = insertShape(I,'Polygon',yellowBoundary,'LineWidth',3,'Color','y');
    end
    
    % Save the image
    if saveFrame
        imwrite(I,['../output/seg_' frame(end-6:end)]);
    end

end
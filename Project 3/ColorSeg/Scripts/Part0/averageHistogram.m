%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code computes average histogram for individual colors and stores
% color distribution as output
% 
% Input:
%   trainFolder --> Location of the training frames
%    cropFolder --> Location of the cropped buoys
%    colorSpace --> Color space to be used
% 
% Submitted by: Ashwin Goyal (UID - 115526297)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function averageHistogram(trainFolder, cropFolder, colorSpace)

    % Read all training image names
    trainFiles = dir([trainFolder '*.jpg']);
    % Read cropped buoy names
    cropFiles = dir([cropFolder '*.jpg']);
    
    % Find the index for each color
    greenIndex = [];
    redIndex = [];
    yellowIndex = [];
    for i = 1:length(cropFiles)
        if cropFiles(i).name(1) == 'G'
            greenIndex = [greenIndex; str2double(cropFiles(i).name(3:5))];
        elseif cropFiles(i).name(1) == 'R'
            redIndex = [redIndex; str2double(cropFiles(i).name(3:5))];
        elseif cropFiles(i).name(1) == 'Y'
            yellowIndex = [yellowIndex; str2double(cropFiles(i).name(3:5))];
        end
    end
    
    % Get the color space to be used
    possibleColorSpace = {'RGB';'HSV';'YIQ';'NTSC';'YCbCr';'Lab'};
    colorMatch = find(strcmpi(colorSpace,possibleColorSpace));
    if isempty(colorMatch) || length(colorMatch) > 1
        disp('The specified color space is incorrect.')
        return
    end
    % Find the color space values and form histogram for each buoy
    greenDist = [];  greenHist = zeros(256,3);  greenCount = 0;
    redDist = [];    redHist = zeros(256,3);    redCount = 0;
    yellowDist = []; yellowHist = zeros(256,3); yellowCount = 0;
    for i = 1:length(trainFiles)
        % Read the image in the format specified
        switch colorMatch
            case 1
                I = imread([trainFolder trainFiles(i).name]);
            case 2
                I = rgb2hsv(imread([trainFolder trainFiles(i).name]));
            case {3,4}
                I = rgb2ntsc(imread([trainFolder trainFiles(i).name]));
            case 5
                I = rgb2ycbcr(imread([trainFolder trainFiles(i).name]));
            case 6
                I = rgb2lab(imread([trainFolder trainFiles(i).name]));
        end
        plane1 = I(:,:,1);
        plane2 = I(:,:,2);
        plane3 = I(:,:,3);
        
        % Check if the image has green buoy
        if any(greenIndex == str2double(trainFiles(i).name(1:3)))
            % Read cropped green buoy image
            I_crop = imread([cropFolder 'G_' trainFiles(i).name]);
            
            % Get color distribution
            ROI = find(I_crop > 0);
            greenDist = [greenDist; plane1(ROI) plane2(ROI) plane3(ROI)];
            
            % Get histogram values
            [countR,~] = imhist(plane1(ROI));
            [countG,~] = imhist(plane2(ROI));
            [countB,~] = imhist(plane3(ROI));
            greenHist = greenHist + length(ROI)*[countR countG countB];
            greenCount = greenCount + length(ROI);
        end
        
        % Check if the image has red buoy
        if any(redIndex == str2double(trainFiles(i).name(1:3)))
            % Read cropped red buoy image
            I_crop = imread([cropFolder 'R_' trainFiles(i).name]);
            
            % Get color distribution
            ROI = find(I_crop > 0);
            redDist = [redDist; plane1(ROI) plane2(ROI) plane3(ROI)];
            
            % Get histogram values
            [countR,~] = imhist(plane1(ROI));
            [countG,~] = imhist(plane2(ROI));
            [countB,~] = imhist(plane3(ROI));
            redHist = redHist + length(ROI)*[countR countG countB];
            redCount = redCount + length(ROI);
        end
        
        % Check if the image has yellow buoy
        if any(yellowIndex == str2double(trainFiles(i).name(1:3)))
            % Read cropped yellow buoy image
            I_crop = imread([cropFolder 'Y_' trainFiles(i).name]);
            
            % Get color distribution
            ROI = find(I_crop > 0);
            yellowDist = [yellowDist; plane1(ROI) plane2(ROI) plane3(ROI)];
            
            % Get histogram values
            [countR,~] = imhist(plane1(ROI));
            [countG,~] = imhist(plane2(ROI));
            [countB,~] = imhist(plane3(ROI));
            yellowHist = yellowHist + length(ROI)*[countR countG countB];
            yellowCount = yellowCount + length(ROI);
        end
    end
    
    % Compute average histogram
    greenHist = (greenHist/greenCount);
    greenHist = greenHist/mean(sum(greenHist));
    redHist = redHist/redCount;
    redHist = redHist/mean(sum(redHist));
    yellowHist = yellowHist/yellowCount;
    yellowHist = yellowHist/mean(sum(yellowHist));
    
    if colorMatch == 1
        % Save image of the histogram for green buoy
        figure('units','normalized','outerposition',[0 0 1 1])
        bar(0:255,greenHist(:,3),'b');
        hold on
        bar(0:255,greenHist(:,1),'r');
        bar(0:255,greenHist(:,2),'g');
        hold off
        xlabel('Intensity')
        ylabel('Frequency')
        title('Normalized Average Histogram for Green Buoy');
        saveas(gcf,'.\ColorSeg\Output\Part0\G_hist.jpg');
        
        % Save image of the histogram for red buoy
        bar(0:255,redHist(:,3),'b');
        hold on
        bar(0:255,redHist(:,2),'g');
        bar(0:255,redHist(:,1),'r');
        hold off
        title('Normalized Average Histogram for Red Buoy');
        xlabel('Intensity')
        ylabel('Frequency')
        saveas(gcf,'.\ColorSeg\Output\Part0\R_hist.jpg');
        
        % Save image of the histogram for yellow buoy
        bar(0:255,yellowHist(:,3),'b');
        hold on
        bar(0:255,yellowHist(:,1),'r');
        bar(0:255,yellowHist(:,2),'g');
        hold off
        title('Normalized Average Histogram for Yellow Buoy');
        xlabel('Intensity')
        ylabel('Frequency')
        saveas(gcf,'.\ColorSeg\Output\Part0\Y_hist.jpg');
    end
    
    greenDist = double(greenDist); redDist = double(redDist); yellowDist = double(yellowDist);
    greenHist = double(greenHist); redHist = double(redHist); yellowHist = double(yellowHist);
    % Save the color distributions
    save(['.\ColorSeg\Output\Part0\colorDistributions_' colorSpace '.mat'],'greenDist','redDist','yellowDist')
    save(['.\ColorSeg\Output\Part0\colorHistograms_' colorSpace '.mat'],'greenHist','redHist','yellowHist')

end
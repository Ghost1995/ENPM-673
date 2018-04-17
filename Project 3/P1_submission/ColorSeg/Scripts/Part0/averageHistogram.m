%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code computes average histogram for individual colors or gives color
% distribution as output
% 
% Input:
%   colorSpace --> Color space to be used
% 
% Submitted by: Ashwin Goyal (UID - 115526297)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout = averageHistogram(colorSpace)

    % Define the folder of training set
    trainFolder = '..\..\Images\TrainingSet\Frames\';
    % Read all training image names
    trainFiles = dir([trainFolder '*.jpg']);
    
    % Define the folder of cropped buoys
    cropFolder = '..\..\Images\TrainingSet\CroppedBuoys\';
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
    
    % Find the color space values for each buoy
    greenBuoy = [];
    redBuoy = [];
    yellowBuoy = [];
    for i = 1:length(trainFiles)
        % Read the image in the format specified
        if strcmp(colorSpace,'RGB')
            I = imread([trainFolder trainFiles(i).name]);
        elseif strcmp(colorSpace,'HSV')
            I = rgb2hsv(imread([trainFolder trainFiles(i).name]));
        elseif strcmp(colorSpace,'YIQ')
            I = rgb2ntsc(imread([trainFolder trainFiles(i).name]));
        elseif strcmp(colorSpace,'YCbCr')
            I = rgb2ycbcr(imread([trainFolder trainFiles(i).name]));
        elseif strcmp(colorSpace,'Lab')
            I = rgb2lab(imread([trainFolder trainFiles(i).name]));
        else
            disp('Color Space incorrect.')
            return
        end
        plane1 = I(:,:,1);
        plane2 = I(:,:,2);
        plane3 = I(:,:,3);
        
        % Check if the image has green buoy, and add elements if it exists
        if any(greenIndex == str2double(trainFiles(i).name(1:3)))
            I_crop = imread([cropFolder 'G_' trainFiles(i).name]);
            ROI = find(I_crop > 0);
            greenBuoy = [greenBuoy; plane1(ROI) plane2(ROI) plane3(ROI)];
        end
        % Check if the image has red buoy, and add elements if it exists
        if any(redIndex == str2double(trainFiles(i).name(1:3)))
            I_crop = imread([cropFolder 'R_' trainFiles(i).name]);
            ROI = find(I_crop > 0);
            redBuoy = [redBuoy; plane1(ROI) plane2(ROI) plane3(ROI)];
        end
        % Check if the image has yellow buoy, and add elements if it exists
        if any(yellowIndex == str2double(trainFiles(i).name(1:3)))
            I_crop = imread([cropFolder 'Y_' trainFiles(i).name]);
            ROI = find(I_crop > 0);
            yellowBuoy = [yellowBuoy; plane1(ROI) plane2(ROI) plane3(ROI)];
        end
    end
    
    % Either save image or return the values of all colour buoys
    if nargout == 0
        figure
        scatter3(greenBuoy(:,1),greenBuoy(:,2),greenBuoy(:,3),'.')
        title('Color Distribubtion for Green Buoy');
        if strcmp(colorSpace,'RGB')
            xlabel('Red (R)');
            ylabel('Green (G)');
            zlabel('Blue (B)');
        elseif strcmp(colorSpace,'HSV')
            xlabel('Hue (H)');
            ylabel('Saturation (S)');
            zlabel('Value (V)');
        elseif strcmp(colorSpace,'YIQ')
            xlabel('Luminance (Y)');
            ylabel('Hue (I)');
            zlabel('Saturation (Q)');
        elseif strcmp(colorSpace,'YCbCr')
            xlabel('Luminance (Y)');
            ylabel('Chrominance (Cb)');
            zlabel('Chrominance (Cr)');
        elseif strcmp(colorSpace,'Lab')
            xlabel('Lightness (L)');
            ylabel('Color Component (a)');
            zlabel('Color Component (b)');
        end
        saveas(gcf,'../../Output/Part0/G_hist.jpg');
        
        figure
        scatter3(redBuoy(:,1),redBuoy(:,2),redBuoy(:,3),'.')
        title('Color Distribubtion for Red Buoy');
        if strcmp(colorSpace,'RGB')
            xlabel('Red (R)');
            ylabel('Green (G)');
            zlabel('Blue (B)');
        elseif strcmp(colorSpace,'HSV')
            xlabel('Hue (H)');
            ylabel('Saturation (S)');
            zlabel('Value (V)');
        elseif strcmp(colorSpace,'YIQ')
            xlabel('Luminance (Y)');
            ylabel('Hue (I)');
            zlabel('Saturation (Q)');
        elseif strcmp(colorSpace,'YCbCr')
            xlabel('Luminance (Y)');
            ylabel('Chrominance (Cb)');
            zlabel('Chrominance (Cr)');
        elseif strcmp(colorSpace,'Lab')
            xlabel('Lightness (L)');
            ylabel('Color Component (a)');
            zlabel('Color Component (b)');
        end
        saveas(gcf,'../../Output/Part0/R_hist.jpg');
        
        figure
        scatter3(yellowBuoy(:,1),yellowBuoy(:,2),yellowBuoy(:,3),'.')
        title('Color Distribubtion for Yellow Buoy');
        if strcmp(colorSpace,'RGB')
            xlabel('Red (R)');
            ylabel('Green (G)');
            zlabel('Blue (B)');
        elseif strcmp(colorSpace,'HSV')
            xlabel('Hue (H)');
            ylabel('Saturation (S)');
            zlabel('Value (V)');
        elseif strcmp(colorSpace,'YIQ')
            xlabel('Luminance (Y)');
            ylabel('Hue (I)');
            zlabel('Saturation (Q)');
        elseif strcmp(colorSpace,'YCbCr')
            xlabel('Luminance (Y)');
            ylabel('Chrominance (Cb)');
            zlabel('Chrominance (Cr)');
        elseif strcmp(colorSpace,'Lab')
            xlabel('Lightness (L)');
            ylabel('Color Component (a)');
            zlabel('Color Component (b)');
        end
        saveas(gcf,'../../Output/Part0/Y_hist.jpg');
    else
        varargout{1} = greenBuoy;
        if nargout > 1
            varargout{2} = redBuoy;
            if nargout > 2
                varargout{3} = yellowBuoy;
            end
        end
    end

    
%     % Initialize centers for the bar graphs
%     centerRGB = linspace(0,255,256);
%     centerRGB = centerRGB + (centerRGB(2) - centerRGB(1))/2;
%     centerHSV = linspace(0,1,256);
%     centerHSV = centerHSV + (centerHSV(2) - centerHSV(1))/2;
    
%     % Compute average histogram for green buoys
%     count_R = 0; count_G = 0; count_B = 0;
%     count_H = 0; count_S = 0; count_V = 0;
%     totalCount = 0;
%     for i = 1:length(greenIndex)
%         I = imread([cropFolder trainFiles(i).name]);
%         [counts,~] = imhist(I(:,:,1));
%         count_R = count_R + sum(counts)*counts;
%         [counts,~] = imhist(I(:,:,2));
%         count_G = count_G + sum(counts)*counts;
%         [counts,~] = imhist(I(:,:,3));
%         count_B = count_B + sum(counts)*counts;
%         I = rgb2hsv(I);
%         [counts,~] = imhist(I(:,:,1));
%         count_H = count_H + sum(counts)*counts;
%         [counts,~] = imhist(I(:,:,2));
%         count_S = count_S + sum(counts)*counts;
%         [counts,~] = imhist(I(:,:,3));
%         count_V = count_V + sum(counts)*counts;
%         totalCount = totalCount + numel(I)/3;
%     end
%     % If output required, give the data for average histogram
%     if nargout > 0
%         varargout{1} = cat(3,[count_R/totalCount, count_G/totalCount, count_B/totalCount],...
%                              [count_H/totalCount, count_S/totalCount, count_V/totalCount]);
%     else % Save average histogram
%         figure
%         bar(centerRGB,count_R/totalCount,'r');
%         hold on
%         bar(centerRGB,count_G/totalCount,'g');
%         bar(centerRGB,count_B/totalCount,'b');
%         hold off
%         xlabel('Intensity')
%         ylabel('Frequency')
%         title('Average Histogram for Green Buoys in RGB')
%         saveas(gcf,'../../Output/Part0/G_hist (RGB).jpg')
%         figure
%         bar(centerHSV,count_H/totalCount,'r');
%         hold on
%         bar(centerHSV,count_S/totalCount,'g');
%         bar(centerHSV,count_V/totalCount,'b');
%         hold off
%         xlabel('Intensity')
%         ylabel('Frequency')
%         title('Average Histogram for Green Buoys in HSV')
%         saveas(gcf,'../../Output/Part0/G_hist (HSV).jpg')
%     end
%     
%     % Compute average histogram for red buoys
%     count_R = 0; count_G = 0; count_B = 0;
%     count_H = 0; count_S = 0; count_V = 0;
%     totalCount = 0;
%     for i = 1:length(redIndex)
%         I = imread([cropFolder trainFiles(i).name]);
%         [counts,~] = imhist(I(:,:,1));
%         count_R = count_R + sum(counts)*counts;
%         [counts,~] = imhist(I(:,:,2));
%         count_G = count_G + sum(counts)*counts;
%         [counts,~] = imhist(I(:,:,3));
%         count_B = count_B + sum(counts)*counts;
%         I = rgb2hsv(I);
%         [counts,~] = imhist(I(:,:,1));
%         count_H = count_H + sum(counts)*counts;
%         [counts,~] = imhist(I(:,:,2));
%         count_S = count_S + sum(counts)*counts;
%         [counts,~] = imhist(I(:,:,3));
%         count_V = count_V + sum(counts)*counts;
%         totalCount = totalCount + numel(I)/3;
%     end
%     % If output required, give the data for average histogram
%     if nargout > 1
%         varargout{2} = cat(3,[count_R/totalCount, count_G/totalCount, count_B/totalCount],...
%                              [count_H/totalCount, count_S/totalCount, count_V/totalCount]);
%     else % Save average histogram
%         figure
%         bar(centerRGB,count_R/totalCount,'r');
%         hold on
%         bar(centerRGB,count_G/totalCount,'g');
%         bar(centerRGB,count_B/totalCount,'b');
%         hold off
%         xlabel('Intensity')
%         ylabel('Frequency')
%         title('Average Histogram for Red Buoys in RGB')
%         saveas(gcf,'../../Output/Part0/R_hist (RGB).jpg')
%         figure
%         bar(centerHSV,count_H/totalCount,'r');
%         hold on
%         bar(centerHSV,count_S/totalCount,'g');
%         bar(centerHSV,count_V/totalCount,'b');
%         hold off
%         xlabel('Intensity')
%         ylabel('Frequency')
%         title('Average Histogram for Red Buoys in HSV')
%         saveas(gcf,'../../Output/Part0/R_hist (HSV).jpg')
%     end
%     
%     % Compute average histogram for yellow buoys
%     count_R = 0; count_G = 0; count_B = 0;
%     count_H = 0; count_S = 0; count_V = 0;
%     totalCount = 0;
%     for i = 1:length(yellowIndex)
%         I = imread([cropFolder trainFiles(i).name]);
%         [counts,~] = imhist(I(:,:,1));
%         count_R = count_R + sum(counts)*counts;
%         [counts,~] = imhist(I(:,:,2));
%         count_G = count_G + sum(counts)*counts;
%         [counts,~] = imhist(I(:,:,3));
%         count_B = count_B + sum(counts)*counts;
%         I = rgb2hsv(I);
%         [counts,~] = imhist(I(:,:,1));
%         count_H = count_H + sum(counts)*counts;
%         [counts,~] = imhist(I(:,:,2));
%         count_S = count_S + sum(counts)*counts;
%         [counts,~] = imhist(I(:,:,3));
%         count_V = count_V + sum(counts)*counts;
%         totalCount = totalCount + numel(I)/3;
%     end
%     % If output required, give the data for average histogram
%     scatter3(count_R/totalCount,count_G/totalCount,count_B/totalCount,'.');
%     if nargout > 2
%         varargout{3} = cat(3,[count_R/totalCount, count_G/totalCount, count_B/totalCount],...
%                              [count_H/totalCount, count_S/totalCount, count_V/totalCount]);
%     else % Save average histogram
%         figure
%         bar(centerRGB,count_R/totalCount,'r');
%         hold on
%         bar(centerRGB,count_G/totalCount,'g');
%         bar(centerRGB,count_B/totalCount,'b');
%         hold off
%         xlabel('Intensity')
%         ylabel('Frequency')
%         title('Average Histogram for Yellow Buoys in RGB')
%         saveas(gcf,'../../Output/Part0/Y_hist (RGB).jpg')
%         figure
%         bar(centerHSV,count_H/totalCount,'r');
%         hold on
%         bar(centerHSV,count_S/totalCount,'g');
%         bar(centerHSV,count_V/totalCount,'b');
%         hold off
%         xlabel('Intensity')
%         ylabel('Frequency')
%         title('Average Histogram for Yellow Buoys in HSV')
%         saveas(gcf,'../../Output/Part0/Y_hist (HSV).jpg')
%     end

end
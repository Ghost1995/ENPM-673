%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code computes average histogram for individual colors
% 
% Input:
%   imageFolder --> Location of the cropped images of the buoy
% 
% Submitted by: Ashwin Goyal (UID - 115526297)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout = averageHistogram

    % Define the folder of cropped buoys
    imageFolder = '..\..\Images\TrainingSet\CroppedBuoys\';
    
    % Read image names
    imgFiles = dir([imageFolder '*.jpg']);
    
    % Find the index for each colour
    green = [];
    red = [];
    yellow = [];
    for i = 1:length(imgFiles)
        if imgFiles(i).name(1) == 'G'
            green = [green; i];
        elseif imgFiles(i).name(1) == 'R'
            red = [red; i];
        elseif imgFiles(i).name(1) == 'Y'
            yellow = [yellow; i];
        end
    end
    
    % Initialize centers for the bar graphs
    centerRGB = linspace(0,255,256);
    centerRGB = centerRGB + (centerRGB(2) - centerRGB(1))/2;
    centerHSV = linspace(0,1,256);
    centerHSV = centerHSV + (centerHSV(2) - centerHSV(1))/2;
    
    % Compute average histogram for green buoys
    count_R = 0; count_G = 0; count_B = 0;
    count_H = 0; count_S = 0; count_V = 0;
    totalCount = 0;
    for i = 1:length(green)
        I = imread([imageFolder imgFiles(i).name]);
        [counts,~] = imhist(I(:,:,1));
        count_R = count_R + sum(counts)*counts;
        [counts,~] = imhist(I(:,:,2));
        count_G = count_G + sum(counts)*counts;
        [counts,~] = imhist(I(:,:,3));
        count_B = count_B + sum(counts)*counts;
        I = rgb2hsv(I);
        [counts,~] = imhist(I(:,:,1));
        count_H = count_H + sum(counts)*counts;
        [counts,~] = imhist(I(:,:,2));
        count_S = count_S + sum(counts)*counts;
        [counts,~] = imhist(I(:,:,3));
        count_V = count_V + sum(counts)*counts;
        totalCount = totalCount + numel(I)/3;
    end
    % If output required, give the data for average histogram
    if nargout > 0
        varargout{1} = cat(3,[count_R/totalCount, count_G/totalCount, count_B/totalCount],...
                             [count_H/totalCount, count_S/totalCount, count_V/totalCount]);
    else % Save average histogram
        figure
        bar(centerRGB,count_R/totalCount,'r');
        hold on
        bar(centerRGB,count_G/totalCount,'g');
        bar(centerRGB,count_B/totalCount,'b');
        hold off
        xlabel('Intensity')
        ylabel('Frequency')
        title('Average Histogram for Green Buoys in RGB')
        saveas(gcf,'../../Output/Part0/G_hist (RGB).jpg')
        figure
        bar(centerHSV,count_H/totalCount,'r');
        hold on
        bar(centerHSV,count_S/totalCount,'g');
        bar(centerHSV,count_V/totalCount,'b');
        hold off
        xlabel('Intensity')
        ylabel('Frequency')
        title('Average Histogram for Green Buoys in HSV')
        saveas(gcf,'../../Output/Part0/G_hist (HSV).jpg')
    end
    
    % Compute average histogram for red buoys
    count_R = 0; count_G = 0; count_B = 0;
    count_H = 0; count_S = 0; count_V = 0;
    totalCount = 0;
    for i = 1:length(red)
        I = imread([imageFolder imgFiles(i).name]);
        [counts,~] = imhist(I(:,:,1));
        count_R = count_R + sum(counts)*counts;
        [counts,~] = imhist(I(:,:,2));
        count_G = count_G + sum(counts)*counts;
        [counts,~] = imhist(I(:,:,3));
        count_B = count_B + sum(counts)*counts;
        I = rgb2hsv(I);
        [counts,~] = imhist(I(:,:,1));
        count_H = count_H + sum(counts)*counts;
        [counts,~] = imhist(I(:,:,2));
        count_S = count_S + sum(counts)*counts;
        [counts,~] = imhist(I(:,:,3));
        count_V = count_V + sum(counts)*counts;
        totalCount = totalCount + numel(I)/3;
    end
    % If output required, give the data for average histogram
    if nargout > 1
        varargout{2} = cat(3,[count_R/totalCount, count_G/totalCount, count_B/totalCount],...
                             [count_H/totalCount, count_S/totalCount, count_V/totalCount]);
    else % Save average histogram
        figure
        bar(centerRGB,count_R/totalCount,'r');
        hold on
        bar(centerRGB,count_G/totalCount,'g');
        bar(centerRGB,count_B/totalCount,'b');
        hold off
        xlabel('Intensity')
        ylabel('Frequency')
        title('Average Histogram for Red Buoys in RGB')
        saveas(gcf,'../../Output/Part0/R_hist (RGB).jpg')
        figure
        bar(centerHSV,count_H/totalCount,'r');
        hold on
        bar(centerHSV,count_S/totalCount,'g');
        bar(centerHSV,count_V/totalCount,'b');
        hold off
        xlabel('Intensity')
        ylabel('Frequency')
        title('Average Histogram for Red Buoys in HSV')
        saveas(gcf,'../../Output/Part0/R_hist (HSV).jpg')
    end
    
    % Compute average histogram for yellow buoys
    count_R = 0; count_G = 0; count_B = 0;
    count_H = 0; count_S = 0; count_V = 0;
    totalCount = 0;
    for i = 1:length(yellow)
        I = imread([imageFolder imgFiles(i).name]);
        [counts,~] = imhist(I(:,:,1));
        count_R = count_R + sum(counts)*counts;
        [counts,~] = imhist(I(:,:,2));
        count_G = count_G + sum(counts)*counts;
        [counts,~] = imhist(I(:,:,3));
        count_B = count_B + sum(counts)*counts;
        I = rgb2hsv(I);
        [counts,~] = imhist(I(:,:,1));
        count_H = count_H + sum(counts)*counts;
        [counts,~] = imhist(I(:,:,2));
        count_S = count_S + sum(counts)*counts;
        [counts,~] = imhist(I(:,:,3));
        count_V = count_V + sum(counts)*counts;
        totalCount = totalCount + numel(I)/3;
    end
    % If output required, give the data for average histogram
    scatter3(count_R/totalCount,count_G/totalCount,count_B/totalCount,'.');
    if nargout > 2
        varargout{3} = cat(3,[count_R/totalCount, count_G/totalCount, count_B/totalCount],...
                             [count_H/totalCount, count_S/totalCount, count_V/totalCount]);
    else % Save average histogram
        figure
        bar(centerRGB,count_R/totalCount,'r');
        hold on
        bar(centerRGB,count_G/totalCount,'g');
        bar(centerRGB,count_B/totalCount,'b');
        hold off
        xlabel('Intensity')
        ylabel('Frequency')
        title('Average Histogram for Yellow Buoys in RGB')
        saveas(gcf,'../../Output/Part0/Y_hist (RGB).jpg')
        figure
        bar(centerHSV,count_H/totalCount,'r');
        hold on
        bar(centerHSV,count_S/totalCount,'g');
        bar(centerHSV,count_V/totalCount,'b');
        hold off
        xlabel('Intensity')
        ylabel('Frequency')
        title('Average Histogram for Yellow Buoys in HSV')
        saveas(gcf,'../../Output/Part0/Y_hist (HSV).jpg')
    end

end
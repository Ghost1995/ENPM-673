%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is Homework 0. The objective of this assignment is to identify the
% different coloured pins on a white background. If possible, also identify
% the white and transparent pins.
% 
% Submitted by: Ashwin Goyal
% University ID: 115526297
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read Image
I = imread('../input/TestImgResized.jpg');
Ioriginal = I;

% Sharpen Image
I = imsharpen(I);

% RGB to Grayscale
Igray = rgb2gray(I);

% Edge Detection
Idetect = edge(Igray);

% Morphological Operators
Idetect = bwmorph(Idetect,'fill');
Idetect = bwmorph(Idetect,'clean');
Idetect = bwmorph(Idetect,'close');

% Identify all Pins
count = 0;
pin = 0;
notPin = [];
whiteSpace = [];
indexUsed = [];
for i=1:size(Idetect,1)
    for j=1:size(Idetect,2)
        if ~isempty(indexUsed)
            if any((i>=indexUsed(:,1))&(j>=indexUsed(:,2))&(i<=indexUsed(:,3))&(j<=indexUsed(:,4)))
                continue;
            end
        end
        if Idetect(i,j)
            for k=i:i+40
                if (k > size(Idetect,1))
                    continue;
                end
                for l=j-40:j+40
                    if (l < 1)||(l > size(Idetect,2))
                        continue;
                    end
                    if Idetect(k,l)
                        count = count+1;
                        whiteSpace(count,:) = [k l];
                    end
                end
            end
            pin = pin+1;
            indexUsed(pin,:) = [min(whiteSpace) max(whiteSpace)];
            if (count < 100)
                notPin = [notPin;pin];
            end
            count = 0;
            whiteSpace = [];
        end
    end
end

figure
imshow(Ioriginal)
hold on

% Count Colours
redCount = 0;
greenCount = 0;
blueCount = 0;
yellowCount = 0;
whiteCount = 0;
for i=1:size(indexUsed,1)
    if ~any(i == notPin)
        Icrop = imcrop(I,[indexUsed(i,2) indexUsed(i,1) indexUsed(i,4)-indexUsed(i,2) indexUsed(i,3)-indexUsed(i,1)]);
        Ibin = imbinarize(Icrop);
        if sum(sum(Ibin(:,:,1))) > 0.75*numel(Ibin(:,:,1))
            isRed = true;
        else
            isRed = false;
        end
        if sum(sum(Ibin(:,:,2))) > 0.61*numel(Ibin(:,:,2))
            isGreen = true;
        else
            isGreen = false;
        end
        if sum(sum(Ibin(:,:,3))) > 0.75*numel(Ibin(:,:,3))
            isBlue = true;
        else
            isBlue = false;
        end
        x = [indexUsed(i,2); indexUsed(i,4); indexUsed(i,4); indexUsed(i,2)];
        y = [indexUsed(i,1); indexUsed(i,1); indexUsed(i,3); indexUsed(i,3)];
        if (isRed)&&(~isGreen)&&(~isBlue)
            redCount = redCount + 1;
            patch(x,y,'red','FaceAlpha',0.5);
        elseif (~isRed)&&(isGreen)&&(~isBlue)
            greenCount = greenCount + 1;
            patch(x,y,'green','FaceAlpha',0.5);
        elseif (~isRed)&&(~isGreen)&&(isBlue)
            blueCount = blueCount + 1;
            patch(x,y,'blue','FaceAlpha',0.5);
        elseif (isRed)&&(isGreen)&&(~isBlue)
            yellowCount = yellowCount + 1;
            patch(x,y,'yellow','FaceAlpha',0.5);
        else
            whiteCount = whiteCount + 1;
            patch(x,y,'white','FaceAlpha',0.5);
        end
    end
end
% Save the figure
saveas(gcf,'../output/Identified_Pins.jpg')
disp(['There are ' num2str(redCount) ' red pins.']);
disp(['There are ' num2str(greenCount) ' green pins.']);
disp(['There are ' num2str(blueCount) ' blue pins.']);
disp(['There are ' num2str(yellowCount) ' yellow pins.']);
disp(['There are a total of ' num2str(redCount + greenCount + blueCount + yellowCount) ' coloured pins.']);
disp(['There are ' num2str(whiteCount) ' white or transparent pins.']);

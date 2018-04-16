%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code computes average histogram for individual colors
% 
% Input:
%   imageFolder --> Location of the cropped images of the buoy
% 
% Submitted by: Ashwin Goyal (UID - 115526297)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function varargout = segment1D(img)

    % Define the folder of cropped buoys
    imageFolder = '..\..\Images\TrainingSet\CroppedBuoys\';
    
    % Read image names
    imgFiles = dir([imageFolder '*.jpg']);
    
    % Compute average color histogram for all images
    [greenHist,redHist,yellowHist] = averageHistogram;
    
    % Create 1-D gaussians
    value = (0.5:255.5)';
    [muG,sigmaG] = normfit(value,0.05,zeros(size(value)),greenHist);
    [muR,sigmaR] = normfit(value,0.05,zeros(size(value)),redHist);
    [muY,sigmaY] = normfit(value,0.05,zeros(size(value)),yellowHist);
    totalHist = (sum(greenHist)*greenHist + sum(redHist)*redHist + sum(yellowHist)*yellowHist)/sum(greenHist + redHist + yellowHist);
    [mu,sigma] = normfit(value,0.05,zeros(size(value)),totalHist);
    
    % Save gaussians being used
    normG = normpdf(value,muG,sigmaG);
    normR = normpdf(value,muR,sigmaR);
    normY = normpdf(value,muY,sigmaY);
    norm = normpdf(value,mu,sigma);
    plot(value,normG);
    saveas(gcf,'../../Output/Part0/G_gauss1D.jpg')
    plot(value,normR);
    saveas(gcf,'../../Output/Part0/R_gauss1D.jpg')
    plot(value,normY);
    saveas(gcf,'../../Output/Part0/Y_gauss1D.jpg')
    plot(value,norm);
    saveas(gcf,'../../Output/Part0/gauss1D.jpg')
    
    
% end
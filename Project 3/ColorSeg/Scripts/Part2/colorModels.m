%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code computes N-D gaussian models for each buoy
% 
% Input:
%   colorSpace --> Color space to be used
%     plotPath --> Path at which the images need to be saved
%         Nmin --> Minimum number of gaussians to be tested
%         Nmax --> Maximum number of gaussians to be tested
%            D --> The dimension of gaussian to be used
% 
% Submitted by: Ashwin Goyal (UID - 115526297)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function colorModels(colorSpace, plotPath, Nmin, Nmax, D)

    % Get the color space to be used
    possibleColorSpace = {'RGB';'HSV';'YCbCr'};
    colorMatch = find(strcmpi(colorSpace,possibleColorSpace));
    if isempty(colorMatch) || length(colorMatch) > 1
        disp('The specified color space is not supported.')
        return
    end
    
    greenHist = []; redHist = []; yellowHist = [];
    greenDist = []; redDist = []; yellowDist = [];
    try
        load(['.\ColorSeg\Output\Part0\colorHistograms_' colorSpace '.mat'],'greenHist','redHist','yellowHist')
        load(['.\ColorSeg\Output\Part0\colorDistributions_' colorSpace '.mat'],'greenDist','redDist','yellowDist')
    catch
        disp('Color Data not found. First compute them using averageHistogram.m')
        return
    end
    if colorMatch == 2
        greenDist = 255*greenDist; redDist = 255*redDist; yellowDist = 255*yellowDist;
    end
    
    if D == 1
        figure('units','normalized','outerposition',[0 0 1 1])
        for N = Nmin:Nmax
            bar(0:255,greenHist(:,1),'r');
            hold on
            greenDist_1D = sort(greenDist(:,1));
            for i = 1:10
                [greenObj_1D,isConverged] = EM(greenDist_1D,N);
                if isConverged
                    plot(greenDist_1D,gauss(greenObj_1D,greenDist_1D),'k')
                    xlabel('Intensity')
                    ylabel('Frequency')
                    title(['Red Color Band of Green Buoy represented by ' num2str(N) ' Gaussians'])
                    saveas(gcf,[plotPath 'greenGauss_red_' num2str(N) '.jpg']);
                    break;
                end
            end
            hold off

            bar(0:255,greenHist(:,2),'g');
            hold on
            greenDist_1D = sort(greenDist(:,2));
            for i = 1:10
                [greenObj_1D,isConverged] = EM(greenDist_1D,N);
                if isConverged
                    plot(greenDist_1D,gauss(greenObj_1D,greenDist_1D),'k')
                    xlabel('Intensity')
                    ylabel('Frequency')
                    title(['Green Color Band of Green Buoy represented by ' num2str(N) ' Gaussians'])
                    saveas(gcf,[plotPath 'greenGauss_green_' num2str(N) '.jpg']);
                    break;
                end
            end
            hold off

            bar(0:255,greenHist(:,3),'b');
            hold on
            greenDist_1D = sort(greenDist(:,3));
            for i = 1:10
                [greenObj_1D,isConverged] = EM(greenDist_1D,N);
                if isConverged
                    plot(greenDist_1D,gauss(greenObj_1D,greenDist_1D),'k')
                    xlabel('Intensity')
                    ylabel('Frequency')
                    title(['Blue Color Band of Green Buoy represented by ' num2str(N) ' Gaussians'])
                    saveas(gcf,[plotPath 'greenGauss_blue_' num2str(N) '.jpg']);
                    break;
                end
            end
            hold off

            bar(0:255,redHist(:,1),'r');
            hold on
            redDist_1D = sort(redDist(:,1));
            for i = 1:10
                [redObj_1D,isConverged] = EM(redDist_1D,N);
                if isConverged
                    plot(redDist_1D,gauss(redObj_1D,redDist_1D),'k')
                    xlabel('Intensity')
                    ylabel('Frequency')
                    title(['Red Color Band of Red Buoy represented by ' num2str(N) ' Gaussians'])
                    saveas(gcf,[plotPath 'redGauss_red_' num2str(N) '.jpg']);
                    break;
                end
            end
            hold off

            bar(0:255,redHist(:,2),'g');
            hold on
            redDist_1D = sort(redDist(:,2));
            for i = 1:10
                [redObj_1D,isConverged] = EM(redDist_1D,N);
                if isConverged
                    plot(redDist_1D,gauss(redObj_1D,redDist_1D),'k')
                    xlabel('Intensity')
                    ylabel('Frequency')
                    title(['Green Color Band of Red Buoy represented by ' num2str(N) ' Gaussians'])
                    saveas(gcf,[plotPath 'redGauss_green_' num2str(N) '.jpg']);
                    break;
                end
            end
            hold off

            bar(0:255,redHist(:,3),'b');
            hold on
            redDist_1D = sort(redDist(:,3));
            for i = 1:10
                [redObj_1D,isConverged] = EM(redDist_1D,N);
                if isConverged
                    plot(redDist_1D,gauss(redObj_1D,redDist_1D),'k')
                    xlabel('Intensity')
                    ylabel('Frequency')
                    title(['Blue Color Band of Red Buoy represented by ' num2str(N) ' Gaussians'])
                    saveas(gcf,[plotPath 'redGauss_blue_' num2str(N) '.jpg']);
                    break;
                end
            end
            hold off

            bar(0:255,yellowHist(:,1),'r');
            hold on
            yellowDist_1D = sort(yellowDist(:,1));
            for i = 1:10
                [yellowObj_1D,isConverged] = EM(yellowDist_1D,N);
                if isConverged
                    plot(yellowDist_1D,gauss(yellowObj_1D,yellowDist_1D),'k')
                    xlabel('Intensity')
                    ylabel('Frequency')
                    title(['Red Color Band of Yellow Buoy represented by ' num2str(N) ' Gaussians'])
                    saveas(gcf,[plotPath 'yellowGauss_red_' num2str(N) '.jpg']);
                    break;
                end
            end
            hold off

            bar(0:255,yellowHist(:,2),'g');
            hold on
            yellowDist_1D = sort(yellowDist(:,2));
            for i = 1:10
                [yellowObj_1D,isConverged] = EM(yellowDist_1D,N);
                if isConverged
                    plot(yellowDist_1D,gauss(yellowObj_1D,yellowDist_1D),'k')
                    xlabel('Intensity')
                    ylabel('Frequency')
                    title(['Green Color Band of Yellow Buoy represented by ' num2str(N) ' Gaussians'])
                    saveas(gcf,[plotPath 'yellowGauss_green_' num2str(N) '.jpg']);
                    break;
                end
            end
            hold off

            bar(0:255,yellowHist(:,3),'b');
            hold on
            yellowDist_1D = sort(yellowDist(:,3));
            for i = 1:10
                [yellowObj_1D,isConverged] = EM(yellowDist_1D,N);
                if isConverged
                    plot(yellowDist_1D,gauss(yellowObj_1D,yellowDist_1D),'k')
                    xlabel('Intensity')
                    ylabel('Frequency')
                    title(['Blue Color Band of Yellow Buoy represented by ' num2str(N) ' Gaussians'])
                    saveas(gcf,[plotPath 'yellowGauss_blue_' num2str(N) '.jpg']);
                    break;
                end
            end
            hold off
        end
    elseif D == 2
        rgGreenHist = zeros(256); gbGreenHist = zeros(256); rbGreenHist = zeros(256);
        rgRedHist = zeros(256); gbRedHist = zeros(256); rbRedHist = zeros(256);
        rgYellowHist = zeros(256); gbYellowHist = zeros(256); rbYellowHist = zeros(256);
        for i = 0:255
            if any(greenDist(:,1) >= i-0.5 & greenDist(:,1) < i+0.5)
                rows = greenDist(:,1) >= i-0.5 & greenDist(:,1) < i+0.5;
                for j = 0:255
                    rgGreenHist(i+1,j+1) = sum(greenDist(rows,2) >= j-0.5 & greenDist(rows,2) < j+0.5);
                    rbGreenHist(i+1,j+1) = sum(greenDist(rows,3) >= j-0.5 & greenDist(rows,3) < j+0.5);
                end
            end
            if any(greenDist(:,2) >= i-0.5 & greenDist(:,2) < i+0.5)
                rows = greenDist(:,2) >= i-0.5 & greenDist(:,2) < i+0.5;
                for j = 0:255
                    gbGreenHist(i+1,j+1) = sum(greenDist(rows,3) >= j-0.5 & greenDist(rows,3) < j+0.5);
                end
            end
            if any(redDist(:,1) >= i-0.5 & redDist(:,1) < i+0.5)
                rows = redDist(:,1) >= i-0.5 & redDist(:,1) < i+0.5;
                for j = 0:255
                    rgRedHist(i+1,j+1) = sum(redDist(rows,2) >= j-0.5 & redDist(rows,2) < j+0.5);
                    rbRedHist(i+1,j+1) = sum(redDist(rows,3) >= j-0.5 & redDist(rows,3) < j+0.5);
                end
            end
            if any(redDist(:,2) >= i-0.5 & redDist(:,2) < i+0.5)
                rows = redDist(:,2) >= i-0.5 & redDist(:,2) < i+0.5;
                for j = 0:255
                    gbRedHist(i+1,j+1) = sum(redDist(rows,3) >= j-0.5 & redDist(rows,3) < j+0.5);
                end
            end
            if any(yellowDist(:,1) >= i-0.5 & yellowDist(:,1) < i+0.5)
                rows = yellowDist(:,1) >= i-0.5 & yellowDist(:,1) < i+0.5;
                for j = 0:255
                    rgYellowHist(i+1,j+1) = sum(yellowDist(rows,2) >= j-0.5 & yellowDist(rows,2) < j+0.5);
                    rbYellowHist(i+1,j+1) = sum(yellowDist(rows,3) >= j-0.5 & yellowDist(rows,3) < j+0.5);
                end
            end
            if any(yellowDist(:,2) >= i-0.5 & yellowDist(:,2) < i+0.5)
                rows = yellowDist(:,2) >= i-0.5 & yellowDist(:,2) < i+0.5;
                for j = 0:255
                    gbYellowHist(i+1,j+1) = sum(yellowDist(rows,3) >= j-0.5 & yellowDist(rows,3) < j+0.5);
                end
            end
        end
        rgGreenHist = rgGreenHist/length(greenDist);
        gbGreenHist = gbGreenHist/length(greenDist);
        rbGreenHist = rbGreenHist/length(greenDist);
        rgRedHist = rgRedHist/length(redDist);
        gbRedHist = gbRedHist/length(redDist);
        rbRedHist = rbRedHist/length(redDist);
        rgYellowHist = rgYellowHist/length(yellowDist);
        gbYellowHist = gbYellowHist/length(yellowDist);
        rbYellowHist = rbYellowHist/length(yellowDist);
        
        figure('units','normalized','outerposition',[0 0 1 1])
        for N = Nmin:Nmax
            bar3(0:255,rgGreenHist','y');
            hold on
            greenDist_2D = greenDist(:,1:2);
            for i = 1:10
                [greenObj_2D,isConverged] = EM(greenDist_2D,N);
                if isConverged
                    plot3(greenDist_2D(:,1),greenDist_2D(:,2),gauss(greenObj_2D,greenDist_2D),'k*')
                    xlabel('Red Intensity')
                    ylabel('Green Intensity')
                    zlabel('Frequency')
                    title(['Red-Green Color Band of Green Buoy represented by ' num2str(N) ' Gaussians'])
                    view(-45,10)
                    saveas(gcf,[plotPath 'greenGauss_red-green_' num2str(N) '.jpg']);
                    break;
                end
            end
            hold off

            bar3(0:255,gbGreenHist','c');
            hold on
            greenDist_2D = greenDist(:,2:3);
            for i = 1:10
                [greenObj_2D,isConverged] = EM(greenDist_2D,N);
                if isConverged
                    plot3(greenDist_2D(:,1),greenDist_2D(:,2),gauss(greenObj_2D,greenDist_2D),'k*')
                    xlabel('Green Intensity')
                    ylabel('Blue Intensity')
                    zlabel('Frequency')
                    title(['Green-Blue Color Band of Green Buoy represented by ' num2str(N) ' Gaussians'])
                    view(-45,10)
                    saveas(gcf,[plotPath 'greenGauss_green-blue_' num2str(N) '.jpg']);
                    break;
                end
            end
            hold off
            
            bar3(0:255,rbGreenHist','m');
            hold on
            greenDist_2D = greenDist(:,[1,3]);
            for i = 1:10
                [greenObj_2D,isConverged] = EM(greenDist_2D,N);
                if isConverged
                    plot3(greenDist_2D(:,1),greenDist_2D(:,2),gauss(greenObj_2D,greenDist_2D),'k*')
                    xlabel('Red Intensity')
                    ylabel('Blue Intensity')
                    zlabel('Frequency')
                    title(['Red-Blue Color Band of Green Buoy represented by ' num2str(N) ' Gaussians'])
                    view(-45,10)
                    saveas(gcf,[plotPath 'greenGauss_red-blue_' num2str(N) '.jpg']);
                    break;
                end
            end
            hold off

            bar3(0:255,rgRedHist','y');
            hold on
            redDist_2D = redDist(:,1:2);
            for i = 1:10
                [redObj_2D,isConverged] = EM(redDist_2D,N);
                if isConverged
                    plot3(redDist_2D(:,1),redDist_2D(:,2),gauss(redObj_2D,redDist_2D),'k*')
                    xlabel('Red Intensity')
                    ylabel('Green Intensity')
                    zlabel('Frequency')
                    title(['Red-Green Color Band of Red Buoy represented by ' num2str(N) ' Gaussians'])
                    view(-45,10)
                    saveas(gcf,[plotPath 'redGauss_red-green_' num2str(N) '.jpg']);
                    break;
                end
            end
            hold off

            bar3(0:255,gbRedHist','c');
            hold on
            redDist_2D = redDist(:,2:3);
            for i = 1:10
                [redObj_2D,isConverged] = EM(redDist_2D,N);
                if isConverged
                    plot3(redDist_2D(:,1),redDist_2D(:,2),gauss(redObj_2D,redDist_2D),'k*')
                    xlabel('Green Intensity')
                    ylabel('Blue Intensity')
                    zlabel('Frequency')
                    title(['Green-Blue Color Band of Red Buoy represented by ' num2str(N) ' Gaussians'])
                    view(-45,10)
                    saveas(gcf,[plotPath 'redGauss_green-blue_' num2str(N) '.jpg']);
                    break;
                end
            end
            hold off
            
            bar3(0:255,rbRedHist','m');
            hold on
            redDist_2D = redDist(:,[1,3]);
            for i = 1:10
                [redObj_2D,isConverged] = EM(redDist_2D,N);
                if isConverged
                    plot3(redDist_2D(:,1),redDist_2D(:,2),gauss(redObj_2D,redDist_2D),'k*')
                    xlabel('Red Intensity')
                    ylabel('Blue Intensity')
                    zlabel('Frequency')
                    title(['Red-Blue Color Band of Red Buoy represented by ' num2str(N) ' Gaussians'])
                    view(-45,10)
                    saveas(gcf,[plotPath 'redGauss_red-blue_' num2str(N) '.jpg']);
                    break;
                end
            end
            hold off

            bar3(0:255,rgYellowHist','y');
            hold on
            yellowDist_2D = yellowDist(:,1:2);
            for i = 1:10
                [yellowObj_2D,isConverged] = EM(yellowDist_2D,N);
                if isConverged
                    plot3(yellowDist_2D(:,1),yellowDist_2D(:,2),gauss(yellowObj_2D,yellowDist_2D),'k*')
                    xlabel('Red Intensity')
                    ylabel('Green Intensity')
                    zlabel('Frequency')
                    title(['Red-Green Color Band of Yellow Buoy represented by ' num2str(N) ' Gaussians'])
                    view(-45,10)
                    saveas(gcf,[plotPath 'yellowGauss_red-green_' num2str(N) '.jpg']);
                    break;
                end
            end
            hold off

            bar3(0:255,gbYellowHist','c');
            hold on
            yellowDist_2D = yellowDist(:,2:3);
            for i = 1:10
                [yellowObj_2D,isConverged] = EM(yellowDist_2D,N);
                if isConverged
                    plot3(yellowDist_2D(:,1),yellowDist_2D(:,2),gauss(yellowObj_2D,yellowDist_2D),'k*')
                    xlabel('Green Intensity')
                    ylabel('Blue Intensity')
                    zlabel('Frequency')
                    title(['Green-Blue Color Band of Yellow Buoy represented by ' num2str(N) ' Gaussians'])
                    view(-45,10)
                    saveas(gcf,[plotPath 'yellowGauss_green-blue_' num2str(N) '.jpg']);
                    break;
                end
            end
            hold off
            
            bar3(0:255,rbYellowHist','m');
            hold on
            yellowDist_2D = yellowDist(:,[1,3]);
            for i = 1:10
                [yellowObj_2D,isConverged] = EM(yellowDist_2D,N);
                if isConverged
                    plot3(yellowDist_2D(:,1),yellowDist_2D(:,2),gauss(yellowObj_2D,yellowDist_2D),'k*')
                    xlabel('Red Intensity')
                    ylabel('Blue Intensity')
                    zlabel('Frequency')
                    title(['Red-Blue Color Band of Yellow Buoy represented by ' num2str(N) ' Gaussians'])
                    view(-45,10)
                    saveas(gcf,[plotPath 'yellowGauss_red-blue_' num2str(N) '.jpg']);
                    break;
                end
            end
            hold off
        end
    end

end

function N = gauss(gmObj, X)
% This function computes N(x|mu,sigma) for N-D Gaussian

    mean = gmObj.mu;
    sigma = gmObj.Sigma;
    mixtureCoeff = gmObj.ComponentProportion;
    N = 0;
    for i = 1:length(mixtureCoeff)
        N = N + mixtureCoeff(i)*(1/(2*pi)^(size(X,2)/2))*(1/sqrt(det(sigma(:,:,i))))*exp(sum(-0.5*((X - mean(i,:))/sigma(:,:,i)).*(X - mean(i,:)),2));
    end

end
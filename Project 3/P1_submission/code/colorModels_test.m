%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code computes N-D gaussian models for each buoy
% 
% Input:
%   colorSpace --> Color space to be used
%     plotPath --> Path at which the images need to be saved
%         Nmax --> Maximum number of gaussians to be tested
%            D --> The dimension of gaussian to be used
% 
% Submitted by: Ashwin Goyal (UID - 115526297)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function colorModels_test(colorSpace, plotPath, Nmax, D)

    greenHist = []; redHist = []; yellowHist = [];
    greenDist = []; redDist = []; yellowDist = [];
    try
        load(['..\output\colorHistograms_' colorSpace '.mat'],'greenHist','redHist','yellowHist')
        load(['..\output\colorDistributions_' colorSpace '.mat'],'greenDist','redDist','yellowDist')
    catch
        disp('Color Data not found. First compute them using averageHistogram.m')
        return
    end

    if D == 1
        for N = 1:Nmax
            figure('units','normalized','outerposition',[0 0 1 1])
            bar(0:255,greenHist(:,1),'r');
            hold on
            greenDist_1D = double(sort(greenDist(:,1)));
            [greenObj_1D,isConverged] = EM(greenDist_1D,N);
            if isConverged
                greenY_1D = pdf(greenObj_1D,greenDist_1D);
                plot(greenDist_1D,greenY_1D,'k')
                xlabel('Intensity')
                ylabel('Frequency')
                title(['Red Color Band of Green Buoy represented by ' num2str(N) ' Gaussians'])
                saveas(gcf,[plotPath 'greenGauss_red_' num2str(N) '.jpg']);
            end
            hold off

            bar(0:255,greenHist(:,2),'g');
            hold on
            greenDist_1D = double(sort(greenDist(:,2)));
            [greenObj_1D,isConverged] = EM(greenDist_1D,N);
            if isConverged
                greenY_1D = pdf(greenObj_1D,greenDist_1D);
                plot(greenDist_1D,greenY_1D,'k')
                xlabel('Intensity')
                ylabel('Frequency')
                title(['Green Color Band of Green Buoy represented by ' num2str(N) ' Gaussians'])
                saveas(gcf,[plotPath 'greenGauss_green_' num2str(N) '.jpg']);
            end
            hold off

            bar(0:255,greenHist(:,3),'b');
            hold on
            greenDist_1D = double(sort(greenDist(:,3)));
            [greenObj_1D,isConverged] = EM(greenDist_1D,N);
            if isConverged
                greenY_1D = pdf(greenObj_1D,greenDist_1D);
                plot(greenDist_1D,greenY_1D,'k')
                xlabel('Intensity')
                ylabel('Frequency')
                title(['Blue Color Band of Green Buoy represented by ' num2str(N) ' Gaussians'])
                saveas(gcf,[plotPath 'greenGauss_blue_' num2str(N) '.jpg']);
            end
            hold off

            figure('units','normalized','outerposition',[0 0 1 1])
            bar(0:255,redHist(:,1),'r');
            hold on
            redDist_1D = double(sort(redDist(:,1)));
            [redObj_1D,isConverged] = EM(redDist_1D,N);
            if isConverged
                redY_1D = pdf(redObj_1D,redDist_1D);
                plot(redDist_1D,redY_1D,'k')
                xlabel('Intensity')
                ylabel('Frequency')
                title(['Red Color Band of Red Buoy represented by ' num2str(N) ' Gaussians'])
                saveas(gcf,[plotPath 'redGauss_red_' num2str(N) '.jpg']);
            end
            hold off

            bar(0:255,redHist(:,2),'g');
            hold on
            redDist_1D = double(sort(redDist(:,2)));
            [redObj_1D,isConverged] = EM(redDist_1D,N);
            if isConverged
                redY_1D = pdf(redObj_1D,redDist_1D);
                plot(redDist_1D,redY_1D,'k')
                xlabel('Intensity')
                ylabel('Frequency')
                title(['Green Color Band of Red Buoy represented by ' num2str(N) ' Gaussians'])
                saveas(gcf,[plotPath 'redGauss_green_' num2str(N) '.jpg']);
            end
            hold off

            bar(0:255,redHist(:,3),'b');
            hold on
            redDist_1D = double(sort(redDist(:,3)));
            [redObj_1D,isConverged] = EM(redDist_1D,N);
            if isConverged
                redY_1D = pdf(redObj_1D,redDist_1D);
                plot(redDist_1D,redY_1D,'k')
                xlabel('Intensity')
                ylabel('Frequency')
                title(['Blue Color Band of Red Buoy represented by ' num2str(N) ' Gaussians'])
                saveas(gcf,[plotPath 'redGauss_blue_' num2str(N) '.jpg']);
            end
            hold off

            figure('units','normalized','outerposition',[0 0 1 1])
            bar(0:255,yellowHist(:,1),'r');
            hold on
            yellowDist_1D = double(sort(yellowDist(:,1)));
            [yellowObj_1D,isConverged] = EM(yellowDist_1D,N);
            if isConverged
                yellowY_1D = pdf(yellowObj_1D,yellowDist_1D);
                plot(yellowDist_1D,yellowY_1D,'k')
                xlabel('Intensity')
                ylabel('Frequency')
                title(['Red Color Band of Yellow Buoy represented by ' num2str(N) ' Gaussians'])
                saveas(gcf,[plotPath 'yellowGauss_red_' num2str(N) '.jpg']);
            end
            hold off

            bar(0:255,yellowHist(:,2),'g');
            hold on
            yellowDist_1D = double(sort(yellowDist(:,2)));
            [yellowObj_1D,isConverged] = EM(yellowDist_1D,N);
            if isConverged
                yellowY_1D = pdf(yellowObj_1D,yellowDist_1D);
                plot(yellowDist_1D,yellowY_1D,'k')
                xlabel('Intensity')
                ylabel('Frequency')
                title(['Green Color Band of Yellow Buoy represented by ' num2str(N) ' Gaussians'])
                saveas(gcf,[plotPath 'yellowGauss_green_' num2str(N) '.jpg']);
            end
            hold off

            bar(0:255,yellowHist(:,3),'b');
            hold on
            yellowDist_1D = double(sort(yellowDist(:,3)));
            [yellowObj_1D,isConverged] = EM(yellowDist_1D,N);
            if isConverged
                yellowY_1D = pdf(yellowObj_1D,yellowDist_1D);
                plot(yellowDist_1D,yellowY_1D,'k')
                xlabel('Intensity')
                ylabel('Frequency')
                title(['Blue Color Band of Yellow Buoy represented by ' num2str(N) ' Gaussians'])
                saveas(gcf,[plotPath 'yellowGauss_blue_' num2str(N) '.jpg']);
            end
            hold off
        end
    elseif D == 2
        rgGreenHist = zeros(256); gbGreenHist = zeros(256); rbGreenHist = zeros(256);
        rgRedHist = zeros(256); gbRedHist = zeros(256); rbRedHist = zeros(256);
        rgYellowHist = zeros(256); gbYellowHist = zeros(256); rbYellowHist = zeros(256);
        for i = 0:255
            if any(greenDist(:,1) == i)
                for j = 0:255
                    rgGreenHist(i+1,j+1) = sum(greenDist(greenDist(:,1) == i,2) == j);
                    rbGreenHist(i+1,j+1) = sum(greenDist(greenDist(:,1) == i,3) == j);
                end
            end
            if any(greenDist(:,2) == i)
                for j = 0:255
                    gbGreenHist(i+1,j+1) = sum(greenDist(greenDist(:,2) == i,3) == j);
                end
            end
            if any(redDist(:,1) == i)
                for j = 0:255
                    rgRedHist(i+1,j+1) = sum(redDist(redDist(:,1) == i,2) == j);
                    rbRedHist(i+1,j+1) = sum(redDist(redDist(:,1) == i,3) == j);
                end
            end
            if any(redDist(:,2) == i)
                for j = 0:255
                    gbRedHist(i+1,j+1) = sum(redDist(redDist(:,2) == i,3) == j);
                end
            end
            if any(yellowDist(:,1) == i)
                for j = 0:255
                    rgYellowHist(i+1,j+1) = sum(yellowDist(yellowDist(:,1) == i,2) == j);
                    rbYellowHist(i+1,j+1) = sum(yellowDist(yellowDist(:,1) == i,3) == j);
                end
            end
            if any(yellowDist(:,2) == i)
                for j = 0:255
                    gbYellowHist(i+1,j+1) = sum(yellowDist(yellowDist(:,2) == i,3) == j);
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
        for N = 1:Nmax
            figure('units','normalized','outerposition',[0 0 1 1])
            bar3(0:255,rgGreenHist','y');
            hold on
            greenDist_2D = double(greenDist(:,1:2));
            [greenObj_2D,isConverged] = EM(greenDist_2D,N);
            if isConverged
                greenY_2D = pdf(greenObj_2D,greenDist_2D);
                plot3(greenDist_2D(:,1),greenDist_2D(:,2),greenY_2D,'k*')
                xlabel('Red Intensity')
                ylabel('Green Intensity')
                zlabel('Frequency')
                title(['Red-Green Color Band of Green Buoy represented by ' num2str(N) ' Gaussians'])
                view(-45,10)
                saveas(gcf,[plotPath 'greenGauss_red-green_' num2str(N) '.jpg']);
            end
            hold off

            bar3(0:255,gbGreenHist','c');
            hold on
            greenDist_2D = double(greenDist(:,2:3));
            [greenObj_2D,isConverged] = EM(greenDist_2D,N);
            if isConverged
                greenY_2D = pdf(greenObj_2D,greenDist_2D);
                plot3(greenDist_2D(:,1),greenDist_2D(:,2),greenY_2D,'k*')
                xlabel('Green Intensity')
                ylabel('Blue Intensity')
                zlabel('Frequency')
                title(['Green-Blue Color Band of Green Buoy represented by ' num2str(N) ' Gaussians'])
                view(-45,10)
                saveas(gcf,[plotPath 'greenGauss_green-blue_' num2str(N) '.jpg']);
            end
            hold off
            
            bar3(0:255,rbGreenHist','m');
            hold on
            greenDist_2D = double(greenDist(:,[1,3]));
            [greenObj_2D,isConverged] = EM(greenDist_2D,N);
            if isConverged
                greenY_2D = pdf(greenObj_2D,greenDist_2D);
                plot3(greenDist_2D(:,1),greenDist_2D(:,2),greenY_2D,'k*')
                xlabel('Red Intensity')
                ylabel('Blue Intensity')
                zlabel('Frequency')
                title(['Red-Blue Color Band of Green Buoy represented by ' num2str(N) ' Gaussians'])
                view(-45,10)
                saveas(gcf,[plotPath 'greenGauss_red-blue_' num2str(N) '.jpg']);
            end
            hold off

            figure('units','normalized','outerposition',[0 0 1 1])
            bar3(0:255,rgRedHist','y');
            hold on
            redDist_2D = double(redDist(:,1:2));
            [redObj_2D,isConverged] = EM(redDist_2D,N);
            if isConverged
                redY_2D = pdf(redObj_2D,redDist_2D);
                plot3(redDist_2D(:,1),redDist_2D(:,2),redY_2D,'k*')
                xlabel('Red Intensity')
                ylabel('Green Intensity')
                zlabel('Frequency')
                title(['Red-Green Color Band of Red Buoy represented by ' num2str(N) ' Gaussians'])
                view(-45,10)
                saveas(gcf,[plotPath 'redGauss_red-green_' num2str(N) '.jpg']);
            end
            hold off

            bar3(0:255,gbRedHist','c');
            hold on
            redDist_2D = double(redDist(:,2:3));
            [redObj_2D,isConverged] = EM(redDist_2D,N);
            if isConverged
                redY_2D = pdf(redObj_2D,redDist_2D);
                plot3(redDist_2D(:,1),redDist_2D(:,2),redY_2D,'k*')
                xlabel('Green Intensity')
                ylabel('Blue Intensity')
                zlabel('Frequency')
                title(['Green-Blue Color Band of Red Buoy represented by ' num2str(N) ' Gaussians'])
                view(-45,10)
                saveas(gcf,[plotPath 'redGauss_green-blue_' num2str(N) '.jpg']);
            end
            hold off
            
            bar3(0:255,rbRedHist','m');
            hold on
            redDist_2D = double(redDist(:,[1,3]));
            [redObj_2D,isConverged] = EM(redDist_2D,N);
            if isConverged
                redY_2D = pdf(redObj_2D,redDist_2D);
                plot3(redDist_2D(:,1),redDist_2D(:,2),redY_2D,'k*')
                xlabel('Red Intensity')
                ylabel('Blue Intensity')
                zlabel('Frequency')
                title(['Red-Blue Color Band of Red Buoy represented by ' num2str(N) ' Gaussians'])
                view(-45,10)
                saveas(gcf,[plotPath 'redGauss_red-blue_' num2str(N) '.jpg']);
            end
            hold off

            figure('units','normalized','outerposition',[0 0 1 1])
            bar3(0:255,rgYellowHist','y');
            hold on
            yellowDist_2D = double(yellowDist(:,1:2));
            [yellowObj_2D,isConverged] = EM(yellowDist_2D,N);
            if isConverged
                yellowY_2D = pdf(yellowObj_2D,yellowDist_2D);
                plot3(yellowDist_2D(:,1),yellowDist_2D(:,2),yellowY_2D,'k*')
                xlabel('Red Intensity')
                ylabel('Green Intensity')
                zlabel('Frequency')
                title(['Red-Green Color Band of Yellow Buoy represented by ' num2str(N) ' Gaussians'])
                view(-45,10)
                saveas(gcf,[plotPath 'yellowGauss_red-green_' num2str(N) '.jpg']);
            end
            hold off

            bar3(0:255,gbYellowHist','c');
            hold on
            yellowDist_2D = double(yellowDist(:,2:3));
            [yellowObj_2D,isConverged] = EM(yellowDist_2D,N);
            if isConverged
                yellowY_2D = pdf(yellowObj_2D,yellowDist_2D);
                plot3(yellowDist_2D(:,1),yellowDist_2D(:,2),yellowY_2D,'k*')
                xlabel('Green Intensity')
                ylabel('Blue Intensity')
                zlabel('Frequency')
                title(['Green-Blue Color Band of Yellow Buoy represented by ' num2str(N) ' Gaussians'])
                view(-45,10)
                saveas(gcf,[plotPath 'yellowGauss_green-blue_' num2str(N) '.jpg']);
            end
            hold off
            
            bar3(0:255,rbYellowHist','m');
            hold on
            yellowDist_2D = double(yellowDist(:,[1,3]));
            [yellowObj_2D,isConverged] = EM(yellowDist_2D,N);
            if isConverged
                yellowY_2D = pdf(yellowObj_2D,yellowDist_2D);
                plot3(yellowDist_2D(:,1),yellowDist_2D(:,2),yellowY_2D,'k*')
                xlabel('Red Intensity')
                ylabel('Blue Intensity')
                zlabel('Frequency')
                title(['Red-Blue Color Band of Yellow Buoy represented by ' num2str(N) ' Gaussians'])
                view(-45,10)
                saveas(gcf,[plotPath 'yellowGauss_red-blue_' num2str(N) '.jpg']);
            end
            hold off
        end
    end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code implements the first question of the mid-term.
% 
% Submitted by: Ashwin Goyal (UID - 115526297)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define the coordinates of the cube
X = [0 0 0; 1 0 0; 1 1 0; 0 1 0; 0 1 1; 0 0 1; 1 0 1; 1 1 1]';

% Define calibration matrix
K = [800, 0, 250; 0, 800, 250; 0, 0, 1];

% Define translation
T = [0; 0; 5]; % assuming a 5 units translation along Z

% Define rotation
tx = 20*pi/180;
R = [1 0 0; 0 cos(tx) -sin(tx); 0 sin(tx) cos(tx)];

% Get the pixel coordinates
x = project(X,R,T,K);

% Plot the image
f = imshow(ones(uint16(max(max(x,[],2))+50)));
axis on
hold on
plot(x(1,1:4),x(2,1:4),'r-');
plot(x(1,[1,4]),x(2,[1,4]),'r-');
plot(x(1,5:8),x(2,5:8),'g-');
plot(x(1,[5,8]),x(2,[5,8]),'g-');
plot(x(1,[1,6]),x(2,[1,6]),'b-');
plot(x(1,[2,7]),x(2,[2,7]),'b-');
plot(x(1,[3,8]),x(2,[3,8]),'b-');
plot(x(1,[4,5]),x(2,[4,5]),'b-');
hold off

% Save the figure
saveas(gca,'../output/Image_of_Cube.jpg')
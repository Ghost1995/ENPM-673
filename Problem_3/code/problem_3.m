%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code implements the third question of the mid-term.
% 
% Submitted by: Ashwin Goyal (UID - 115526297)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define rotation matrix for rotation about X-axis
tx = pi/4;
Rx = [1 0 0; 0 cos(tx) -sin(tx); 0 sin(tx) cos(tx)];

% Define rotation matrix for rotation about Z-axis
tz = pi/6;
Rz = [cos(tz) -sin(tz) 0; sin(tz) cos(tz) 0; 0 0 1];

% Get equivalent rotation matrix for rotations around a stationary frame
R_fixed = Rx*Rz;
save('../output/R_fixed.mat','R_fixed');

% Get equivalent rotation matrix for rotations around a moving frame
R_moving = Rz*Rx;
save('../output/R_moving.mat','R_moving');
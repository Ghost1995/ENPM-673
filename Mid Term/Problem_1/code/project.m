%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code generates pixel coordinates of projected points in a image.
% 
% Input:
%   X --> 3 × n vector of the coordinates of 3D points
%   R --> rotation matrix for the camera coordinate frame with respect to
%         the world coordinate frame
%   T --> translation matrix for the camera coordinate frame with respect 
%         to the world coordinate frame
%   K --> matrix of intrinsic image parameters
% Output:
%   x --> pixel coordinates of the projected points in the image
% 
% Submitted by: Ashwin Goyal (UID - 115526297)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function x = project(X, R, T, K)

% Compute projective pixel coordinates
x = K*(R*X + T*ones(1,size(X,2)));

% Compute final pixel coordinates
x = x(1:2,:)./x(3,:);

end
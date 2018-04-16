%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code implements the fourth question of the mid-term.
% 
% Submitted by: Ashwin Goyal (UID - 115526297)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

I_left = imread('../input/tsukuba_l.png');
I_right = imread('../input/tsukuba_r.png');

% Compute disparity using SSD
disp('Solving Problem 4 using SSD');
disparity_map_ssd(I_left,I_right,7);

% Compute disparity using SAD
disp('Solving Problem 4 using SAD');
disparity_map_sad(I_left,I_right,7);

% Compute disparity using SSD with sub-pixel refinement
disp('Solving Problem 4 using SSD with Sub-fixel Refinement');
disparity_map_ssd_subpixel(I_left,I_right,3);

% Compute disparity using SAD with sub-pixel refinement
disp('Solving Problem 4 using SAD with Sub-fixel Refinement');
disparity_map_sad_subpixel(I_left,I_right,9);

% Compute disparity using median filtering after SSD
disp('Solving Problem 4 using Median Filtering after SSD');
disparity_map_median(I_left,I_right,7,9);

% Compute disparity using energy minimization after SSD
disp('Solving Problem 4 using Energy Minimization');
disparity_map_energy_minimization(I_left,I_right,3);

% Compute disparity using adaptive window after SSD
disp('Solving Problem 4 using Adaptive Window');
disparity_map_ssd_adaptive_window(I_left,I_right,9);
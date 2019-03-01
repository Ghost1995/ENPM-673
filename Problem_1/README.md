# Problem 1

**Problem Statement**

*Given:* A unit cube placed in the origin of the world coordinate system; camera is translated along z-axis by
some amount and rotated around x-axis by 20 degrees; calibration matrix K = [800, 0, 250; 0, 800, 250; 0, 0, 1]

*To Code:* A MATLAB function x = project(X, R, T, K) which takes image coordinates of 3-D points in
the world coordinate frame as input and generates pixel coordinates of the projected points in the image,
assuming that (R, T) is the displacement of camera coordinate frame with respect to world frame, K is the
matrix of intrinsic image parameters, and X is a 3 x n vector of the coordinates of 3D points.

**Running the Code**

Simply run *problem_1.m* in the code folder.

**Result**

The result is an image of the cube when camera is translated by 5 units.

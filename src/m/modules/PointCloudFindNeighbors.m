function [flags] = PointCloudFindNeighbors(x,y,mindistance,multithread);
%POINTCLOUDFINDNEIGHBORS - Flag points that are too near one another, within an array of point coordinates
%
%	   Usage:
%	      [flags]=PointCloudFindNeighbors(x,y,mindistance,multithread);
%
%	      x,y: list of points.
%	      mindistance: minimum distance that should exist between points in the cloud.
%	      multithread: run multithreaded or not. with multithreads, flags can get 1 and 2 values in duplicates.
%	      flags: array of flags (flag==1 means point is within mindistance of another point)

% Check usage
if nargin~=4
	help PointCloudFindNeighbors
	error('Wrong usage (see above)');
end

% Call mex module
[flags] = PointCloudFindNeighbors_matlab(x,y,mindistance,multithread);

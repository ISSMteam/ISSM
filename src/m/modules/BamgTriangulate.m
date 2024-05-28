function index = BamgTriangulate(x,y);
%BAMGTRIANGULATE - Delaunay Triangulation of a list of points
%
%   Usage:
%      index=BamgTriangulate(x,y)
%
%      index: index of the triangulation
%      x,y: coordinates of the nodes

%Check usage
if nargin~=2
	help BamgTriangulate
	error('Wrong usage (see above)');
end

%Call mex module
[index] = BamgTriangulate_matlab(x,y);

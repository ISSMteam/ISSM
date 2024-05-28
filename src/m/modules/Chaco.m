function [assgn] = Chaco(A,vwgts,ewgts,x,y,z,options,nparts,goal);
%CHACO
%
%   Usage:
%      [assgn] = Chaco(A,vwgts,ewgts,x,y,z,options,nparts,goal);
%
%   A:			Input adjacency matrix
%   vwgts:		weights for all vertices
%   ewgts:		weights for all edges
%   x,y,z:		coordinates for inertial method
%   options:		architecture and partitioning options
%   nparts:		number of parts options
%   goal:		desired set sizes

% Check usage
if nargin~=9
	help Chaco
	error('Wrong usage (see above)');
end

% Call mex module
[assgn] = Chaco_matlab(A,vwgts,ewgts,x,y,z,options,nparts,goal);

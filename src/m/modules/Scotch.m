function [maptab] = Scotch(varargin);
%SCOTCH - Scotch partitioner
%
%   Usage:
%      [maptab]=Scotch(adjmat,vertlb,vertwt,edgewt,archtyp,archpar,Scotch-specific parameters);

% Check usage
if nargin<6
	help Scotch
	error('Wrong usage (see above)');
end

% Call mex module
[maptab]=Scotch_matlab(varargin{:});

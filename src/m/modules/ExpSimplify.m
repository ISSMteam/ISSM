function ExpSimplify(varargin);
%EXPSIMPLIFY - Simplify an Exp contour
%
%   Usage:
%      ExpSimplify(expfile,tol);
%
%   expfile:	name of the expfile
%   tol:	tolerance (maximal euclidean distance allowed between the new line and a vertex)
%   min:	minimum number of vertices to save contours in an exp file (default is 3) [OPTIONAL]
%
%   Example:
%      ExpSimplify('file.exp',100);
%      ExpSimplify('file.exp',100,'min','4');

% Check usage
if nargin~=2 && nargin~=4
	help ExpSimplify
	error('Wrong usage (see above)');
end

% Call mex module
switch nargin
	case 2
		ExpSimplify_matlab(varargin{1},varargin{2});
	case 4
		ExpSimplify_matlab(varargin{1},varargin{2},varargin{3},varargin{4});
	otherwise
		error('Exp2Kml not supported yet');
end


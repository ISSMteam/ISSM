function data_prime=InterpFromMesh2d(varargin);
%INTERPFROMMESH2D
%
%   Usage:
%      data_prime=InterpFromMesh2d(index,x,y,data,x_prime,y_prime);
%      OR
%      data_prime=InterpFromMesh2d(index,x,y,data,x_prime,y_prime,default_value);
%      OR
%      data_prime=InterpFromMesh2d(index,x,y,data,x_prime,y_prime,default_value,contourname);
%
%   index:	index of the mesh where data is defined
%   x,y:	coordinates of the nodes where data is defined
%   data:	vector holding the data to be interpolated onto the points
%   x_prime,y_prime:	coordinates of the mesh vertices onto which we interpolate
%   default_value:	a scalar or vector of size length(x_prime)
%   contourname:	linear interpolation will happen on all x_interp,y_interp inside the contour,
%      default value will be adopted on the rest of the mesh.
%
%   data_prime:	vector of prime interpolated data

% Check usage
if nargin~=6 && nargin~=7 && nargin~=8
	help InterpFromMesh2d
	error('Wrong usage (see above)');
end

% Call mex module
switch nargin
	case 6
		data_prime=InterpFromMesh2d_matlab(varargin{1},varargin{2},varargin{3},varargin{4},varargin{5},varargin{6});
	case 7
		data_prime=InterpFromMesh2d_matlab(varargin{1},varargin{2},varargin{3},varargin{4},varargin{5},varargin{6},varargin{7});
	case 8 
		data_prime=InterpFromMesh2d_matlab(varargin{1},varargin{2},varargin{3},varargin{4},varargin{5},varargin{6},varargin{7},varargin{8});
	otherwise
		% NOTE: Should never get here because of previous check
		error('InterpFromMesh2d not supported');
end




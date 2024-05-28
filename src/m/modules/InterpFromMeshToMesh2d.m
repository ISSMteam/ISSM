function data_interp = InterpFromMeshToMesh2d(varargin);
%INTERPFROMMESHTOMESH2D - Interpolation from a 2d triangular mesh onto a list of points
%
%   Usage:
%      data_interp=InterpFromMeshToMesh2d(index,x,y,data,x_interp,y_interp);
%      or data_interp=InterpFromMeshToMesh2d(index,x,y,data,x_interp,y_interp,OPTIONS);
%
%   index             : index of the mesh where data is defined (e.g. md.mesh.elements)
%   x,y               : coordinates of the nodes where data is defined
%   data              : matrix holding the data to be interpolated onto the mesh (one column per field)
%   x_interp,y_interp : coordinates of the points onto which we interpolate
%   data_interp       : vector of mesh interpolated data
%   Available options:
%      default:	default value if point is outsite of triangulation (instead of linear interpolation)
%
%   Example:
%      load('temperature.mat');
%      md.initialization.temperature=InterpFromMeshToMesh2d(index,x,y,temperature,md.mesh.x,md.mesh.y);
%      md.initialization.temperature=InterpFromMeshToMesh2d(index,x,y,temperature,md.mesh.x,md.mesh.y,'default',253);

% Check usage
if nargin~=6 && nargin~=8
	help InterpFromMeshToMesh2d
	error('Wrong usage (see above)');
end

% Call mex module
data_interp=InterpFromMeshToMesh2d_matlab(varargin{:});

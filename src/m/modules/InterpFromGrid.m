function dataout = InterpFromGrid(x, y, data, x_interp,y_interp, method)
%INTERPFROMGRID - Interpolation from a grid onto a list of points (faster than InterpFromGridToMesh)
%
%   Usage:
%      dataout = InterpFromGrid2(x, y, data, x_interp,y_interp)
%      dataout = InterpFromGrid2(x, y, data, x_interp,y_interp, method)
%
%   data:    matrix holding the data to be interpolated onto the mesh
%   x,y:     coordinates of matrix data (x and y must be in increasing order)
%   x_interp,y_interp: coordinates of the points onto which we interpolate
%   method:  interpolation method ('linear'/default, 'nearest', 'cubic')
%
%	 Example:
%      md.inversion.vx_obs=InterpFromGrid(x, y, vx, md.mesh.x,md.mesh.y);

% Call mex module
if nargin==5
	dataout = InterpFromGrid_matlab(x, y, data, x_interp, y_interp);
elseif nargin==6
	dataout = InterpFromGrid_matlab(x, y, data, x_interp, y_interp, method);
else
	help InterpFromGrid
	error('Wrong usage (see above)');
end

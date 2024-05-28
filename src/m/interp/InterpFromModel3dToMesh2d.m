function data_prime = InterpFromModel3dToMesh2d(md,data,x_prime,y_prime,sigma,default_value);
%INTERPFROMMODEL3DTOMESH2D - Interpolation from a 3d hexahedron mesh onto a list of 2d points
%
%   Usage:
%      md:  3d model holding the data to be interpolated onto the 2d mesh
%      data:	matrix holding the data to be interpolated onto the 2d mesh
%      x_prime,y_prime:	coordinates of the points onto which we interpolate
%      sigma:  scaled vertical coordinates from which the data will be interpolated (from base=0; from surface=1, NaN=vertical average of data)
%      default_value:	default value if no data is found (holes)
%      data_prime:	vector of mesh interpolated data
%
%   Example:
%      load('temperature.mat');
%
%      % interpolating the temperature from the base of a 3d model:
%      md.initialization.temperature=InterpFromModel3dToMesh2d(md3d,temperature,md.mesh.x,md.mesh.y,0,253);
%
%      % interpolating the temperature from the surface of a 3d model:
%      md.initialization.temperature=InterpFromModel3dToMesh2d(md3d,temperature,md.mesh.x,md.mesh.y,1,253);
%
%      % averaging the temperature over the vertical layers and then interpolating onto the 2d mesh:
%      md.initialization.temperature=InterpFromModel3dToMesh2d(md3d,temperature,md.mesh.x,md.mesh.y,NaN,253);

% Check usage
if nargin~=6
	help InterpFromModel3dToMesh2d
	error('Wrong usage (see above)');
end

if dimension(md.mesh)~=3
	error('Model should be 3d');
end

if (length(data)~=md.mesh.numberofelements & length(data)~=md.mesh.numberofvertices),
   error('Data not supported yet');
end

if sigma<0 | sigma>1
	help InterpFromModel3dToMesh2d
	error('Wrong value for sigma. It should be between 0 and 1, or NaN (see above)');
end

if length(x_prime)~=length(y_prime)
	error('x and y should have the same size')
end

% First, check if a vertical average should be performed. If yes, perform a interpolation from a 2d mesh  
if isnan(sigma),
	% average data and then interpolate onto the 2d mesh
	averaged_data = DepthAverage(md,data);
	data_prime = InterpFromMeshToMesh2d(md.mesh.elements2d,md.mesh.x2d,md.mesh.y2d,averaged_data,x_prime,y_prime,'default',default_value);
else
	% Ok, perform an interpolation from a 3d mesh into a 2d mesh 
	% Scaling the z coordinate (3d mesh)
	alpha = (md.mesh.z-md.geometry.base)./md.geometry.thickness;

	if  alpha<0 | alpha>1
		error('Wrong value for alpha. Check the geometry of your model');
	end

	% Building the z coordinate for the 2d mesh
	z_prime = sigma*ones(size(x_prime));

	% Adjusting such that the 2d mesh is inside the 3d mesh
	if sigma==0,
		z_prime = z_prime+eps;
	elseif sigma==1,
		z_prime = z_prime-eps;
	end

	% Now, call the 3d interpolation
	data_prime = InterpFromMeshToMesh3d(md.mesh.elements,md.mesh.x,md.mesh.y,alpha,data,x_prime,y_prime,z_prime,default_value);

end

function data_prime = InterpFromModelToModel3d(md1,data,md2,default_value);
%INTERPFROMMODELTOMODEL3D - Interpolation from a 3d hexahedron mesh onto another 3d hexahedron mesh
%
%   Usage:
%      md1:  3d model holding the data to be interpolated onto another 3d model
%      data:  matrix holding the data to be interpolated onto another 3d mesh
%      md2:  3d model for which the data will be interpolated
%      default_value:  default value if no data is found (holes)
%      data_prime:  vector of 3d mesh interpolated data
%
%   Example:
%
%      % interpolating the temperature from 3d mesh model:
%      md2.initialization.temperature = InterpFromModelToModel3d(md1,md1.results.ThermalSolution.Temperature,md2,253);
%

% Check usage
if nargin~=4
   help InterpFromModelToModel3d
   error('Wrong usage (see above)');
end

if (length(data)~=md1.mesh.numberofelements & length(data)~=md1.mesh.numberofvertices),
   error('Data not supported yet');
end

if (dimension(md1.mesh)~=3 | dimension(md2.mesh)~=3)
      error('Both models should be 3d');
end

% Scaling the vertical coordinates:
sigma1 = (md1.mesh.z-md1.geometry.base)./md1.geometry.thickness;
sigma2 = (md2.mesh.z-md2.geometry.base)./md2.geometry.thickness;

% Adjusting sigma2 such that mesh 2 is inside mesh 1
pos = find(sigma2==0);
sigma2(pos) = sigma2(pos)+eps; 
pos = find(sigma2==1);
sigma2(pos) = sigma2(pos)-eps;

% Now, perform the interpolation
data_prime = InterpFromMeshToMesh3d(md1.mesh.elements,md1.mesh.x,md1.mesh.y,sigma1,data,md2.mesh.x,md2.mesh.y,sigma2,default_value); 

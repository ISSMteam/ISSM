function projection_value=project2d(md3d,value,layer)
%PROJECT2D - returns the value of a field for a given layer of the mesh
%
%   returns the value of a vector for a given layer from extruded mesh onto the 2d mesh 
%   used to do the extrusion. This function is used to compare values between different
%   layers of a 3d mesh.
%
%   Usage:
%      projection_value=project2d(md3d,value,layer)
%
%   Example:
%      vel2=project2d(md3d,md3d.initialization.vel,2);
%      returns the velocity of the second layer (1 is the base)

%some checks on list of arguments
if ((nargin~=3) ),
	help project2d
	error('project2d error message');
end

if ~strcmp(md3d.mesh.domaintype,'3D');
	error('wrong model type ... should be ''3d''');
end

if ((layer<1) | (layer>md3d.mesh.numberoflayers)),
	error(['layer must be between 1 and ' num2str(md3d.mesh.numberoflayers)]);
end

if numel(value)==1
	projection_value=value;
elseif size(value,1)==md3d.mesh.numberofvertices,
	projection_value=value((layer-1)*md3d.mesh.numberofvertices2d+1:layer*md3d.mesh.numberofvertices2d,:);
elseif size(value,1)==md3d.mesh.numberofvertices+1
	projection_value=[value((layer-1)*md3d.mesh.numberofvertices2d+1:layer*md3d.mesh.numberofvertices2d,:); value(end,:)];
elseif size(value,1)==md3d.mesh.numberofelements
	projection_value=value((layer-1)*md3d.mesh.numberofelements2d+1:layer*md3d.mesh.numberofelements2d,:);
else
	error('Dimensions not supported yet');
end

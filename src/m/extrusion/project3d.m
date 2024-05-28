function projected_vector=project3d(md,varargin)
%PROJECT3D - vertically project a vector from 2d mesh
%
%   vertically project a vector from 2d mesh (split in noncoll and coll areas) into a 3d mesh.
%   This vector can be a node vector of size (md.mesh.numberofvertices2d,N/A) or an 
%   element vector of size (md.mesh.numberofelements2d,N/A). 
%   arguments: 
%      'vector': 2d vector
%      'type': 'element' or 'node' or 'poly'
%   options: 
%      'layer' a layer number where vector should keep its values. If not specified, all layers adopt the 
%             value of the 2d vector.
%      'padding': default to 0 (value adopted by other 3d layers not being projected0
%		 'degree': degree of polynomials when extrude from bottom to the top
%
%   Egs:
%      extruded_vector=project3d(md,'vector',vector2d,'type','node','layer',1,'padding',NaN);
%      extruded_vector=project3d(md,'vector',vector2d,'type','element','padding',0);
%      extruded_vector=project3d(md,'vector',vector2d,'type','node');

%some regular checks
if nargin==0,
	help project3d
	error('bad usage');
end
if ~strcmp(elementtype(md.mesh),'Penta')
	error('input model is not 3d');
end

%retrieve parameters from options.
options      = pairoptions(varargin{:});
vector2d     = getfieldvalue(options,'vector');     %mandatory
type         = getfieldvalue(options,'type');       %mandatory
layer        = getfieldvalue(options,'layer',0);    %optional (do all layers otherwise)
paddingvalue = getfieldvalue(options,'padding',0);  %0 by default
polyexponent = getfieldvalue(options,'degree',0);   %0 by default, 0-degree polynomial

if length(vector2d)==1,
	projected_vector=vector2d;

elseif strcmpi(type,'node'),
	%Initialize 3d vector
	if size(vector2d,1)==md.mesh.numberofvertices2d
		projected_vector=paddingvalue*ones(md.mesh.numberofvertices, size(vector2d,2));
	elseif size(vector2d,1)==md.mesh.numberofvertices2d+1
		projected_vector=paddingvalue*ones(md.mesh.numberofvertices+1,size(vector2d,2));
		projected_vector(end,:)=vector2d(end,:);
		vector2d=vector2d(1:end-1,:);
	else
		error('vector length not supported')
	end

	%Fill in
	if layer==0,
		for i=1:md.mesh.numberoflayers,
			projected_vector(((i-1)*md.mesh.numberofvertices2d+1):(i*md.mesh.numberofvertices2d),:)=vector2d;
		end
	else
		projected_vector(((layer-1)*md.mesh.numberofvertices2d+1):(layer*md.mesh.numberofvertices2d),:)=vector2d;
	end

elseif strcmpi(type,'element'),

	%Initialize 3d vector
	if size(vector2d,1)==md.mesh.numberofelements2d
		projected_vector=paddingvalue*ones(md.mesh.numberofelements,  size(vector2d,2));
	elseif size(vector2d,1)==md.mesh.numberofelements2d+1
		projected_vector=paddingvalue*ones(md.mesh.numberofelements+1,size(vector2d,2));
		projected_vector(end,:)=vector2d(end,:);
		vector2d=vector2d(1:end-1,:);
	else
		error('vector length not supported')
	end

	%Fill in
	if layer==0,
		for i=1:(md.mesh.numberoflayers-1),
			projected_vector( ((i-1)*md.mesh.numberofelements2d+1):(i*md.mesh.numberofelements2d),:)=vector2d;
		end
	else
		projected_vector( ((layer-1)*md.mesh.numberofelements2d+1):(layer*md.mesh.numberofelements2d),:)=vector2d;
	end

elseif strcmpi(type,'poly'), % interpolate values from 0 to 1 with a polynomial degree n
	%Initialize 3d vector
	if size(vector2d,1)==md.mesh.numberofvertices2d
		projected_vector=paddingvalue*ones(md.mesh.numberofvertices, size(vector2d,2));
	elseif size(vector2d,1)==md.mesh.numberofvertices2d+1
		projected_vector=paddingvalue*ones(md.mesh.numberofvertices+1,size(vector2d,2));
		projected_vector(end,:)=vector2d(end,:);
		vector2d=vector2d(1:end-1,:);
	else
		error('vector length not supported')
	end

	polycoeff = [0:1./(md.mesh.numberoflayers-1):1];

	%Fill in
	if layer==0,
		for i=1:md.mesh.numberoflayers,
			projected_vector(((i-1)*md.mesh.numberofvertices2d+1):(i*md.mesh.numberofvertices2d),:)=vector2d*(1-(1-polycoeff(i)).^polyexponent);
		end
	else
		projected_vector(((layer-1)*md.mesh.numberofvertices2d+1):(layer*md.mesh.numberofvertices2d),:)=vector2d*(1-(1-polycoeff(layer)).^polyexponent);
	end

else
	error('project3d error message: unknown projection type');
end

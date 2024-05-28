function average=averaging(md,data,iterations,varargin)
%AVERAGING - smooths the input over the mesh
%
%   This routine takes a list over the elements or the nodes in input
%   and return a list over the nodes.
%   For each iterations it computes the average over each element (average 
%   of the vertices values) and then computes the average over each node
%   by taking the average of the element around a node weighted by the
%   elements volume
%   For 3d mesh, a last argument can be added to specify the layer to be averaged on.
%
%   Usage:
%      smoothdata=averaging(md,data,iterations)
%      smoothdata=averaging(md,data,iterations,layer)
%
%   Examples:
%      velsmoothed=averaging(md,md.initialization.vel,4);
%      pressure=averaging(md,md.initialization.pressure,0);
%      temperature=averaging(md,md.initialization.temperature,1,1);

if ((nargin~=4) & (nargin~=3)),
	error('averaging error message: wrong number of arguments');
end
if (length(data)~=md.mesh.numberofelements & length(data)~=md.mesh.numberofvertices),
	error('averaging error message: data not supported yet');
end
if dimension(md.mesh)==3 & nargin==4,
	if varargin{1}<=0 | varargin{1}>md.mesh.numberoflayers,
		error('layer should be between 1 and md.mesh.numberoflayers');
	end
	layer=varargin{1};
else
	layer=0;
end

%initialization
if layer==0,
	weights=zeros(md.mesh.numberofvertices,1);
	data=data(:);
else 
	weights=zeros(md.mesh.numberofvertices2d,1);
	data=data((layer-1)*md.mesh.numberofvertices2d+1:layer*md.mesh.numberofvertices2d,:);
end

%load some variables (it is much faster if the variables are loaded from md once for all)
if layer==0,
	index=md.mesh.elements;
	numberofnodes=md.mesh.numberofvertices;
	numberofelements=md.mesh.numberofelements;
else
	index=md.mesh.elements2d;
	numberofnodes=md.mesh.numberofvertices2d;
	numberofelements=md.mesh.numberofelements2d;
end

%build some variables
line=index(:);
if dimension(md.mesh)==3 & layer==0,
	rep=6;
	areas=GetAreas(index,md.mesh.x,md.mesh.y,md.mesh.z);
elseif dimension(md.mesh)==2,
	rep=3;
	areas=GetAreas(index,md.mesh.x,md.mesh.y);
else
	rep=3;
	if isa(md.mesh,'mesh3dsurface'),
		areas=GetAreas3DTria(md.mesh.elements,md.mesh.x,md.mesh.y,md.mesh.z);
	else
		areas=GetAreas(index,md.mesh.x,md.mesh.y);
	end
end
summation=1/rep*ones(rep,1);
linesize=rep*numberofelements;

%update weights that hold the volume of all the element holding the node i
weights=sparse(line,ones(linesize,1),repmat(areas,rep,1),numberofnodes,1);

%initialization
if length(data)==numberofelements
	average_node=sparse(line,ones(linesize,1),repmat(areas.*data,rep,1),numberofnodes,1);
	average_node=average_node./weights;
else
	average_node=data;
end

%loop over iteration
for i=1:iterations
	average_el=average_node(index)*summation;
	average_node=sparse(line,ones(linesize,1),repmat(areas.*average_el,rep,1),numberofnodes,1);
	average_node=average_node./weights;
end

%return output as a full matrix (C code does not like sparse matrices)
average=full(average_node);

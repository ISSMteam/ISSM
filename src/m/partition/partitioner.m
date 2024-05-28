function [partitionvector,md]=partitioner(md,varargin)
%PARTITIONER - partition mesh 
%
%   List of options to partitioner: 
%
%   package: 'chaco', 'metis' or 'scotch'
%   npart: number of partitions.
%   weighting: 'on' or 'off': default off
%   section:  1 by defaults(1=bisection, 2=quadrisection, 3=octasection)
%   recomputeadjacency:  'on' by default (set to 'off' to compute existing one)
%   type: 'node' or 'element' partition vector (default to 'node') 
%   Output: partitionvector: the partition vector
%   
%   Usage:
%      partitionvector=partitioner(md,'package','chaco','npart',100,'weighting','on');
%      [partitionvector,md]=partitioner(md,'package','chaco','npart',100,'weighting','on');
%

%get options: 
options=pairoptions(varargin{:});

%set defaults
options=addfielddefault(options,'package','chaco');
options=addfielddefault(options,'npart',10);
options=addfielddefault(options,'weighting','on');
options=addfielddefault(options,'section',1);
options=addfielddefault(options,'recomputeadjacency','on');
options=addfielddefault(options,'type','node');

%get package: 
package=getfieldvalue(options,'package');
npart=getfieldvalue(options,'npart');
recomputeadjacency=getfieldvalue(options,'recomputeadjacency');
vectortype=getfieldvalue(options,'type');

if(dimension(md.mesh)==3),
	%partitioning essentially happens in 2D. So partition in 2D, then 
	%extrude the partition vector vertically. 
	md3d=md; %save for later
	md.mesh.elements=md.mesh.elements2d;
	md.mesh.x=md.mesh.x2d;
	md.mesh.y=md.mesh.y2d;
	md.mesh.numberofvertices=md.mesh.numberofvertices2d;
	md.mesh.numberofelements=md.mesh.numberofelements2d;
	md.qmu.vertex_weight=[];
	md.mesh.vertexconnectivity=[];
	recomputeadjacency='on';
end

%adjacency matrix if needed:
if strcmpi(recomputeadjacency,'on'),
	md=adjacency(md);
else
	disp('skipping adjacency matrix computation as requested in the options');
end

if strcmpi(package,'chaco'),

	if strcmpi(vectortype,'element')
		error(['partitioner error message: package ' package ' does not allow element partitions.']);
	else

		%  default method (from chaco.m)
		method=[1 1 0 0 1 1 50 0 .001 7654321]';
		method(1)=3;    %  global method (3=inertial (geometric))
		method(3)=0;    %  vertex weights (0=off, 1=on)

		%specify bisection
		method(6)=getfieldvalue(options,'section');%  ndims (1=bisection, 2=quadrisection, 3=octasection)

		%are we using weights? 
		if strcmpi(getfieldvalue(options,'weighting'),'on'),
			weights=floor(md.qmu.vertex_weight/min(md.qmu.vertex_weight));
			method(3)=1;
		else 
			weights=[];
		end

		%  partition into nparts
		if isa(md.mesh,'mesh2d'),
			part=Chaco(md.qmu.adjacency,weights,[],md.mesh.x,md.mesh.y,zeros(md.mesh.numberofvertices,1),method,npart,[])'+1; %index partitions from 1 up. like metis.
		else
			part=Chaco(md.qmu.adjacency,weights,[],md.mesh.x, md.mesh.y,md.mesh.z,method,npart,[])'+1; %index partitions from 1 up. like metis.
		end

	end

elseif strcmpi(package,'scotch'),

	if strcmpi(vectortype,'element')
		error(['partitioner error message: package ' package ' does not allow element partitions.']);
	else
		%are we using weights? 
		if strcmpi(getfieldvalue(options,'weighting'),'on'),
			weights=floor(md.qmu.vertex_weight/min(md.qmu.vertex_weight));
		else
			weights=[];
		end
		maptab=Scotch(md.qmu.adjacency,[],weights,[],'cmplt',[npart]);

		part=maptab(:,2)+1;%index partitions from 1 up. like metis.
	end

elseif strcmpi(package,'linear'),

	if strcmpi(vectortype,'element')
		part=1:1:md.mesh.numberofelements;
		disp('Linear partitioner requesting partitions on elements');
	else
		part=1:1:md.mesh.numberofvertices;
	end

elseif strcmpi(package,'metis'),

	if strcmpi(vectortype,'element')
		error(['partitioner error message: package ' package ' does not allow element partitions.']);
	else
		[element_partitioning,part]=MeshPartition(md,md.qmu.numberofpartitions);
	end

else

	error(['partitioner error message: could not find ' package ' partitioner']);
	help partitioner
end

%extrude if we are in 3D:
if dimension(md.mesh)==3,
	md3d.qmu.vertex_weight=md.qmu.vertex_weight;
	md3d.qmu.adjacency=md.qmu.adjacency;
	md=md3d;
	if strcmpi(vectortype,'element')
		part=project3d(md,'vector',part','type','element');
	else
		part=project3d(md,'vector',part','type','node');
	end
end

if size(part,1)==1
	part=part';
end

%output 
partitionvector=part;

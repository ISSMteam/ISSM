function md=modelmerge3d(md1,md2,varargin)
%MODELMERGE  - merge two models by merging their meshes
%
%   Usage:
%      md=modelmerge(md1,md2);
	
	%process options: 
	options=pairoptions(varargin{:});
	
	tolerance=getfieldvalue(options,'tolerance',1e-4);
	
	md=model();

	%first ,copy md1 mesh into md.mesh to initialize, and additional classes:
	md.mesh=md1.mesh;
	md.private=md1.private;

	%some initialization: 
	elements1=md1.mesh.elements;
	x1=md1.mesh.x;
	y1=md1.mesh.y;
	z1=md1.mesh.z;
	nods1=md1.mesh.numberofvertices;
	nel1=md1.mesh.numberofelements;

	elements2=md2.mesh.elements;
	x2=md2.mesh.x;
	y2=md2.mesh.y;
	z2=md2.mesh.z;
	nods2=md2.mesh.numberofvertices;
	nel2=md2.mesh.numberofelements;

	%offset elements2 by nods1: 
	elements2=elements2+nods1;

	%go into the vertices on boundary of mesh 1 and figure out which ones are common with mesh2: 
	verticesonboundary=find(md1.mesh.vertexonboundary);

	for i=1:length(verticesonboundary),
		node1=verticesonboundary(i);
		xnode1=x1(node1);
		ynode1=y1(node1);
		znode1=z1(node1);

		%is there another node with these coordinates in mesh 2?
		ind=find(sqrt((x2-xnode1).^2+(y2-ynode1).^2+(z2-znode1).^2)<tolerance);
		if length(ind)>1,
			disp('should reduce the tolerance, several vertices picked up!');
		end
		if ~isempty(ind),
			x2(ind)=NaN;
			y2(ind)=NaN;
			z2(ind)=NaN;
			pos=find(elements2==(ind+nods1));
			elements2(pos)=node1;
		end
	end
	%go through elements2 and drop counter on each vertex that is above the x2 and y2 vertices being dropped:
	indices_nan=isnan(x2);
	while(~isempty(indices_nan)),
		% Use the index of the first instance of 'nan' value to remove that element from 'x2', 'y2', and 'z2'
		index_nan=indices_nan(1);
		pos=find(elements2>(index_nan+nods1));
		elements2(pos)=elements2(pos)-1;
		x2(index_nan)=[];
		y2(index_nan)=[];
		z2(index_nan)=[];

		% Check again in 'x2' for instances of 'nan'
		indices_nan=find(isnan(x2));
	end

	%merge elements: 
	elements=[elements1;elements2];

	%merge vertices: 
	x=[x1;x2]; 
	y=[y1;y2];
	z=[z1;z2];

	%output: 
	md.mesh.x=x;
	md.mesh.y=y;
	md.mesh.z=z;
	md.mesh.elements=elements;
	md.mesh.numberofvertices=length(x);
	md.mesh.numberofelements=size(elements,1);

	%connectivities: 
	md.mesh.vertexconnectivity=NodeConnectivity(md.mesh.elements,md.mesh.numberofvertices);
	md.mesh.elementconnectivity=ElementConnectivity(md.mesh.elements,md.mesh.vertexconnectivity);

	%find segments: 
	md.mesh.segments=findsegments(md);

	%vertex on boundary: 
	md.mesh.vertexonboundary=zeros(md.mesh.numberofvertices,1);
	md.mesh.vertexonboundary(md.mesh.segments(:,1:2))=1;

	%some checks: 
	if max(md.mesh.elements)>md.mesh.numberofvertices, 
		error('issue in modelmerge, one of the element ids is > number of vertices!');
	end

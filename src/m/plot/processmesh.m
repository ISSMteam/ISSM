function [x y z elements is2d isplanet]=processmesh(md,data,options)
%PROCESSMESH - process mesh to be plotted
%
%   Usage:
%      [x y z elements is2d]=processmesh(md,data,options)
%
%   See also: PLOTMODEL, PROCESSDATA

%some checks
if md.mesh.numberofvertices==0,
	error('plot error message: mesh is empty')
end
if md.mesh.numberofvertices==md.mesh.numberofelements
	error(['plot error message: the number of elements is the same as the number of nodes...']);
end

%special case for mesh 2dvertical
if strcmp(domaintype(md.mesh),'2Dvertical'),
	[x y z elements is2d isplanet] = processmesh(md.mesh,options);
	return;
end

%special case for mesh 3dsurface
if strcmp(domaintype(md.mesh),'3Dsurface'),
	[x y z elements is2d isplanet] = processmesh(md.mesh,options);
	if strcmpi(getfieldvalue(options,'coord','xy'),'latlon') | strcmpi(getfieldvalue(options,'coord','xy'),'latlong'),
		x0=md.mesh.long;
		y0=md.mesh.lat;
		%add row at lat=90 and lat=-90
		add=[(-180:.1:-1e-5)';(1e-5:.1:180)'];
		nadd=length(add);
		xextra=[add;add];
		yextra=[90*ones(nadd,1); -90*ones(nadd,1)];
		x=[x0;xextra];
		y=[y0;yextra];

		if strcmpi(getfieldvalue(options,'coordcent','atlantic'),'pacific'),
			pos=find(x>0);  x(pos)=-360+x(pos);
		end
		elements=delaunay(x,y);
		z=x; z(:)=0;
	end
	return;
end

if isprop(md.mesh,'elements2d'), elements2d=md.mesh.elements2d; end

if exist(options,'amr'),
	step = getfieldvalue(options,'amr');
	x = md.results.TransientSolution(step).MeshX;
	y = md.results.TransientSolution(step).MeshY;
	elements = md.results.TransientSolution(step).MeshElements;
else
	elements=md.mesh.elements;
	if ~strcmpi(getfieldvalue(options,'coord','xy'),'latlon') &  ~strcmpi(getfieldvalue(options,'coord','xy'),'latlong') ,
		x=md.mesh.x;
		if isprop(md.mesh,'x2d'), x2d=md.mesh.x2d; end
		y=md.mesh.y;
		if isprop(md.mesh,'y2d'), y2d=md.mesh.y2d; end
	else
		x=md.mesh.long;
		y=md.mesh.lat;
		if strcmpi(getfieldvalue(options,'coordcent','atlantic'),'pacific'),
			pos=find(x>0);  x(pos)-360+x(pos);
		end
	end
end

if isprop(md.mesh,'z'),
	z=md.mesh.z;
else
	z=zeros(size(x));
end
z=getfieldvalue(options,'z',z);
if ischar(z),
	z=md.(z);
end

%is it a 2d plot?
if md.mesh.dimension()==2 || getfieldvalue(options,'layer',0)>=1 || getfieldvalue(options,'depthaverage',0)
	is2d=1;
else
	is2d=0;
end

%layer projection? 
if getfieldvalue(options,'layer',0)>=1 || getfieldvalue(options,'depthaverage',0)
	if strcmpi(getfieldvalue(options,'coord','xy'),'latlon'),
		error('processmesh error message: cannot work with 3D meshes for now');
	end
	%we modify the mesh temporarily to a 2d mesh from which the 3d mesh was extruded
	x=x2d;
	y=y2d;
	z=zeros(size(x2d));
	elements=elements2d;
end

%units
if exist(options,'unit'),
	unit=getfieldvalue(options,'unit');
	x=x*unit;
	y=y*unit;
	z=z*unit;
end

%Quiver plot for elements?
if size(data,2)>1 && size(data,1)==size(elements,1)
	x = mean(x(elements),2);
	y = mean(y(elements),2);
end

if isa(md,'planet'),
	isplanet=1;
else
	isplanet=0;
end

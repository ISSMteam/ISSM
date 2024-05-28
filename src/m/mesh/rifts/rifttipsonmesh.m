function tips=rifttipsonmesh(md,riftoutline)
%RIFTTIPSONMESH: identify, using a rift outline, the nodes that are tips of 
%                rifts.

%read rift file according to its extension: 
[path,name,ext]=fileparts(riftoutline);
if strcmp(ext,'.exp'),
	rifts=expread(riftoutline);
elseif strcmp(ext,'.shp'),
	rifts=shpread(riftoutline);
else
	error(['bamg error message: file ' riftoutline ' format not supported (.shp or .exp)']);
end

tips=[];

for i=1:length(rifts),
	rift=rifts(i);

	x_tip=rift.x(1);
	y_tip=rift.y(1);

	index=find_point(md.mesh.x,md.mesh.y,x_tip,y_tip);
	tips(end+1)=index;

	x_tip=rift.x(end);
	y_tip=rift.y(end);

	index=find_point(md.mesh.x,md.mesh.y,x_tip,y_tip);
	tips(end+1)=index;

end

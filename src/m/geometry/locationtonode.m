function nodes=locationtonode(md,location,varargin)
%LOCATIONTONODE - find, given a string location (ex: 'LA', 'new york', the nearest node on a mesh3dsurface.
%
%   Usage:
%      node=locationnode(md,'LA');
%      nodes=locationnode(md,'LA',100); %option, specify a distance (km) around the location that will return nodes
%
%   See also: geoCode (in externalpackages), find_point

	if nargin==3,
		radius=varargin{1};
	else
		radius=0;
	end

	coords=geoCode(location,'osm');  
	latny=coords(1); longny=coords(2);
	node=find_point(md.mesh.lat,md.mesh.long,latny,longny);

	if radius>0,
		distance=sqrt( (md.mesh.x-md.mesh.x(node)).^2 + (md.mesh.y-md.mesh.y(node)).^2 + (md.mesh.z-md.mesh.z(node)).^2);
		nodes=find(distance<radius*1000);
	else
		nodes=node;
	end
end

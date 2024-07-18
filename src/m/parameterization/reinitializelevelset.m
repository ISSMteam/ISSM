function levelsetnew = reinitializelevelset(md,levelset)
%REINITIALIZELEVELSET - reinitialize levelset as a signed distance function
%
%   Usage:
%      levelsetnew = reinitializelevelset(md,levelset)
%
%   Example:
%      md.mask.ice_levelset = reinitializelevelset(md, md.mask.ice_levelset);

% if md is 3d, levelset should be projected on a 2d mesh 

if isempty(levelset), error('levelset provided is empty'); end
if dimension(md.mesh)==3,
   if length(levelset)~=md.mesh.numberofvertices2d, error('levelset provided should be specified at the 2d vertices of the mesh'); end
else
   if length(levelset)~=md.mesh.numberofvertices, error('levelset provided should be specified at the vertices of the mesh'); end
end

%First: extract segments
contours=isoline(md,levelset,'value',0);

%Now, make this a distance field (might not be closed)
levelsetnew=abs(ExpToLevelSet(md.mesh.x,md.mesh.y,contours)); % levelsetnew comes on the 3d vertices, if mesh is 3d

%Finally, change sign
pos=find(levelset<0); % if mesh is 3d, it refers to the vertices on the base
if dimension(md.mesh)==3
	for i=1:md.mesh.numberoflayers
		pos3d=pos+(i-1)*md.mesh.numberofvertices2d;
		levelsetnew(pos3d)=-levelsetnew(pos3d);
	end	
else
	levelsetnew(pos)=-levelsetnew(pos);
end

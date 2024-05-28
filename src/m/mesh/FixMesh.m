function  [index2 x2 y2 value2 newpos]=FixMesh(index,x,y,value)
% FIXMESH - FixMesh fix mesh with broken triangles, orphan vertices, etc ...
%
%   Usage: 
%      [index2 x2 y2 value2]=FixMesh(index,x,y,value)
%      where index,x,y is a delaunay triangulation, 
%      value is a field on the input triangulation, with values at the vertices
%      index2,x2,y2,value2 is the repaired triangulation, with new values on new vertices
%
%

%duplicate inputs
index2=index;
x2=x;
y2=y;
value2=value;
newpos=1:length(x);

%First, look for orphan vertices, and take them out.
flags=zeros(length(x2),1); flags(index2)=1;
orphans=find(flags==0);

while ~isempty(orphans),

	%take the first orphan, the lower numbered, and take it out
	orphan=orphans(1);

	%first x,y,value
	x2(orphan)=[];
	y2(orphan)=[];
	value2(orphan)=[];
	newpos(orphan)=[];

	%now, the index:
	pos=find(index2>orphan); index2(pos)=index2(pos)-1;

	%look again for orphans on new mesh
	flags=zeros(length(x2),1);flags(index2)=1;
	orphans=find(flags==0);
end

%Check all triangles are well oriented.
aires=GetAreas(index2,x2,y2);
pos=find(aires<0);
temp=index2(pos,1);
index2(pos,1)=index2(pos,2);
index2(pos,2)=temp;

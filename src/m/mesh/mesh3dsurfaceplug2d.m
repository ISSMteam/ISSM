function mh=mh3dsurfaceplug2d(mh,mh2,flags,segments,xsegs,ysegs,varargin)
%MESH3DSURFACEPLUG2D - plug 2d mesh into a 3D surface mesh
%
%   Usage:
%      mh=mesh3dsurfaceplug2d(mh,mh2);
%

	%First process options
	options=pairoptions(varargin{:});

	%Remove the elements that are flagged: 
	pos=find(flags); 
	mh.elements(pos,:)=[];
	mh.numberofelements=size(mh.elements,1);

	%Offset mh2.elements by number of vertices in the 3D structure: 
	mh2.elements=mh2.elements+mh.numberofvertices;
	mh2.segments(:,1:2)=mh2.segments(:,1:2)+mh.numberofvertices;
	mh2.segments(:,3)=mh2.segments(:,3)+mh.numberofelements;

	%The segments of md2 and the outer segments of md are identical. Go into  the elements of 
	%dmd2 and set them to their md equivalent: 
	for i=1:length(mh2.segments),
		node2=mh2.segments(i,1);
		%this node2 has an equivalent on the segments  of md: 
		for j=1:length(segments),
			node1=segments(j,1);
			if mh2.x(node2-mh.numberofvertices) == xsegs(j) &&  mh2.y(node2-mh.numberofvertices) == ysegs(j),
				%go into the mesh of md2, and replace by node1.
				pos=find(mh2.elements==node2); mh2.elements(pos)=node1;
				segs=mh2.segments(:,1:2); pos=find(segs==node2); segs(pos)=node1; mh2.segments(:,1:2)=segs;
				break;
			end
		end
	end

	%Do the merge: 
	mh.elements=[mh.elements;mh2.elements];
	mh.lat=[mh.lat;mh2.lat];
	mh.long=[mh.long;mh2.long];
	mh.segments=[mh.segments;mh2.segments];
	mh.numberofvertices=length(mh.lat);
	mh.numberofelements=size(mh.elements,1);
	
	%Remove orphans:
	lat=mh.lat; long=mh.long; 
	elements=mh.elements; segments=mh.segments;
	orphan=find(~ismember([1:length(lat)],sort(unique(elements(:)))));
	for i=1:length(orphan),
		%disp('WARNING: removing orphans');
		%get rid of the orphan node i
		%update lat and long
		lat=[lat(1:orphan(i)-(i-1)-1); lat(orphan(i)-(i-1)+1:end)];
		long=[long(1:orphan(i)-(i-1)-1); long(orphan(i)-(i-1)+1:end)];
		%update elements
		pos=find(elements>orphan(i)-(i-1));
		elements(pos)=elements(pos)-1;
		%update segments
		pos1=find(segments(:,1)>orphan(i)-(i-1));
		pos2=find(segments(:,2)>orphan(i)-(i-1));
		segments(pos1,1)=segments(pos1,1)-1;
		segments(pos2,2)=segments(pos2,2)-1;
	end
	
	mh.elements=elements;
	mh.lat=lat;
	mh.long=long;
	mh.segments=segments;
	mh.numberofelements=length(mh.elements);
	mh.numberofvertices=length(mh.lat);

	%reconstruct x,y and z:
	R=mh.r(1);
	mh.x = R .* cosd(mh.lat) .* cosd(mh.long);
	mh.y = R .* cosd(mh.lat) .* sind(mh.long);
	mh.z = R .* sind(mh.lat);
	mh.r=sqrt(mh.x.^2+mh.y.^2+mh.z.^2);

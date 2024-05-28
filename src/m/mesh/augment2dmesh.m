function [mh,internalsegments,internalelements]=augment2dmesh(mh,mhband,varargin)
%AUGMENT2DMESH - augment mh mesh with a band around it (provided by mhband)
%
%   Usage:
%      mh=augment2dmesh(mh,mhband);
%      [mh,internalsegments]=augment2dmesh(mh,mhband);
%
%   Example: 
%      md.mesh=augment2dmesh(md.mesh,md2.mesh);
%

%First process options
options=pairoptions(varargin{:});

%Offset the mesh band elements: 
mhband.elements=mhband.elements+mh.numberofvertices;
mhband.segments(:,1:2)=mhband.segments(:,1:2)+mh.numberofvertices;
mhband.segments(:,3)=mhband.segments(:,3)+mh.numberofelements;

bandsegments=mhband.segments; %keep  a copy around for when user requests this be returned. 

%The innner segments of mhband and the outer segments of mh are identical. Go into  the elements of 
%mhband and set them to their md1 equivalent: 
tol=getfieldvalue(options,'tolerance',1); %1 meter 
for i=1:length(mhband.segments),
	node2=mhband.segments(i,1);
	%this node2 has an equivalent on the segments  of mdh: 
	for j=1:length(mh.segments),
		node1=mh.segments(j,1);
		%if mhband.x(node2-mh.numberofvertices) == mh.x(node1) &&  mhband.y(node2-mh.numberofvertices) == mh.y(node1),
		if sqrt((mhband.x(node2-mh.numberofvertices) - mh.x(node1))^2 + (mhband.y(node2-mh.numberofvertices) - mh.y(node1))^2)<tol,
			%go into the mesh of mhband, and replace by node1.
			pos=find(mhband.elements==node2); mhband.elements(pos)=node1;
			segs=mhband.segments(:,1:2); pos=find(segs==node2); segs(pos)=node1; mhband.segments(:,1:2)=segs; bandsegments(:,1:2)=segs;
			%flag the segment for future erasal: 
			pos=find(mhband.segments(:,1)==node1 | mhband.segments(:,2)==node1); mhband.segments(pos,3)=NaN; 
			break;
		end
	end
end

%Do the merge: 
internalelements=mh.elements;
mh.elements=[mh.elements;mhband.elements];
mh.x=[mh.x;mhband.x];
mh.y=[mh.y;mhband.y];
if ~isnan(mh.lat), mh.lat=[mh.lat;mhband.lat]; end
if ~isnan(mh.long), mh.long=[mh.long;mhband.long];  end
mh.segments=[mhband.segments]; 

%segments that are internal mght have been requested: 
if strcmpi(getfieldvalue(options,'internalsegments','off'),'on'),
	pos=find(isnan(mh.segments(:,3)));
	internalsegments=bandsegments(pos,:);
end

%get rid of segments that are internal: 
pos=find(~isnan(mh.segments(:,3)));
mh.segments=mh.segments(pos,:); 

%Remove orphans:
x=mh.x; y=mh.y; lat=mh.lat; long=mh.long; 
elements=mh.elements; segments=mh.segments;
orphan=find(~ismember([1:length(x)],sort(unique(elements(:)))));
for i=1:length(orphan),
	%disp('WARNING: removing orphans');
	%get rid of the orphan node i
	%update x and y
	x=[x(1:orphan(i)-(i-1)-1); x(orphan(i)-(i-1)+1:end)];
	y=[y(1:orphan(i)-(i-1)-1); y(orphan(i)-(i-1)+1:end)];
	if ~isnan(lat), lat=[lat(1:orphan(i)-(i-1)-1); lat(orphan(i)-(i-1)+1:end)]; end
		if ~isnan(long), long=[long(1:orphan(i)-(i-1)-1); long(orphan(i)-(i-1)+1:end)]; end
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
mh.x=x;
mh.y=y;
mh.lat=lat;
mh.long=long;
mh.segments=segments;
mh.numberofelements=length(mh.elements);
mh.numberofvertices=length(mh.x);

%vertex on boundaries: 
mh.vertexonboundary=zeros(mh.numberofvertices,1);
mh.vertexonboundary(mh.segments(:,1:2))=1;

%connectivities: 
mh.vertexconnectivity=NodeConnectivity(mh.elements,mh.numberofvertices);
mh.elementconnectivity=ElementConnectivity(mh.elements,mh.vertexconnectivity);

%return: 
if nargout==1,
	varargout{1}=mh;
elseif nargout==2,
	varargout{1}=mh;
	varargout{2}=internalsegments;
elseif nargout==3,
	varargout{1}=mh;
	varargout{2}=internalsegments;
	varargout{3}=internalelements;
else 
	error('unknow number of output arguments')
end

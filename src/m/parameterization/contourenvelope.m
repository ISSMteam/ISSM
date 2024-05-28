function segments=contourenvelope(mh,varargin)
%CONTOURENVELOPE - build a set of segments enveloping a contour .exp
%
%   Usage:
%      segments=contourenvelope(mh,varargin)
%
%   Example:
%      segments=contourenvelope(mh,'Stream.exp');
%      segments=contourenvelope(mh);

%some checks
if nargin>2,
	help contourenvelope
	error('contourenvelope error message: bad usage');
end
if nargin==2,
	flags=varargin{1};

	if ischar(flags),
		file=flags;
		if ~exist(file),
			error(['contourenvelope error message: file ' file ' not found']);
		end
		isfile=1;
	elseif isnumeric(flags),
		%do nothing for now
		isfile=0;
	else
		error('contourenvelope error message: second argument should be a file or an elements flag');
	end
end

%Now, build the connectivity tables for this mesh.
%Computing connectivity
if isnan(mh.vertexconnectivity),
	mh.vertexconnectivity=NodeConnectivity(mh.elements,mh.numberofvertices);
end
if isnan(mh.elementconnectivity),
	mh.elementconnectivity=ElementConnectivity(mh.elements,mh.vertexconnectivity);
end

%get nodes inside profile
mesh.elementconnectivity=mh.elementconnectivity;
if dimension(mh)==2,
	mesh.elements=mh.elements;
	mesh.x=mh.x;
	mesh.y=mh.y;
	mesh.numberofvertices=mh.numberofvertices;
	mesh.numberofelements=mh.numberofelements;
else
	mesh.elements=mh.elements2d;
	mesh.x=mh.x2d;
	mesh.y=mh.y2d;
	mesh.numberofvertices=mh.numberofvertices2d;
	mesh.numberofelements=mh.numberofelements2d;
end

if nargin==2,

	if isfile,
		%get flag list of elements and nodes inside the contour
		nodein=ContourToMesh(mesh.elements,mesh.x,mesh.y,file,'node',1);
		elemin=(sum(nodein(mesh.elements),2)==size(mesh.elements,2));
		%modify element connectivity
		elemout=find(~elemin);
		mesh.elementconnectivity(elemout,:)=0;
		mesh.elementconnectivity(find(ismember(mesh.elementconnectivity,elemout)))=0;
	else
		%get flag list of elements and nodes inside the contour
		nodein=zeros(mesh.numberofvertices,1);
		elemin=zeros(mesh.numberofelements,1);

		pos=find(flags);
		elemin(pos)=1;
		nodein(mesh.elements(pos,:))=1;

		%modify element connectivity
		elemout=find(~elemin);
		mesh.elementconnectivity(elemout,:)=0;
		mesh.elementconnectivity(find(ismember(mesh.elementconnectivity,elemout)))=0;
	end
end

%Find element on boundary
%First: find elements on the boundary of the domain
flag=mesh.elementconnectivity;
if nargin==2,
	flag(find(flag))=elemin(flag(find(flag)));
end
elementonboundary=double(prod(flag,2)==0 & sum(flag,2)>0);

%Find segments on boundary
pos=find(elementonboundary);
num_segments=length(pos);
segments=zeros(num_segments*3,3);
count=1;

for i=1:num_segments,
	el1=pos(i);
	els2=mesh.elementconnectivity(el1,find(mesh.elementconnectivity(el1,:)));
	if length(els2)>1,
		flag=intersect(intersect(mesh.elements(els2(1),:),mesh.elements(els2(2),:)),mesh.elements(el1,:));
		nods1=mesh.elements(el1,:);
		nods1(find(nods1==flag))=[];
		segments(count,:)=[nods1 el1];

		ord1=find(nods1(1)==mesh.elements(el1,:));
		ord2=find(nods1(2)==mesh.elements(el1,:));

		%swap segment nodes if necessary
		if ( (ord1==1 & ord2==2) | (ord1==2 & ord2==3) | (ord1==3 & ord2==1) ),
			temp=segments(count,1);
			segments(count,1)=segments(count,2);
			segments(count,2)=temp;
		end
		segments(count,1:2)=fliplr(segments(count,1:2));
		count=count+1;
	else
		nods1=mesh.elements(el1,:);
		flag=setdiff(nods1,mesh.elements(els2,:));
		for j=1:3,
			nods=nods1; nods(j)=[];
			if any(ismember(flag,nods)),
				segments(count,:)=[nods el1];
				ord1=find(nods(1)==mesh.elements(el1,:));
				ord2=find(nods(2)==mesh.elements(el1,:));
				if ( (ord1==1 & ord2==2) | (ord1==2 & ord2==3) | (ord1==3 & ord2==1) ),
					temp=segments(count,1);
					segments(count,1)=segments(count,2);
					segments(count,2)=temp;
				end
				segments(count,1:2)=fliplr(segments(count,1:2));
				count=count+1;
			end
		end
	end
end
segments=segments(1:count-1,:);

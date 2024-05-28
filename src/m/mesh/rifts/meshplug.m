function [elements,x,y,z,numberofelements,numberofnodes,elconv,nodeconv,elconv2,nodeconv2]=meshplug(elements,x,y,z,elements2,x2,y2,z2,extractednodes,extractedelements,domain)
%MESHPLUG - embed mesh into another one
%     See also meshaddrifts

%initialize elconv,nodeconv conversion tables from md mesh to new md mesh
elconv=1:size(elements,1); elconv=elconv';
nodeconv=1:size(x,1); nodeconv=nodeconv';

%take away old elements in area of interest: 
elements(extractedelements,:)=[];
element_offset=size(elements,1);

%update elconv after having extracted the area of interest elements
temp_elconv=elconv; temp_elconv(extractedelements)=[];
temp_elconvnum=1:length(temp_elconv);
elconv(temp_elconv)=temp_elconvnum;
elconv(extractedelements)=NaN;

%initialize elconv2 and nodeconv2, conversion tables from md2 mesh to new md mesh
elconv2=1:size(elements2,1);elconv2=elconv2'+element_offset;
nodeconv2=(size(x,1)+1):(size(x,1)+size(x2,1)); nodeconv2=nodeconv2';

extractednodes_minusborder=extractednodes;
extractednodes_minusborder(domain)=[];

x(extractednodes_minusborder)=NaN;
y(extractednodes_minusborder)=NaN;

%now, plug md2 mesh: 

%first, offset all ids of md2 mesh
elements2=elements2+length(x);

%NaN border nodes in second mesh
x2(1:length(domain))=NaN;
y2(1:length(domain))=NaN;

%redirect border nodes in elements2  to elements
for i=1:length(domain),
	pos=find(elements2==(i+length(x)));
	elements2(pos)=extractednodes(domain(i));
end

%same deal for nodeconv2:
for i=1:length(domain),
	nodeconv2(i)=extractednodes(domain(i));
end

%plug elements
elements=[elements;elements2];

%now, increase number of nodes
x=[x; x2];
y=[y; y2];
z=[z; z2];

%now, get rid of NaN in x:
while  ~isempty(find(isnan(x))),

	pos=find(isnan(x));
	node=pos(1);

	%collapse node
	x(node)=[];
	y(node)=[];
	z(node)=[];

	%renumber all nodes > node in elements
	pos=find(elements>node);
	elements(pos)=elements(pos)-1;

	%same deal for nodeconv2: 
	pos=find(nodeconv2>node);
	nodeconv2(pos)=nodeconv2(pos)-1;

end

numberofnodes=length(x);
numberofelements=length(elements);

%finish nodeconv: 
temp_nodeconv=nodeconv;  temp_nodeconv(extractednodes_minusborder)=[];
temp_nodeconvnum=1:length(temp_nodeconv);
nodeconv(temp_nodeconv)=temp_nodeconvnum;
nodeconv(extractednodes_minusborder)=NaN;

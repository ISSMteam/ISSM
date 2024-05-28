function md=argusmesh(md,infile)
%ARGUSMESH - load an Argus mesh onto a model
%
%   Convert an Argus mesh contained in a file into
%   fields needed for the mesh in a model md.
%
%   Usage:
%      md=argusmesh(md,infile)
%
%   Example:
%     md=argusmesh(md,'Domain.exp')

%some argument check: 
if nargin~=2 | nargout~=1,
	help argustomodel;
	error('argustomodel error message: bad usage');
end

%determine root of infile: strip extension
[a,root,b,c]=fileparts(infile);

%inform user we start the script: 
disp(['   Translating argus file ''' infile ''' into matlab model object']);

%open infile: 
fileid=fopen(infile,'r');
if fileid==-1,
	error(['Could not open file ' infile  ' for reading']);
end

%Read first line of the argus mesh: node and element parameters
[buffer,bytecount]=fscanf(fileid,'%i %i %i %i',[1 4]);
if bytecount~=4, 
	error(['Problem reading ' infile ' file at line #1']);
end
nel=buffer(1);
nods=buffer(2);
num_element_parameters=buffer(3);
num_node_parameters=buffer(4);
disp(['      argus model '''   root ''' contains ' num2str(nel) ' elements and ' num2str(nods) ' nodes.']);

%initialize elements and nodes
elements=zeros(nel,3);
element_parameters=zeros(nel,num_element_parameters);
x=zeros(nods,1);
y=zeros(nods,1);
z=zeros(nods,1);
node_parameters=zeros(nods,num_node_parameters);

%read nodes:
format_string='%s %i %f %f ';
for n=1:num_node_parameters,
	format_string=[format_string ' %i '];
end

for n=1:nods,
	[buffer,bytecount]=fscanf(fileid,format_string,[1,num_node_parameters+4]);
	x(n)=buffer(3);
	y(n)=buffer(4);
	node_parameters(n,:)=buffer(5:length(buffer));
end

%read elements: 
format_string='%s %i %i %i %i';
for n=1:num_element_parameters,
	format_string=[format_string ' %i '];
end
for n=1:nel,
	[buffer,bytecount]=fscanf(fileid,format_string,[1,num_element_parameters+5]);
	elements(n,:)=buffer(3:5);
	element_parameters(n,:)=buffer(6:length(buffer));
end

%Create a name and a note for this model: 
notes=['Model created by Argus from input file: ' infile ' and parameter file: ' root '.par on: ' date];
name=root;

%Finally, use model constructor to build a complete model: 
md.mesh=mesh2d();
md.mesh.elements=elements;
md.mesh.x=x;
md.mesh.y=y;
md.mesh.numberofvertices=size(md.mesh.x,1);
md.mesh.numberofelements=size(md.mesh.elements,1);
md=addnote(md,notes);

%Add segments and nodes on boundary
md.mesh.segments=findsegments(md);
md.mesh.vertexonboundary=zeros(md.mesh.numberofvertices,1);
md.mesh.vertexonboundary(md.mesh.segments(:,1))=1;
md.mesh.vertexonboundary(md.mesh.segments(:,2))=1;

function [index,X,Y,Z,S,data_interp]=SectionValues(md,data,infile,resolution)
%SECTIONVALUES - compute the value of a field on a section
%
%   This routine gets the value of a given field of the model on points
%   given by filname (Argus type file)
%
%   Usage:
%      [elements,x,y,z,s,data]=SectionValues(md,data,filename,resolution)
%      [elements,x,y,z,s,data]=SectionValues(md,data,profile_structure,resolution)

%check what we have for profile as input
if ischar(infile),
	%read infile:
	profile=expread(infile);
	nods=profile.nods;
	x=profile.x;
	y=profile.y;
else
	%read infile:
	nods=infile.nods;
	x=infile.x;
	y=infile.y;
end

%get the specified resolution
if isnumeric(resolution(1))
	res_h=resolution(1);
else
	error('SectionValues error message: wrong resolution type. Resolution must be an array [horizontal_resolution vertical_resolution]')
end
if dimension(md.mesh)==3
	if (length(resolution)==2 & isnumeric(resolution(2)))
		res_v=resolution(2);
	else
		error('SectionValues error message: wrong resolution type. Resolution must be an array [horizontal_resolution vertical_resolution]')
	end
end

%initialization
X=[]; %X-coordinate
Y=[]; %Y-coordinate
S=0;  %curvilinear coordinate

for i=1:nods-1

	x_start=x(i);
	x_end=x(i+1);
	y_start=y(i);
	y_end=y(i+1);
	s_start=S(end);

	length_segment=sqrt((x_end-x_start)^2+(y_end-y_start)^2);
	portion=ceil(length_segment/res_h);

	x_segment=zeros(portion,1);
	y_segment=zeros(portion,1);
	s_segment=zeros(portion,1);

	for j=1:portion
		x_segment(j)=x_start+(j-1)*(x_end-x_start)/portion;
		y_segment(j)=y_start+(j-1)*(y_end-y_start)/portion;
		s_segment(j)=s_start+j*length_segment/portion;
	end

	%plug into X and Y
	X=[X;x_segment];
	Y=[Y;y_segment];
	S=[S;s_segment];
end
X(end+1)=x(nods);
Y(end+1)=y(nods);

%Number of nodes:
numberofnodes=size(X,1);

%Compute Z
Z=zeros(numberofnodes,1);

%New mesh and Data interpolation
if (dimension(md.mesh)==2)

	%Interpolation of data on specified points
	data_interp=InterpFromMesh2d(md.mesh.elements,md.mesh.x,md.mesh.y,data,X,Y);
	%data_interp=InterpFromMeshToMesh2d(md.mesh.elements,md.mesh.x,md.mesh.y,data,X,Y);
	%data_interp=griddata(md.mesh.x,md.mesh.y,data,X,Y);

	%Compute index
	index=[1:1:(numberofnodes-1);2:1:numberofnodes]';

else

	%vertically extrude mesh

	%Get base and surface for each 2d point, offset to make sure that it is inside the glacier system
	offset=10^-3;
	base=InterpFromMeshToMesh2d(md.mesh.elements2d,md.mesh.x2d,md.mesh.y2d,project2d(md,md.geometry.base,1),X,Y)+offset;
	surface=InterpFromMeshToMesh2d(md.mesh.elements2d,md.mesh.x2d,md.mesh.y2d,project2d(md,md.geometry.surface,1),X,Y)-offset;

	%Some useful parameters
	layers=ceil(mean(md.geometry.thickness)/res_v);
	nodesperlayer=numberofnodes;
	nodestot=nodesperlayer*layers;
	elementsperlayer=nodesperlayer-1;
	elementstot=(nodesperlayer-1)*(layers-1);

	%initialization
	X3=zeros(nodesperlayer*layers,1); Y3=zeros(nodesperlayer*layers,1); Z3=zeros(nodesperlayer*layers,1); S3=zeros(nodesperlayer*layers,1); index3=zeros(elementstot,4);

	%Get new coordinates in 3d
	for i=1:layers
		X3(i:layers:end)=X;
		Y3(i:layers:end)=Y;
		Z3(i:layers:end)=base+(i-1)*(surface-base)/(layers-1);
		S3(i:layers:end)=S;

		if i<layers %Build index3 with quads
			index3((i-1)*elementsperlayer+1:i*elementsperlayer,:)=[i:layers:nodestot-layers; i+1:layers:nodestot-layers; i+layers+1:layers:nodestot; i+layers:layers:nodestot]';
		end
	end

	%Interpolation of data on specified points
	data_interp=InterpFromMeshToMesh3d(md.mesh.elements,md.mesh.x,md.mesh.y,md.mesh.z,data,X3,Y3,Z3,NaN);

	%build outputs
	X=X3; Y=Y3; Z=Z3;  S=S3; index=index3;
end

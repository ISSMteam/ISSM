function data_out=InterpFromFile(x,y,filename,default_value)
%INTERPFROMFILE - load data and interpolate on the given nodes
%
%   load a matlab file (extension .mat) which holds 3 or 4 variables
%   and interpolate the data on the mesh and plug it onto the model.
%
%   o 3 variables
%     - a vector x (if the name of the variable do not begin with "x", an error can appear)
%     - a vector y (if the name of the variable do not begin with "y", an error can appear)
%     - a vector or matrix data (if the name of the variable do not begin with the field name, an error can appear)
%   o 4 variables
%     - a vector x (if the name of the variable do not begin with "x", an error can appear)
%     - a vector y (if the name of the variable do not begin with "y", an error can appear)
%     - a matrix with 3 columns (if the name of the variable do not begin with "index" or "elements", an error can appear)
%     - a vector data (if the name of the variable do not begin with the field name, an error can appear)
%
%   Usage:
%      data=InterpFromFile(x,y,filename,default_value);
%
%   Example:
%      md.geometry.surface=InterpFromFile(md.mesh.x,md.mesh.y,'surfacefile.mat',0);
%
%   See also: PLUGVELOCITIES, INTERPFROMGRID, INTERPFROMMESH2D, INTERPFROMMESH3D

%some checks
if nargin~=4 | nargout~=1
	help InterpFromFile
	error('plugdata error message: bad usage');
end
if ~exist(filename)
	error(['plugdata error message: file ' filename  ' does not exist']);
end
if length(x)~=length(y),
	error('plugdata error message: x and y should have the same length');
end

%load file
Names=FieldFindVarNames(filename);
Data=load(filename);
disp('WARNING: function deprecated, replace InterpFromFile by the following command:');
if strcmpi(Names.interp,'node'),
	disp(['   data=InterpFromGridToMesh(' Names.xname ',' Names.yname ',' Names.dataname ',x,y,default_value);']);
	data_out=InterpFromGridToMesh(Data.(Names.xname),Data.(Names.yname),Data.(Names.dataname),x,y,default_value);
else
	disp(['   data=InterpFromMeshToMesh2d(' Names.indexname ',' Names.xname ',' Names.yname ',' Names.dataname ',x,y);']);
	data_out=InterpFromMeshToMesh2d(Data.(Names.indexname),Data.(Names.xname),Data.(Names.yname),Data.(Names.dataname),x,y);
end

end
function Names=FieldFindVarNames(filename)
%FIELDFINDVARNAMES - find names of variables in a data set file
%
%   This routines looks at the variables contained in a file and finds out
%   the names of the variables that are needed for an interpolation (x,y,data)
%   or (index,x,y,data)
%
%   Usage:
%      Names=FieldFindVarNames(filename)
%
%   Example:
%      Names=FieldFindVarNames('thickness.mat')
%
%   See also: INTERPFROMFILE, GRIDDATA

%some checks
if nargin~=1 | nargout~=1
	help FieldFindVarNames
	error('FieldFindVarNames error message: bad usage');
end
if ~exist(filename)
	error(['FieldFindVarNames error message: file ' filename  ' does not exist']);
end

%Get variables
A=whos('-file',filename);

%find x,y,vx and vy
xenum=NaN; yenum=NaN; dataenum=NaN; indexenum=NaN;
if length(A)==3,
	isnode=1;
	for i=1:3
		if strcmpi(A(i).name(1),'x');
			xenum=i;
		elseif strcmpi(A(i).name(1),'y');
			yenum=i;
		elseif (strncmpi(A(i).name,filename,3) | strncmpi(A(i).name,'data',4)),
			dataenum=i;
		else
			%nothing
		end
	end
elseif length(A)==4,
	isnode=0;
	for i=1:4
		if strcmpi(A(i).name(1),'x');
			xenum=i;
		elseif strcmpi(A(i).name(1),'y');
			yenum=i;
		elseif (strncmpi(A(i).name,'index',5) | strncmpi(A(i).name,'elements',7));
			indexenum=i;
		elseif (strncmpi(A(i).name,filename,3) | strncmpi(A(i).name,'data',4)),
			dataenum=i;
		else
			%nothing
		end
	end
else
	error(['FieldFindVarNames error message: file ' filename  ' not supported yet (it should hold 3 variables x,y and data (for nodes) OR 4 variables  x,y,index and data (for mesh))']);
end

%2: if only one item is missing, find it by elimination
if ~isnode,
	pos=find(isnan([xenum yenum indexenum dataenum]));
	if length(pos)==1,
		list=[xenum yenum indexenum dataenum]; list(pos)=[];
		if pos==1,
			xenum=setdiff(1:4,list);
		elseif pos==2,
			yenum=setdiff(1:4,list);
		elseif pos==3,
			indexenum=setdiff(1:4,list);
		elseif pos==4,
			dataenum=setdiff(1:4,list);
		end
	end
else
	pos=find(isnan([xenum yenum dataenum]));
	if length(pos)==1,
		list=[xenum yenum indexenum dataenum]; list(pos)=[];
		if pos==1,
			xenum=setdiff(1:3,list);
		elseif pos==2,
			yenum=setdiff(1:3,list);
		elseif pos==3,
			dataenum=setdiff(1:3,list);
		end
	end
end

%assum that we have found at least xenum and yenum
if ( isnan(xenum) | isnan(yenum))
	error(['FieldFindVarNames error message: file ' filename  ' not supported yet (the coordinates vectors should be named x and y)']);
end

%find index
if (~isnode & isnan(indexenum)),
	for i=1:4
		lengthi=min(A(i).size);
		if (lengthi==3),
			indexenum=i;
		end
	end
	if isnan(indexenum),
		error(['FieldFindVarNames error message: file ' filename  ' not supported yet (index not found)']);
	end
end

%4: last chance
if ~isnode,
	pos=find(isnan([xenum yenum indexenum dataenum]));
	if length(pos)==1,
		list=[xenum yenum indexenum dataenum]; list(pos)=[];
		if pos==1,
			xenum=setdiff(1:4,list);
		elseif pos==2,
			yenum=setdiff(1:4,list);
		elseif pos==3,
			indexenum=setdiff(1:4,list);
		elseif pos==4,
			dataenum=setdiff(1:4,list);
		end
	end
else
	pos=find(isnan([xenum yenum dataenum]));
	if length(pos)==1,
		list=[xenum yenum indexenum dataenum]; list(pos)=[];
		if pos==1,
			xenum=setdiff(1:3,list);
		elseif pos==2,
			yenum=setdiff(1:3,list);
		elseif pos==3,
			dataenum=setdiff(1:3,list);
		end
	end
end

%last check
if isnan(dataenum)
	error(['FieldFindVarNames error message: file ' filename  ' not supported yet (data not found)']);
end

%create output
Names=struct();
Names.xname=A(xenum).name;
Names.yname=A(yenum).name;
Names.dataname=A(dataenum).name;
if ~isnode,
	Names.indexname=A(indexenum).name; 
	Names.interp='mesh';
else
	Names.interp='node';
end
end

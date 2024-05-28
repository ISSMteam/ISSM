function Struct=expread(filename,varargin)
%EXPREAD - read a file exp and build a Structure
%
%   This routine reads a file .exp and build a Structure containing the 
%   fields x and y corresponding to the coordinates, one for the filename of
%   the exp file, for the density, for the nodes, and a field closed to 
%   indicate if the domain is closed. 
%   The first argument is the .exp file to be read and the second one (optional) 
%   indicate if the last point shall be read (1 to read it, 0 not to).
%
%   Usage:
%      Struct=expread(filename)
%
%   Example:
%      Struct=expread('domainoutline.exp')
%      Struct=expread('domainoutline.exp')
%
%   See also EXPDOC, EXPWRITEASVERTICES

%recover options
options=pairoptions(varargin{:});

%some checks
if ~exist(filename),
	error(['expread error message: file ' filename ' not found!']);
end

%initialize number of profile
count=0;
Struct = struct();

%open file
fid=fopen(filename,'r');

%loop over the number of profiles
while (~feof(fid)),

	%update number of profiles
	count=count+1;

	%Get file name
	A=fscanf(fid,'%s %s',2);
	if ~strncmp(A,'##Name:',7), break; end
	if length(A)>7, 
		Struct(count).name=A(8:end);
	else
		Struct(count).name='';
	end

	%Get Icon
	A=fscanf(fid,'%s %s',2);
	if ~strncmp(A,'##Icon:',7), break; end

	%Get Info
	A=fscanf(fid,'%s %s %s %s',4);
	if ~strncmp(A,'#Points',7), break; end

	%Get number of nodes and density
	A=fscanf(fid,'%f %f',[1 2]);
	Struct(count).nods=A(1);
	Struct(count).density=A(2);

	%Get Info
	A=fscanf(fid,'%s %s %s %s',5);
	if ~strncmp(A,'#XposYpos',9), break; end

	%Get Coordinates
	A=fscanf(fid,'%f %f',[2 Struct(count).nods]);
	Struct(count).x=A(1,:)';
	Struct(count).y=A(2,:)';
	if any(isnan(A))
		warning('NaNs found in coordinates, note that some tools like exptool will not work properly with NaNs');
	end

	if(Struct(count).nods~=length(Struct(count).x))error(['Profile ' num2str(count) ' reports incorrect length']); end;

	%Check if closed
	if (Struct(count).nods > 1) && ...
	   (Struct(count).x(end) == Struct(count).x(1)) && ...
	   (Struct(count).y(end) == Struct(count).y(1))
		Struct(count).closed=true;
	else
		Struct(count).closed=false;
	end
end


%close file
fclose(fid);
	
invert=getfieldvalue(options,'invert',0);
if invert,
	for i=1:length(Struct),
		Struct(i).x=flipud(Struct(i).x);
		Struct(i).y=flipud(Struct(i).y);
	end
end

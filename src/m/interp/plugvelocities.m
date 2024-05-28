function md=plugvelocities(md,filename,default_value)
%PLUGVELOCITIES - load velocities on a model
%
%   load a matlab file (extension .mat) which holds 4 variables
%   x,y,vx,vy to be plugged onto the model (or similar names)
%   x and y must be vectors, vx, vy matrices
%
%   Usage:
%      md=plugvelocities(md,filename,default_value)
%
%   Example:
%      md=plugvelocities(md,'velocityfile.mat',0);
%
%   See also: INTERPFROMFILE, GRIDDATA

disp('WARNING: deprecated functions (plugvelocities)');
%some checks
if nargin~=3 | nargout~=1
	help plugvelocities
	error('plugvelocities error message: bad usage');
end
if ~exist(filename)
	error(['plugvelocities error message: file ' filename  ' does not exist']);
end

%load velocities 
Names=VelFindVarNames(filename);
Vel=load(filename);

%Interpolation
if strcmpi(Names.interp,'node'),
	md.inversion.vx_obs=InterpFromGridToMesh(Vel.(Names.xname),Vel.(Names.yname),Vel.(Names.vxname),md.mesh.x,md.mesh.y,default_value);
	md.inversion.vy_obs=InterpFromGridToMesh(Vel.(Names.xname),Vel.(Names.yname),Vel.(Names.vyname),md.mesh.x,md.mesh.y,default_value);
else
	md.inversion.vx_obs=InterpFromMeshToMesh2d(Vel.(Names.indexname),Vel.(Names.xname),Vel.(Names.yname),Vel.(Names.vxname),md.mesh.x,md.mesh.y,default_value);
	md.inversion.vy_obs=InterpFromMeshToMesh2d(Vel.(Names.indexname),Vel.(Names.xname),Vel.(Names.yname),Vel.(Names.vyname),md.mesh.x,md.mesh.y,default_value);
end

md.inversion.vel_obs=sqrt(md.inversion.vx_obs.^2+md.inversion.vy_obs.^2);
md.initialization.vx=md.inversion.vx_obs;
md.initialization.vy=md.inversion.vy_obs;
md.initialization.vel=md.inversion.vel_obs;
end

function Names=VelFindVarNames(filename)
%VELFINDVARNAMES - find names of variables in a velocity data set file
%
%   This routines looks at the variables contained in a file and finds out
%   the names of the variables that are needed for an interpolation (x,y,vx,vy)
%   or (index,x,y,vx,vy)
%
%   Usage:
%      Names=VelFindVarNames(filename)
%
%   Example:
%      Names=VelFindVarNames('velocities.mat')
%
%   See also: INTERPFROMFILE, GRIDDATA

%some checks
if nargin~=1 | nargout~=1
	help VelFindVarNames
	error('VelFindVarNames error message: bad usage');
end
if ~exist(filename)
	error(['VelFindVarNames error message: file ' filename  ' does not exist']);
end

%Get variables
A=whos('-file',filename);

%find x,y,vx and vy
xenum=NaN; yenum=NaN; vxenum=NaN; vyenum=NaN; indexenum=NaN;
if length(A)==4,
	isnode=1;
	for i=1:4
		if strcmpi(A(i).name(1),'x');
			xenum=i;
		elseif strcmpi(A(i).name(1),'y');
			yenum=i;
		else
			if (strcmpi(A(i).name(end),'x') | strncmpi(A(i).name,'vx',2));
				vxenum=i;
			elseif (strcmpi(A(i).name(end),'y') | strncmpi(A(i).name,'vy',2));
				vyenum=i;
			end
		end
	end
elseif length(A)==5,
	isnode=0;
	for i=1:5
		if strcmpi(A(i).name(1),'x');
			xenum=i;
		elseif strcmpi(A(i).name(1),'y');
			yenum=i;
		elseif (strcmpi(A(i).name(1),'index') | strcmpi(A(i).name(1),'elements'));
			indexenum=i;
		else
			if (strcmpi(A(i).name(end),'x') | strncmpi(A(i).name,'vx',2));
				vxenum=i;
			elseif (strcmpi(A(i).name(end),'y') | strncmpi(A(i).name,'vy',2));
				vyenum=i;
			end
		end
	end
else
	error(['VelFindVarNames error message: file ' filename  ' not supported yet (it should hold 4 variables x,y,vx and vy (for nodes) OR 5 variables  x,y,index,vx and vy (for mesh))']);
end

%assum that we have found at least vxenum and vyenum
if ( isnan(vxenum) | isnan(vyenum))
	error(['VelFindVarNames error message: file ' filename  ' not supported yet (the velocities should be named vx and vy)']);
end

%find index
if (~isnode & isnan(indexenum)),
	for i=1:5
		lengthi=min(A(i).size);
		if (lengthi==3),
			indexenum=i;
		end
	end
	if isnan(indexenum),
		error(['VelFindVarNames error message: file ' filename  ' not supported yet (index not found)']);
	end
end

%find x y
if (isnan(xenum) | isnan(yenum))

	%check the size
	if A(vxenum).size(1)==A(vxenum).size(2),
		error(['VelFindVarNames error message: file ' filename  ' not supported (velocities is a square matrix, save x and y with another name)']);
	end
	if ~(A(vxenum).size(1)==A(vyenum).size(1) & A(vxenum).size(2)==A(vyenum).size(2)),
		error(['VelFindVarNames error message: file ' filename  ' not supported (vx and vy matrices do not have the same size)']);
	end

	%find xenum and yenum
	for i=1:4
		lengthi=max(A(i).size);
		if ((i~=vxenum) & (lengthi==A(vxenum).size(1) | lengthi==A(vxenum).size(1)+1)),
			yenum=i;
		elseif ((i~=vxenum) & (lengthi==A(vxenum).size(2) | lengthi==A(vxenum).size(2)+1)),
			xenum=i;
		end
	end

	%last check
	if (isnan(xenum) | isnan(yenum))
		error(['plugdata error message: file ' filename  ' not supported yet']);
	end
end

%create output
Names=struct();
Names.xname=A(xenum).name;
Names.yname=A(yenum).name;
Names.vxname=A(vxenum).name;
Names.vyname=A(vyenum).name;
if ~isnode,
	Names.indexname=A(indexenum).name; 
	Names.interp='mesh';
else
	Names.interp='node';
end
end

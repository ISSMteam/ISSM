function [vxout vyout] = interpJoughin(X,Y,Date),
	%Available dates:
	% 2000 2005 2006 2007 2008

switch oshostname(),
	case {'murdo','thwaites','astrid'}
		if nargin==3,
			rootname = ['/u/astrid-r1b/morlighe/issmjpl/proj-morlighem/DatasetGreenland/Data/Vel/Joughin/' num2str(Date) '/'];
		else
			error('not supported');
		end
	case {'ronne'}
		error('not supported');
	otherwise
		error('machine not supported yet');
end
verbose = 1;

if ~exist(rootname,'dir'),
	error(['file ' rootname ' not found']);
end

rootname = [rootname 'greenland_vel_mosaic500_' num2str(Date) '_' num2str(Date+1)];

if verbose, disp('   -- Joughin: loading vx'); end
[data,R] = geotiffread([rootname '_vx.tif']);
pos=find(data<-10^9); data(pos)=NaN;
data=double(flipud(data));
xdata=R.XLimWorld(1):R.DeltaX:R.XLimWorld(2); xdata=xdata(:);
xdata =(xdata(1:end-1)+xdata(2:end))/2;
ydata=R.YLimWorld(2):R.DeltaY:R.YLimWorld(1); ydata=flipud(ydata(:));
ydata =(ydata(1:end-1)+ydata(2:end))/2;
if verbose, disp('   -- Joughin: interpolating vx'); end
vxout = InterpFromGrid(xdata,ydata,data,X,Y);
vxout = reshape(vxout,size(X,1),size(X,2));

if verbose, disp('   -- Joughin: loading vy'); end
[data,R] = geotiffread([rootname '_vy.tif']);
pos=find(data<-10^9); data(pos)=NaN;
data=double(flipud(data));
xdata=R.XLimWorld(1):R.DeltaX:R.XLimWorld(2); xdata=xdata(:);
xdata =(xdata(1:end-1)+xdata(2:end))/2;
ydata=R.YLimWorld(2):R.DeltaY:R.YLimWorld(1); ydata=flipud(ydata(:));
ydata =(ydata(1:end-1)+ydata(2:end))/2;
if verbose, disp('   -- Joughin: interpolating vy'); end
vyout = InterpFromGrid(xdata,ydata,data,X,Y);
vyout = reshape(vyout,size(X,1),size(X,2));
return

% Get geodat info
if verbose, disp('   -- Joughin: loading geodat info'); end
xd=readgeodat(strcat(rootname,'.vx.geodat'));
xmin=xd(3,1)*1000.+xd(2,1)/2;
xmax=xd(3,1)*1000.+(xd(2,1)-1)*xd(1,1)+xd(2,1)/2;
ymin=xd(3,2)*1000.+xd(2,2)/2;
ymax=xd(3,2)*1000.+(xd(2,2)-1)*xd(1,2)+xd(2,2)/2;
%xmin=xd(3,1)*1000.;
%xmax=xd(3,1)*1000.+(xd(2,1)-1)*xd(1,1);
%ymin=xd(3,2)*1000.;
%ymax=xd(3,2)*1000.+(xd(2,2)-1)*xd(1,2);
xdata=linspace(xmin,xmax,xd(1,1));
ydata=linspace(ymin,ymax,xd(1,2));

% Vx component
if verbose, disp('   -- Joughin: loading vx'); end
fid = fopen(strcat(rootname,'.vx'),'r','ieee-be');
[data,count]=fread(fid,[xd(1,1) xd(1,2)],'float32');
fclose(fid);

if verbose, disp('   -- Joughin: interpolating vx'); end
vxout = InterpFromGrid(xdata,ydata,data',X,Y);
vxout = reshape(vxout,size(X,1),size(X,2));

% Vy component
fid = fopen(strcat(rootname,'.vy'),'r','ieee-be');
[data,count]=fread(fid,[xd(1,1) xd(1,2)],'float32');
fclose(fid);
vyout = InterpFromGrid(xdata,ydata,data',X,Y);
vyout = reshape(vyout,size(X,1),size(X,2));

end

function xgeo=readgeodat(filein)
% Read a geodat file
fid = fopen(filein,'r');
xgeo=zeros(3,2);
i=1;
while ~feof(fid),
	line=fgets(fid);
	[A,count]=sscanf(line,'%f %f',[1 2]);
	if(count == 2) 
		xgeo(i,:)=A;
		i=i+1;
	end
end
fclose(fid);
end

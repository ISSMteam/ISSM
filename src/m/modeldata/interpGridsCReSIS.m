function output = interpGridsCReSIS(X,Y,filename),

%Convert to lat/lon
disp('   -- Griggs2013: converting coordinates');
[LAT,  LON  ] = xy2ll(double(X(:)),double(Y(:)),+1,45,70);

disp(['   -- GridsCReSIS: loading data']);
if ~exist(filename)
	error([filename ' does not exist']);
end
fid   = fopen(filename);
for i=1:6,
	thisline = fgetl(fid);
	dummy    = regexp(thisline,'(\S+)','match');
	if strcmp(dummy{1},'ncols'),       ncols=str2num(dummy{2}); end
	if strcmp(dummy{1},'nrows'),       nrows=str2num(dummy{2}); end
	if strcmp(dummy{1},'xllcorner'),    xllcorner=str2num(dummy{2}); end
	if strcmp(dummy{1},'yllcorner'),    yllcorner=str2num(dummy{2}); end
	if strcmp(dummy{1},'cellsize'),     cellsize=str2num(dummy{2}); end
	if strcmp(dummy{1},'NODATA_value'), nodata=str2num(dummy{2}); end
end
data  = fscanf(fid,'%g %g %g %g %g',[ncols nrows])';
fclose(fid);

xdata=linspace(xllcorner+cellsize/2,xllcorner+cellsize/2+(ncols-1)*cellsize,ncols);
ydata=linspace(yllcorner+cellsize/2,yllcorner+cellsize/2+(nrows-1)*cellsize,nrows);

disp(['   -- GridsCReSIS: interpolating ']);
output = InterpFromGrid(xdata,ydata,data,LAT,LON);
output = reshape(output,size(X,1),size(X,2));

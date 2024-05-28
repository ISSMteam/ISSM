function output = interpBedmachine(X,Y,string,continent),

if strcmpi(continent,'antarctica'),
	ncdate='2014-09-23';
	morlighemnc=['/Users/larour/ModelData2/MMAntarctica2013/AntarcticaMCdataset-' ncdate '.nc'];
else
	%ncdate='2015-10-05';
	ncdate='2017-04-04';
	morlighemnc=['/Users/larour/ModelData2/MMGreenland2013/MCdataset-' ncdate '.nc'];
end

disp(['   -- BedMachine version: ' ncdate]);
xdata = double(ncread(morlighemnc,'x'));
ydata = double(ncread(morlighemnc,'y'));

offset=2;

xmin=min(X(:)); xmax=max(X(:));
posx=find(xdata<=xmax);
id1x=max(1,find(xdata>=xmin,1)-offset);
id2x=min(numel(xdata),posx(end)+offset);

ymin=min(Y(:)); ymax=max(Y(:));
posy=find(ydata>=ymin);
id1y=max(1,find(ydata<=ymax,1)-offset);
id2y=min(numel(ydata),posy(end)+offset);

disp(sprintf('   -- BedMachine %s : loading %s',continent,string));
data  = double(ncread(morlighemnc,string,[id1x id1y],[id2x-id1x+1 id2y-id1y+1],[1 1]))';
xdata=xdata(id1x:id2x);
ydata=ydata(id1y:id2y);
data(find(data==-999999))=NaN;

disp(sprintf('   -- BedMachine %s : interpolating %s',continent,string));
if strcmp(string,'mask') | strcmp(string,'source'),
	%Need nearest neighbor to avoid interpolation between 0 and 2
	if strcmpi(continent,'greenland'), ydata=flipud(ydata); data=flipud(data); end
	output = InterpFromGridToMesh(xdata,ydata,data,double(X),double(Y),0);
else
	ydata=flipud(ydata); data=flipud(data); 
	output = InterpFromGridToMesh(xdata,ydata,data,double(X),double(Y),0);
end

output(find(output==-999999))=NaN;

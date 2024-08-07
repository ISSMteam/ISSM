function md = googlemaps(md,ullat,ullon,lrlat,lrlon,varargin)
%GOOGLEMAPS - Extract image from Google maps for given region
%
%   Usage:
%       md = googlemaps(md)
%       md = googlemaps(md,zoomlevel)
%       md = googlemaps(md,ullat,ullon,lrlat,lrlon)
%       md = googlemaps(md,ullat,ullon,lrlat,lrlon,options)
%
%   - ullat,ullon: Upper Left corner latitude and longitude
%   - lrlat,lrlon: Lower Right corner latitude and longitude
%
%   Available options:
%      - zoomlevel: between 1 and 21 (default dynamically calculated)

%Parse inputs
if nargin<=5,
	options=pairoptions;
else
	options=varargin{:};
	if ~isa(options,'pairoptions'),
		options=pairoptions(varargin{:});
	end
end

%Check that no temp.* exists
if exist('temp.tiff','file'),
	error('File temp.tiff already exists, remove first');
end
if exist('temp.png','file'),
	error('File temp.png already exists, remove first');
end

if nargin==2,
	options=addfielddefault(options,'zoomlevel',ullat);
end

if md.mesh.epsg==0,
	error('md.mesh.epsg not defined');
end
if nargin<3,
	%Get xlim and ylim (used to extract Google maps image)
	xlim=getfieldvalue(options,'xlim',[min(md.mesh.x) max(md.mesh.x)]);
	ylim=getfieldvalue(options,'ylim',[min(md.mesh.y) max(md.mesh.y)]);
	if md.mesh.epsg==3413,
		[latlist lonlist]= xy2ll(...
			[linspace(xlim(1),xlim(2),100) linspace(xlim(2),xlim(2),100) linspace(xlim(2),xlim(1),100) linspace(xlim(1),xlim(1),100)],...
			[linspace(ylim(1),ylim(1),100) linspace(ylim(1),ylim(2),100) linspace(ylim(2),ylim(2),100) linspace(ylim(2),ylim(1),100)],...
			+1,45,70);
	elseif md.mesh.epsg==3031,
		[latlist lonlist]= xy2ll(...
			[linspace(xlim(1),xlim(2),100) linspace(xlim(2),xlim(2),100) linspace(xlim(2),xlim(1),100) linspace(xlim(1),xlim(1),100)],...
			[linspace(ylim(1),ylim(1),100) linspace(ylim(1),ylim(2),100) linspace(ylim(2),ylim(2),100) linspace(ylim(2),ylim(1),100)],...
			-1,0,71);
	elseif md.mesh.epsg==26906, %UTM 6V Columbia Glacier Alaska
		[latlist lonlist]= utm2ll(...
			[linspace(xlim(1),xlim(2),100) linspace(xlim(2),xlim(2),100) linspace(xlim(2),xlim(1),100) linspace(xlim(1),xlim(1),100)],...
			[linspace(ylim(1),ylim(1),100) linspace(ylim(1),ylim(2),100) linspace(ylim(2),ylim(2),100) linspace(ylim(2),ylim(1),100)],...
			6);
	elseif numel(md.mesh.lat)==numel(md.mesh.x),
		latlist = md.mesh.lat; %That might work?
		lonlist = md.mesh.long;
	else
		error('EPSG code not supported yet, and no lat long found in md.mesh');
	end

	%Image corners in lat/long
	ullat = max(latlist); ullon = min(lonlist);
	lrlat = min(latlist); lrlon = max(lonlist);
elseif nargin>1 & nargin<5,
	help googlemaps
	error('Wrong usage');
end

%Get region specific projection parameters
EPSGgoogle = 'EPSG:3785';   % Mercator       http://www.spatialreference.org/ref/epsg/3785/
EPSGlocal  = ['EPSG:' num2str(md.mesh.epsg)];

%Find optimal zoomlevel
if exist(options,'zoomlevel'),
	zoomlevel = getfieldvalue(options,'zoomlevel');
else
	zoomlevel = optimalzoomlevel(ullat,ullon,lrlat,lrlon);
	display(['googlemaps info: default zoomlevel level ' num2str(zoomlevel)]);
end
scale   = 1;
maxsize = 640;
bottom  = 50;

%Read Google Maps key
if exist('~/.googlemapskey');
	key=deblank(fileread('~/.googlemapskey'));
	iskey = true;
else
	%To create an API key, visit: https://developers.google.com/maps/documentation/maps-static/get-api-key
	%Once approved, go to Google Cloud Platform and enable "Maps Static APIs"
	%You will also need to "enable billing fot this project"
	%Then get the key and enter it in ~/.googlemapskey
	disp('It appears that you do not have a Google Maps API key, using MATLAB''s mapping toolbox instead');
	iskey = false;
end

if iskey
	%convert all these coordinates to pixels
	[ulx, uly]= latlontopixels(ullat, ullon, zoomlevel);
	[lrx, lry]= latlontopixels(lrlat, lrlon, zoomlevel);

	%calculate total pixel dimensions of final image
	dx = lrx - ulx;
	dy = uly - lry;

	%calculate rows and columns
	cols = ceil(dx/maxsize);
	rows = ceil(dy/(maxsize-bottom));

	%calculate pixel dimensions of each small image
	width   = ceil(dx/cols);
	height  = ceil(dy/rows);
	heightplus = height + bottom;

	%Initialize final image
	final = zeros(floor(dy),floor(dx),3);%RGB image
	for x=0:cols-1,
		for y=0:rows-1,
			dxn = width  * (0.5 + x);
			dyn = height * (0.5 + y);
			[latn, lonn] = pixelstolatlon(ulx + dxn, uly - dyn - bottom/2, zoomlevel);

			position = [num2str(latn) ',' num2str(lonn)];
			disp(['Google Earth tile: ' num2str(x) '/' num2str(cols-1) ' ' num2str(y) '/' num2str(rows-1) ' (center: ' position ')']);

			%Google maps API: http://developers.google.com/maps/documentation/staticmaps/
			params = [...
				'center=' position ...
				'&zoomlevel=' num2str(zoomlevel)...
				'&size=' num2str(width) 'x' num2str(heightplus)...
				'&maptype=satellite'...
				'&sensor=false'...
				'&scale=' num2str(scale)];
			if iskey,
				params = [params,'&key=' key];
			end
			url = ['http://maps.google.com/maps/api/staticmap?' params];
			count = 0;
			countmax = 10;
			while(true)
				try,
					[X, map]=imread(url,'png');
					break;
				catch me,
					disp(['Failed, trying again... (' num2str(countmax-count) ' more attempts)']);
					count = count+1;
					pause(.3);
					if count>countmax,
						disp('Giving up...');
						rethrow(me);
					end
				end
			end
			X=ind2rgb(X,map);

			indx1 = floor(x*width)+1;
			indx2 = min(floor(dx),floor(x*width)+size(X,2));
			indy1 = floor(y*height)+1;
			indy2 = min(floor(dy),floor(y*height)+size(X,1));
			final(indy1:indy2,indx1:indx2,:)=X(1:indy2-indy1+1,1:indx2-indx1+1,:);
		end
	end

	%prepare coordinate matrix of images
	[gX gY]=meshgrid(ulx:ulx+size(final,2)-1,uly:-1:uly-size(final,1)+1);
	[LAT LON]=pixelstolatlon(gX,gY, zoomlevel);
else
	% Read the basemap image from the web:
	fprintf(['Downloading image from readBasemapImage (zoomlevel=' num2str(zoomlevel) ')... ']);
	[final,R,attrib] = readBasemapImage('satellite', [lrlat ullat], [ullon lrlon], zoomlevel);
	fprintf('done!\n');

	%prepare coordinate matrix of images
	[X,Y] = meshgrid(...
		linspace(R.XWorldLimits(1),R.XWorldLimits(2),R.RasterSize(2)),...
		linspace(R.YWorldLimits(2),R.YWorldLimits(1),R.RasterSize(1)));

	%remove logo
	logopix = 9;
	final(end-logopix:end,:,:) = [];
	X(end-logopix:end,:) = [];
	Y(end-logopix:end,:) = [];

	[LAT LON] = mercator2ll(X, Y);
end

%Convert image Lat/Lon to X/Y
if md.mesh.epsg==3413
	[X Y]=ll2xy(LAT,LON,+1,45,70);
elseif md.mesh.epsg==3031
	[X Y]=ll2xy(LAT,LON,-1,0,71);
elseif md.mesh.epsg==4326
	X=LON;
	Y=LAT;
else
	error('EPSG code not supported yet');
end

%Write image
imwrite(final,'temp.png','png')
[ulmx ulmy]=ll2mercator(ullat,ullon);
[lrmx lrmy]=ll2mercator(lrlat,lrlon);

%Create Geotiff for Mercator projection 
[status,result] = system(['gdal_translate -of Gtiff -co "tfw=yes"  -a_ullr '...
	num2str(ulmx,'%15.8f') ' ' num2str(ulmy,'%15.8f') ' ' num2str(lrmx,'%15.8f') ' ' num2str(lrmy,'%15.8f')...
	' -a_srs "' EPSGgoogle '" "temp.png" "temp.tiff"']);
delete('temp.png');

%If not GDAL, exit
if status~=0,
	disp(result);
	disp('googlemaps info: GDAL not found or not working properly, the Google image will not be transformed');
	md.radaroverlay.pwr=final;
	md.radaroverlay.x=X;
	md.radaroverlay.y=Y;
	return
end

%reproject from mercator (EPSG:3785) to UPS Ant (EPSG:3031)
[status,result] = system(['gdalwarp  -s_srs ' EPSGgoogle ' -t_srs ' EPSGlocal ' temp.tiff temp2.tiff']);
delete('temp.tiff','temp.tfw');

%If previous command failed, exit
if ~isempty(strfind(result,'ERROR')),
	disp(result);
	disp(' ');disp('googlemaps info: GDAL not working properly (missing PROJ.4 library?), Google image will not be transformed');
	disp(result);
	md.radaroverlay.pwr=final;
	md.radaroverlay.x=X;
	md.radaroverlay.y=Y;
	return
end

%Put everything in model
[status output]=system('gdalinfo temp2.tiff | command grep "Upper Left"');
ul = sscanf(output,'Upper Left  (%f, %f)');
[status output]=system('gdalinfo temp2.tiff | command grep "Lower Right"');
lr = sscanf(output,'Lower Right (%f, %f)');
[status output]=system('gdalinfo temp2.tiff | command grep "Size is"');
si = sscanf(output,'Size is %i, %i');
x_m=linspace(ul(1),lr(1),si(1));
y_m=linspace(ul(2),lr(2),si(2)); %We need to reverse y_m because the image is read upside down by matlab
final=imread('temp2.tiff');
delete('temp2.tiff');

md.radaroverlay.pwr=final;
md.radaroverlay.x=x_m;
md.radaroverlay.y=y_m;
end
function [px py]=latlontopixels(lat, lon, zoomlevel),
	EARTH_RADIUS = 6378137;
	EQUATOR_CIRCUMFERENCE = 2 * pi * EARTH_RADIUS;
	INITIAL_RESOLUTION = EQUATOR_CIRCUMFERENCE / 256.0;
	ORIGIN_SHIFT = EQUATOR_CIRCUMFERENCE / 2.0;

	[mx,my]=ll2mercator(lat,lon);
	res = INITIAL_RESOLUTION / (2^zoomlevel);
	px = (mx + ORIGIN_SHIFT) / res;
	py = (my + ORIGIN_SHIFT) / res;
end

function [lat lon]=pixelstolatlon(px, py, zoomlevel),
	EARTH_RADIUS = 6378137;
	EQUATOR_CIRCUMFERENCE = 2 * pi * EARTH_RADIUS;
	INITIAL_RESOLUTION = EQUATOR_CIRCUMFERENCE / 256.0;
	ORIGIN_SHIFT = EQUATOR_CIRCUMFERENCE / 2.0;

	res = INITIAL_RESOLUTION / (2^zoomlevel);
	mx = px * res - ORIGIN_SHIFT;
	my = py * res - ORIGIN_SHIFT;
	[lat lon] = mercator2ll(mx,my);
end
function  zoomlevel = optimalzoomlevel(ullat,ullon,lrlat,lrlon)

	EARTH_RADIUS = 6378137;
	EQUATOR_CIRCUMFERENCE = 2 * pi * EARTH_RADIUS;
	INITIAL_RESOLUTION = EQUATOR_CIRCUMFERENCE / 256.0;

	optimalsize = 1000; %Number of pixels in final image

	[ulmx ulmy]=ll2mercator(ullat,ullon);
	[lrmx lrmy]=ll2mercator(lrlat,lrlon);
	distance = sqrt((lrmx-ulmx)^2 + (lrmy-ulmy)^2);

	zoomlevel1 = floor(log(INITIAL_RESOLUTION*optimalsize/(lrmx-ulmx))/log(2));
	zoomlevel2 = floor(log(INITIAL_RESOLUTION*optimalsize/(ulmy-lrmy))/log(2));

	zoomlevel=max(zoomlevel1,zoomlevel2);

	zoomlevel = min(max(1,zoomlevel),21);
end

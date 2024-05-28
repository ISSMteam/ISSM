function [output] = interpBedmap2(X,Y,string),
%INTERPBEDMAP2 - interpolate bedmap2 data
%
%   Available data:
%      1. bed                          is bed height
%      2. surface                      is surface height
%      3. thickness                    is ice thickness
%      4. icemask_grounded_and_shelves is a mask file showing the grounding line and the extent of the floating ice shelves
%      5. rockmask                     is a mask file showing rock outcrops
%      6. lakemask_vostok              is a mask file showing the extent of the lake cavity of Lake Vostok
%      7. grounded_bed_uncertainty     is the bed uncertainty grid shown in figure 12 of the manuscript
%      8. thickness_uncertainty_5km    is the thickness uncertainty grid shown in figure 11 of the manuscript
%      9. coverage                     is a binary grid showing the distribution of ice thickness data used in the grid of ice thickness
%     10. gl04c_geoid_to_wgs84         is the height conversion values (as floating point) used to convert from WGS84 datum heights to
%                                      g104c geoidal heights (to convert back to WGS84, add this grid)
%
%   Usage:
%      [dataout] = interpBedmap2(X,Y,string)

switch (oshostname()),
	case {'ronne'}
		nc = '/home/ModelData/Antarctica/BedMap2/bedmap2_bin/Bedmap2.nc';
	case {'totten'}
		nc = '/totten_1/ModelData/Antarctica/BedMap2/bedmap2_bin/Bedmap2.nc';
	otherwise
		error('hostname not supported yet');
end
if exist(nc,'file')
	if strcmp(string,'thickness_uncertainty_5km')
		xdata = double(ncread(nc,'x_5km'));
		ydata = double(ncread(nc,'y_5km'));
	else
		xdata = double(ncread(nc,'x'));
		ydata = double(ncread(nc,'y'));
	end

	offset=2;

	xmin=min(X(:)); xmax=max(X(:));
	posx=find(xdata<=xmax);
	id1x=max(1,find(xdata>=xmin,1)-offset);
	id2x=min(numel(xdata),posx(end)+offset);

	ymin=min(Y(:)); ymax=max(Y(:));
	posy=find(ydata>=ymin);
	id1y=max(1,find(ydata<=ymax,1)-offset);
	id2y=min(numel(ydata),posy(end)+offset);

	data  = double(ncread(nc,string,[id1x id1y],[id2x-id1x+1 id2y-id1y+1],[1 1]))';
	xdata=xdata(id1x:id2x);
	ydata=ydata(id1y:id2y);

	if ~strcmp(string,'coverage'),
		data(find(data==-9999))=NaN;
	end

	if strcmpi(string,'icemask_grounded_and_shelves') | strcmpi(string,'rockmask'),
		output = InterpFromGrid(xdata,ydata,data,double(X),double(Y),'nearest');
	else
		output = InterpFromGrid(xdata,ydata,data,double(X),double(Y)); % linear interpolation is default
	end
	return;

%For Eric's computer (using Binary files)
elseif exist('/Users/larour/ModelData/BedMap2/bedmap2_bin/','dir')
	% ================================  OLD ===============================================
	path='/Users/larour/ModelData/BedMap2/bedmap2_bin/'
	if strcmp(string,'gl04c_geoid_to_wgs84'),
		filepath = [path '/gl04c_geiod_to_wgs84.flt'];
	else
		filepath = [path '/bedmap2_' string '.flt'];
	end
	fid=fopen(filepath,'r','l');
	data=fread(fid,[6667,6667],'float32');
	fclose(fid);

	% define grid
	if strcmp(string,'thickness_uncertainty_5km'),
		ncols    =1361;
		nrows    =1361;
		xll      =-3401000;
		yll      =-3402000;
		gridsize =5000;
	else
		ncols    =6667;
		nrows    =6667;
		xll      =-3333000;
		yll      =-3333000;
		gridsize =1000;
	end
	x_m=xll+(0:1:ncols-1)'*gridsize;
	y_m=yll+(0:1:nrows-1)'*gridsize;

	%Change default to NaN
	if ~strcmp(string,'coverage'),
		data(find(data==-9999))=NaN;
	end

	%rotate 90 degrees clockwise
	data = rot90(data);

	%Interpolate
	if strcmpi(string,'icemask_grounded_and_shelves') | strcmpi(string,'rockmask'),
		dataout = InterpFromGrid(x_m,y_m,data,double(X),double(Y),'nearest');
	else
		dataout = InterpFromGrid(x_m,y_m,data,double(X),double(Y));
	end
else
	error('not supported');
end

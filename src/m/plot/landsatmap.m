% Explain
%  upload landsatmap to md.radaroverlay
%
% Usage
%  md = landsatmap(md);
%  md = landsatmap(md,'highres',1);
function md = landsatmap(md,varargin),

% check input variables
if nargin == 1 ,% {{{
	options = pairoptions;
else
	options = pairoptions(varargin{:});
end% }}}

%process mesh and data
[x y z elements is2d isplanet]=processmesh(md,[],options);

%check is2d
if ~is2d,
   error('buildgridded error message: gridded not supported for 3d meshes, project on a layer');
end

% get xlim, and ylim
xlim    = getfieldvalue(options,'xlim',[min(x) max(x)])/getfieldvalue(options,'unit',1);
ylim    = getfieldvalue(options,'ylim',[min(y) max(y)])/getfieldvalue(options,'unit',1);
highres = getfieldvalue(options,'highres',0);

if md.mesh.epsg == 3031 % Antarctica region {{{
	if highres, % high resolution geotiff file
		if 1, disp('   LIMA with geotiff'), % {{{
			disp('WARNING : this image shoud be collected with geocoded tif file');
			% find merged mosaic landsat image {{{
			limapath = {'/drive/project_inwoo/issm/Data/LIMA/AntarcticaLandsat.tif';
				'/home/DATA/ICESHEET/LIMA/AntarcticaLandsat.tif'};
			pos = zeros(length(limapath),1);
			for ii = 1:length(limapath)
				if exist(limapath{ii}), pos(ii) = 1; end
			end
			limapath = limapath{find(pos)};
			fprintf('   LIMA path is %s\n', limapath);
			% }}}

			% read image
			im = imread(limapath);

			% Region of LIMA data set
			info = gdalinfo(limapath); % get geotiff info
			xm = info.xmin + info.dx*[0:info.nx-1];
			ym = info.ymax - info.dy*[0:info.ny-1];

			%disp('   find region of model at LIMA');
			offset = 1e+4;
			posx = find((xm > xlim(1)-offset).* (xm < xlim(2)+offset));
			posy = find((ym > ylim(1)-offset).* (ym < ylim(2)+offset));
		end % }}}
	else
		if 1, disp('   LIMA with jp2'), % {{{
			% find merged mosaic landsat image {{{
			limapath = {'/drive/project_inwoo/issm/Data/LIMA/jp2_100pct/00000-20080314-144756363.jp2';
				'/data/project_inwoo/issm/Data/LIMA/jp2_100pct/00000-20080314-144756363.jp2';
				'/home/DATA/ICESHEET/LIMA/jp2_100pct/00000-20080314-144756363.jp2'};
			pos = zeros(length(limapath),1);
			for ii = 1:length(limapath)
				if exist(limapath{ii}), pos(ii) = 1; end
			end
			
			if sum(pos) == 0,
				fprintf('download website : https://lima.usgs.gov/fullcontinent.php\n');
				error('Landsat image at Antarctic region should be downloaded at above website');
			end
			limapath = limapath{find(pos)};
			fprintf('   LIMA path is %s\n', limapath);
			% }}}

			% Resolution and coordinates of upper left corner:
			xres = 240.010503438;
			yres = -240.000000516;
			xul = -2668154.98388;
			yul = 2362214.96998;

			% Arrays of pixel coordinates:
			xm = [xul:xres:2.813684914643920e+06];
			ym = [yul:yres:-2.294505040031947e+06];

			% reduction level 3 corresponds to very 2^3 = 8 points 
			rlevel = 2;
			xm = xm(1:2^rlevel:end);
			ym = ym(1:2^rlevel:end);
			im = imread(limapath,'reductionlevel',rlevel);

			%disp('   find region of model at LIMA');
			offset = 1e+4;
			posx = find((xm > xlim(1)-offset).* (xm < xlim(2)+offset));
			posy = find((ym > ylim(1)-offset).* (ym < ylim(2)+offset));
		end % }}}
		if 0, disp('   LIMA with reduced tiff'), % {{{
			% find merged mosaic landsat image {{{
			limapath = {'/drive/project_inwoo/issm/Data/LIMA/tiff_90pct/00000-20080319-092059124.tif'};
			pos = zeros(length(limapath),1);
			for ii = 1:length(limapath)
				if exist(limapath{ii}), pos(ii) = 1; end
			end
			
			if sum(pos) == 0,
				fprintf('download website : https://lima.usgs.gov/fullcontinent.php\n');
				error('Landsat image at Antarctic region should be downloaded at above website');
			end
			limapath = limapath{find(pos)};
			fprintf('   LIMA path is %s\n', limapath);
			% }}}

			% read image
			im = imread(limapath);

			% Region of LIMA data set
			info = gdalinfo(limapath); % get geotiff info
			xm = info.xmin + info.dx*[0:info.nx-1];
			ym = info.ymax - info.dy*[0:info.ny-1];

			%disp('   find region of model at LIMA');
			offset = 1e+4;
			posx = find((xm > xlim(1)-offset).* (xm < xlim(2)+offset));
			posy = find((ym > ylim(1)-offset).* (ym < ylim(2)+offset));
		end % }}}
	end

	% update region of radaroverlay
	md.radaroverlay.x = xm(posx);
	md.radaroverlay.y = ym(posy);
	md.radaroverlay.pwr = im(posy, posx,:);

	% }}}
else
	error('Check md.mesh.epsg, available LIMA regeion is at Antarctica (EPSG:3031)');
end


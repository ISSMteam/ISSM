function plot_landsat(md,data,options,plotlines,plotcols,i),
	% Explain
	%  This funtion loads Landsat Image Mosaic Antarctica (LIMA) for background image.
	%
	% Usage
	%  plot_landsat(md,data,options,plotlines,plotcols,i),
	%

	%process mesh and data
	[x2d y2d z2d elements2d is2d isplanet]=processmesh(md,[],options);
	[data datatype]=processdata(md,data,options);

	%check is2d
	if ~is2d,
		error('buildgridded error message: gridded not supported for 3d meshes, project on a layer');
	end

	%Get some options
	transparency = getfieldvalue(options,'transparency',.3);
	highres =  getfieldvalue(options,'highres',0);

	%Get xlim, and ylim
	xlim=getfieldvalue(options,'xlim',[min(x2d) max(x2d)])/getfieldvalue(options,'unit',1);
	ylim=getfieldvalue(options,'ylim',[min(y2d) max(y2d)])/getfieldvalue(options,'unit',1);

	pwr = md.radaroverlay.pwr;
	xm  = md.radaroverlay.x;
	ym  = md.radaroverlay.y;
	nx  = numel(xm);
	ny  = numel(ym);

	if md.mesh.epsg == 3031 & (isempty(md.radaroverlay.pwr) | isempty(md.radaroverlay.x) | isempty(md.radaroverlay.y) | length(size(md.radaroverlay.pwr)) < 3), % Antarctica region {{{
		if highres, 
			disp('   LIMA with geotiff'), % {{{
			disp('WARNING : this image shoud be collected with geocoded tif file');
			% find merged mosaic landsat image {{{
			limapath = {'/drive/project_inwoo/issm/Data/LIMA/AntarcticaLandsat.tif'};
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

			% find region of model at LIMA
			offset = 1e+4;
			posx = find((xm > xlim(1)-offset).* (xm < xlim(2)+offset));
			posy = find((ym > ylim(1)-offset).* (ym < ylim(2)+offset));
			% }}}
		else
			disp('   LIMA with reduced tiff'),
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

			% find region of model at LIMA
			offset = 1e+4;
			posx = find((xm > xlim(1)-offset).* (xm < xlim(2)+offset));
			posy = find((ym > ylim(1)-offset).* (ym < ylim(2)+offset));
		end

		% update region of radaroverlay
		md.radaroverlay.x = xm(posx);
		md.radaroverlay.y = ym(posy);
		md.radaroverlay.pwr = im(posy, posx,:);
		% }}}
	elseif length(size(md.radaroverlay.pwr)) == 3
		% it already contains LIMA image.
	elseif md.mesh.epsg == 3431 & (isempty(md.radaroverlay.pwr) | isempty(md.radaroverlay.x) | isempty(md.radaroverlay.y) | length(size(md.radaroverlay.pwr)) < 3), % Greenladn region 
		error('Greenland region is not yet available.');
	else
		error('Check md.mesh.epsg, available Landsat regeion is at Antarctica (EPSG:3031)');
	end

	%Check dataset
	assert(ndims(pwr) ~= 3, 'Error: Check ndims(md.radaroverlay.pwr) should be equal to 3.');
	assert(size(pwr) == [nx, ny, 3], 'Error: Given md.radaroverlay.pwr shoule be equal to (nx, ny, 3)');

	if any(diff(xm)) < 0
		disp('WARNING: md.radaroverlay.x should be increasing order.')
		xm = flip(xm);
		pwr= flip(pwr, 1);
	end
	if any(diff(ym)) < 0
		disp('WARNING: md.radaroverlay.y should be increasing order.')
		ym = flip(ym);
		pwr= flip(pwr, 2);
	end

	%Process image from model
	final = double(pwr)/double(max(md.radaroverlay.pwr(:))); %rescale between 0 and 1

	%Prepare grid
	if size(md.radaroverlay.x,1)==1 | size(md.radaroverlay.x,2)==1,
		data_grid=InterpFromMeshToGrid(elements2d,x2d/getfieldvalue(options,'unit',1),y2d/getfieldvalue(options,'unit',1),data,xm,ym,NaN);
		%data_grid=InterpFromMeshToGrid(md.mesh.elements,md.mesh.x/getfieldvalue(options,'unit',1),md.mesh.y/getfieldvalue(options,'unit',1),data,x_m,y_m,NaN);
	else
		X = md.radaroverlay.x;
		Y = md.radaroverlay.y;
		data_grid=InterpFromMeshToMesh2d(elements2d,x2d,y2d,data,X(:),Y(:),'default',NaN); data_grid=reshape(data_grid,size(X));
		%data_grid=InterpFromMeshToMesh2d(md.mesh.elements,md.mesh.x,md.mesh.y,data,X(:),Y(:),'default',NaN); data_grid=reshape(data_grid,size(X));
		xm=X(1,:); ym=Y(:,1);
	end

	data_nan=isnan(data_grid);
	if exist(options,'caxis'),
		caxis_opt=getfieldvalue(options,'caxis');
		data_grid(find(data_grid<caxis_opt(1)))=caxis_opt(1);
		data_grid(find(data_grid>caxis_opt(2)))=caxis_opt(2);
		data_min=caxis_opt(1);
		data_max=caxis_opt(2);
	else
		data_min=min(data_grid(:));
		data_max=max(data_grid(:));
	end
	colorm = getcolormap(options);
	image_rgb = ind2rgb(uint16((data_grid - data_min)*(length(colorm)/(data_max-data_min))),colorm);

	alpha=ones(size(data_grid));
	alpha(find(~data_nan))=transparency;
	alpha=repmat(alpha,[1 1 3]);

	final=alpha.*final+(1-alpha).*image_rgb;

	%Select plot area 
	subplotmodel(plotlines,plotcols,i,options);

	h=imagesc(xm*getfieldvalue(options,'unit',1),ym*getfieldvalue(options,'unit',1),final);

	%last step: mesh gridded?
	if exist(options,'edgecolor'),
		A=elements(:,1); B=elements(:,2); C=elements(:,3); 
		patch('Faces',[A B C],'Vertices', [x y z],'FaceVertexCData',data_grid(1)*ones(size(x)),'FaceColor','none','EdgeColor',getfieldvalue(options,'edgecolor'));
	end

	%Apply options
	if ~isnan(data_min),
		options=changefieldvalue(options,'caxis',[data_min data_max]); % force caxis so that the colorbar is ready
	end
	options=addfielddefault(options,'axis','xy equal off'); % default axis
	applyoptions(md,data,options);
	end

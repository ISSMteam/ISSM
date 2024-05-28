function diagnostics(md,varargin) 
%DIAGNOSTICS - output diagnostics for use in qgis
%
%
%   Usage:
%      diagnostics(md,options)
%   
%      where options include: 
%            'path': where are the diagnostics file output
%            'mesh': 0 or 1 (output mesh?, default 0)
%            'gl': 0 or 1 (output grounding line position?, default 0)
%            'vel': defaults to md.initialization.vel (can suplly another field if need be)
%                Suboptions include:
%                 'velposting': resolution at which the tiff file will be printed.
%            'hotspots' (0 or 1, default 0). spots where md.results.TransientSolution(i).Vel>threshold
%                Suboptions include:
%                 'threshold': vel value
%                 'hotspotsi': index into results.TransientSolution
%
%           
%
%   Examples:
%      diagnostics(md,'vel',md.initialization.vel,'gl',1,'hotspots',1,'hotspotsi',10,'threshold',500);


	%process options: 
	options=pairoptions(varargin{:});
	path=getfieldvalue(options,'path','./');

	%mesh:
	if getfieldvalue(options,'mesh',0),
		mesh2shp(md,[path '/mesh']);
	end

	%grounding line : 
	if getfieldvalue(options,'gl',0),
		contours=isoline(md,md.mask.ocean_levelset);
		expwrite(contours,[path '/groundingline.exp']);
		exp2shp([path '/groundingline.shp'],[path '/groundingline.exp']);
	end

	%velocity: 
	if exist(options,'vel'),
		vel=getfieldvalue(options,'vel',md.initialization.vel);
		xposting=getfieldvalue(options,'velposting',500);
		yposting=getfieldvalue(options,'velposting',500);

		xmin=min(md.mesh.x); ymax=max(md.mesh.y); 
		ncols=(max(md.mesh.x)-min(md.mesh.x))/xposting+1;
		nlines=(max(md.mesh.y)-min(md.mesh.y))/yposting+1;
		
		[xm,ym,vel]=InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,vel,xmin,ymax,xposting,yposting,nlines,ncols,0);
		vel=uint16(flipud(vel));

		imwrite(vel,[path '/veltemp.tif'],'tiff');

		string=sprintf('!gdal_translate -a_srs EPSG:3031 -a_ullr %g %g %g %g %s/veltemp.tif %s/vel.tif',...
		min(md.mesh.x),max(md.mesh.y),max(md.mesh.x),min(md.mesh.y),path,path);
		eval(string);
		delete([path '/veltemp.tif']);
	end

	%hot spots: 
	if getfieldvalue(options,'hotspots',0),
		threshold=getfieldvalue(options,'threshold',5000);
		i=getfieldvalue(options,'hotspotsi',length(md.results.TransientSolution));
		pos=find(md.results.TransientSolution(i).Vel>threshold);
		contour.x=md.mesh.x(pos); contour.y=md.mesh.y(pos); contour.density=1;
		expwrite(contour,[path '/hotspots.exp']);
		exp2shp([path '/hotspots.shp'],[path '/hotspots.exp']);
	end

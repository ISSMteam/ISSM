function mh=patchglobe(mh,mh2d,varargin)

	%process options: 
	options=pairoptions(varargin{:});

	%recover basic options:
	bandwidth=getfieldvalue(options,'bandwidth',100000);

	%some checks on the mesh: 
	if (isempty(mh.x) | isempty(mh.y) | isempty(mh.z) | isempty(mh.lat) | isempty(mh.long) | isempty(mh.r) ) 
		error('patchglobe error message: 3D planet mesh has one of the following empty: ''x'',''y'',''z'',''lat'',''long'' or ''r''');
	end
	if (isempty(mh2d.x) | isempty(mh2d.y) | isempty(mh2d.lat) | isempty(mh2d.long)) 
		error('patchglobe error message: 3D planet mesh has one of the following empty: ''x'',''y'',''lat'' or ''long''');
	end

	%give ourselves a unique temporary directory: 
	temproot=tempname; mkdir(temproot);

	%align external segments on 2d model: 
	mh2d.segments=alignsegments(mh2d.segments);

	%figure out domain outline: 
	meshtodomain(mh2d,[temproot '/Patch.exp']);

	%broaden this domain outline: 
	expcoarsen([temproot '/PatchBroad.exp'],[temproot '/Patch.exp'],200000);
	expcontract([temproot '/PatchBroad.exp'],[temproot '/PatchBroad.exp'],-bandwidth);

	%now flag vertices (from mh2d's broad contour that are on the global mesh: do this in the local 2d mesh reference system. 
	[x,y]=CoordTransform(mh.long,mh.lat,'EPSG:4326',['EPSG:' num2str(mh2d.epsg)]);
	flagsnods=ContourToNodes(x,y,[temproot '/PatchBroad.exp'],1);

	%expand flags to any element that touches the contour: 
	pos=find(sum(flagsnods(mh.elements),2));
	flags=zeros(mh.numberofelements,1); flags(pos)=1;

	%need to find the segment enveloppe of these elements:
	mh.segments=contourenvelope(mh,flags);

	%segments need to be ordered in line: 
	mh.segments=alignsegments(mh.segments);

	%x,y for segments: 
	[xsegs,ysegs]=CoordTransform(mh.long(mh.segments(:,1)),mh.lat(mh.segments(:,1)),'EPSG:4326',['EPSG:' num2str(mh2d.epsg)]);

	%create lat,long contour out of these segments:
	meshtodomain(mh,[temproot '/PatchEnveloppe.exp'],'latlong','on');
		
	%get these lat,long transformed to local mesh referencial:
	env=expread([temproot '/PatchEnveloppe.exp']); 
	[env.x,env.y]=CoordTransform(env.x,env.y,'EPSG:4326',['EPSG:' num2str(mh2d.epsg)]);

	%now, create domain outine from broad enveloppe and initial mesh 
	dom=expread([temproot '/Patch.exp']);
	
	%close the contours:
	env(1).x=[env(1).x;env(1).x(1)];
	env(1).y=[env(1).y;env(1).y(1)];
	dom(1).x=[dom(1).x;dom(1).x(1)];
	dom(1).y=[dom(1).y;dom(1).y(1)];

	%flip inner hole: 
	dom(1).x=flipud(dom(1).x);
	dom(1).y=flipud(dom(1).y);

	domain(1)=env; 
	domain(2)=dom;
	expwrite(domain,[temproot '/PatchBand.exp']);

	%plot mesh: 
	if  getfieldvalue(options,'plot',0), expdisp([temproot '/PatchBand.exp']); end

	%mesh: 
	mdb=bamg(model(),'domain',[temproot '/PatchBand.exp'],'MaxCornerAngle',1e-15,'gradation',10000,'NoBoundaryRefinment',1); 
	mhb=mdb.mesh; clear mdb;

	%double check: 
	if length(mhb.segments) ~= (length(dom(1).x)+length(env(1).x)-2),
		error('band mesh not consistent');
	end

	%augment patch with band
	[mhb.long,mhb.lat]=CoordTransform(mhb.x,mhb.y,['EPSG:' num2str(mh2d.epsg)],'EPSG:4326');
	mh2db=augment2dmesh(mh2d,mhb);

	%merge inner band and earth: 
	mh=mesh3dsurfaceplug2d(mh,mh2db,flags,mh.segments,xsegs,ysegs);

	%erase temporary directory: 
	system(['rm -rf ' temproot]);

end

function export_gl(md,varargin)

	templist=plotoptions(varargin{:}); 
	optionslist=templist.list;
	options=optionslist{1};
	options=checkplotoptions(md,options);

	%Setup unique directory in present dir: 
	directory=getfieldvalue(options,'directory','./');
	databasename=getfieldvalue(options,'database','webgl');

	%scaling factor: 
	scaling_factor=getfieldvalue(options,'scaling_factor',50);

	%Deal with title: 
	if exist(options,'title')
		title=getfieldvalue(options,'title');
	else
		title='';
	end

	%initialize model: 
	model.title=title;
	model.initialZoomFactor=getfieldvalue(options,'zoom',-.25);

	%Deal with contour {{{

	contour_lat1=md.mesh.lat(md.mesh.segments(:,1));
	contour_lat2=md.mesh.lat(md.mesh.segments(:,2));
	contour_long1=md.mesh.long(md.mesh.segments(:,1));
	contour_long2=md.mesh.long(md.mesh.segments(:,2));
	contour_surface1=md.geometry.surface(md.mesh.segments(:,1));
	contour_surface2=md.geometry.surface(md.mesh.segments(:,2));

	R1=6371000*ones(length(contour_surface1),1)+scaling_factor*contour_surface1;
	R2=6371000*ones(length(contour_surface2),1)+scaling_factor*contour_surface2;

	contourx1 = R1 .* cosd(contour_lat1) .* cosd(contour_long1);
	contoury1 = R1 .* cosd(contour_lat1) .* sind(contour_long1);
	contourz1 = R1 .* sind(contour_lat1);
	
	contourx2 = R2 .* cosd(contour_lat2) .* cosd(contour_long2);
	contoury2 = R2 .* cosd(contour_lat2) .* sind(contour_long2);
	contourz2 = R2 .* sind(contour_lat2);


	model.contourx1=contourx1;
	model.contoury1=contoury1;
	model.contourz1=contourz1;
	model.contourx2=contourx2;
	model.contoury2=contoury2;
	model.contourz2=contourz2;


	%}}}
%Deal with mesh and results {{{
	
	
	lat=md.mesh.lat;
	long=md.mesh.long;
	surface=md.geometry.surface;
	numberofelements=md.mesh.numberofelements;
	numberofvertices=md.mesh.numberofvertices;

	R=6371000*ones(numberofvertices,1)+scaling_factor*surface;

	x = R .* cosd(lat) .* cosd(long);
	y = R .* cosd(lat) .* sind(long);
	z = R .* sind(lat);


	%Deal with triangulation: 
	model.index=md.mesh.elements;
	model.x=x;
	model.y=y;
	model.z=z;
	model.surface=surface;
	
	%Deal with data: 
	results=struct([]);
	for i=1:length(optionslist),
		options=optionslist{i}; options=checkplotoptions(md,options);
		data=getfieldvalue(options,'data');
		results(i).data=data;
		results(i).caxis=getfieldvalue(options,'caxis',[min(data(:)) max(data(:))]);

		label=getfieldvalue(options,'label','');
		if strcmpi(label,''),
			%create generic label: 
			label=['data' num2str(i)];
		end
		results(i).label=label;

		shortlabel=getfieldvalue(options,'shortlabel','');
		if strcmpi(shortlabel,''),
			%create generic short label: 
			shortlabel=['data' num2str(i)];
		end
		results(i).shortlabel=shortlabel;
		
		if size(data,2)>1,
			time_range=getfieldvalue(options,'time_range',[0 100]);
			results(i).time_range=time_range;
		end

		unit=getfieldvalue(options,'unit','');
		if strcmpi(unit,''),
			%create generic unit: 
			unit='SI';
		end
		results(i).unit=unit;
	end
	model.results=results;
	
	%Write model to javascript database file: 
	writejsfile([directory databasename '.js'],model,databasename);
%}}}

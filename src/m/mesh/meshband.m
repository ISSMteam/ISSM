function band=meshband(mh,outerdomain,varargin)

	%process options: 
	options=pairoptions(varargin{:});

	%some checks on the mesh: 
	if (isempty(mh.x) | isempty(mh.y)  ) 
		error('meshband error message: mesh has one of the following empty: ''x'',''y''');
	end

	%give ourselves a unique temporary directory: 
	temproot=tempname; mkdir(temproot);

	%align external segments on internal mesh: 
	mh.segments=alignsegments(mh.segments);

	%figure out domain outline: 
	meshtodomain(mh,[temproot '/innerdomain.exp']);

	%create domain outine from inner and outer domain
	inner=expread([temproot '/innerdomain.exp']);
	if inner.closed==0,
		inner.nods=inner.nods+1;
		inner.x(end+1)=inner.x(1);
		inner.y(end+1)=inner.y(1);
		inner.closed=1;
		if getfieldvalue(options,'invert',0),
			inner.x=flipud(inner.x);
			inner.y=flipud(inner.y);
		end
	end

	[path,name,ext]=fileparts(outerdomain);
	if strcmpi(ext,'.shp'),
		outer=shpread(outerdomain);
	else
		outer=expread(outerdomain);
	end
	if outer.closed==0,
		outer.nods=outer.nods+1;
		outer.x(end+1)=outer.x(1);
		outer.y(end+1)=outer.y(1);
		outer.closed=1;
	end

	domain(1).x=outer.x; 
	domain(1).y=outer.y; 
	domain(1).density=outer.density; 
	domain(1).nods=outer.nods; 
	domain(1).closed=outer.closed;
	
	
	domain(2).x=inner.x; 
	domain(2).y=inner.y; 
	domain(2).density=inner.density; 
	domain(2).nods=inner.nods; 
	domain(2).closed=inner.closed;


	expwrite(domain,[temproot '/Band.exp']);
	
	if getfieldvalue(options,'plot',0),
		figure(1),clf,hold on,axis image;
		expdisp([temproot '/Band.exp'],'linestyle','k-*');
	end

	%mesh: 
	md=bamg(model(),'domain',[temproot '/Band.exp'],'MaxCornerAngle',1e-15,'gradation',10000); band=md.mesh; clear md;

	if getfieldvalue(options,'plot',0),
		figure(2),clf,trisurf(band.elements,band.x,band.y,band.x),view(2),shading faceted;
		hold on,expdisp([temproot '/Band.exp'],'linestyle','k-*');
	end

	%check that the domain vertices = the number of segments: 
	if abs(length(band.segments) - (length(domain(1).x)+length(domain(2).x)-2))>=2,
		disp(sprintf('band mesh not consistent: %i!=%i+%i\n',length(band.segments), length(domain(1).x), length(domain(2).x)));

		figure(3),clf,expdisp([temproot '/Band.exp'],'linestyle','r*');
		hold on; 
		for i=1:length(band.segments),
			i1=band.segments(i,1); i2=band.segments(i,2);
			plot([band.x(i1) band.x(i2)],[band.y(i1) band.y(i2)],'k*-');
			plot([band.x(i1)+ band.x(i2)]/2,[band.y(i1)+ band.y(i2)]/2,'g*');
		end
		%error('band mesh not consistent');
	end

	%erase temporary directory: 
	system(['rm -rf ' temproot]);

end

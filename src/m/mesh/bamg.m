function md=bamg(md,varargin)
%BAMG - mesh generation
%
%   Available options (for more details see ISSM website http://issm.jpl.nasa.gov/):
%
%   - domain :            followed by an ARGUS file that prescribes the domain outline
%   - holes :             followed by an ARGUS file that prescribes the holes
%   - subdomains :        followed by an ARGUS file that prescribes the list of
%                         subdomains (that need to be inside domain)
%
%   - hmin :              minimum edge length (default is 10^-100)
%   - hmax :              maximum edge length (default is 10^100)
%   - hVertices :         imposed edge length for each vertex (geometry or mesh)
%   - hminVertices :      minimum edge length for each vertex (mesh)
%   - hmaxVertices :      maximum edge length for each vertex (mesh)
%
%   - anisomax :          maximum ratio between the smallest and largest edges (default is 10^30)
%   - coeff :             coefficient applied to the metric (2-> twice as many elements, default is 1)
%   - cutoff :            scalar used to compute the metric when metric type 2 or 3 are applied
%   - err :               error used to generate the metric from a field
%   - errg :              geometric error (default is 0.1)
%   - field :             field of the model that will be used to compute the metric
%                         to apply several fields, use one column per field
%   - gradation :         maximum ratio between two adjacent edges
%   - Hessiantype :       0 -> use double L2 projection (default)
%                         1 -> use Green formula
%   - KeepVertices :      try to keep initial vertices when adaptation is done on an existing mesh (default 1)
%   - NoBoundaryRefinement: do not refine boundary, only follow contour provided (default 0). Allow subdomain boundary refinement though.
%   - NoBoundaryRefinementAllBoundaries: do not refine boundary, only follow contour provided (default 0)
%   - maxnbv :            maximum number of vertices used to allocate memory (default is 10^6)
%   - maxsubdiv :         maximum subdivision of existing elements (default is 10)
%   - metric :            matrix (numberofnodes x 3) used as a metric
%   - Metrictype :        0 -> absolute error          c/(err coeff^2) * Abs(H)        (default)
%                         1 -> relative error          c/(err coeff^2) * Abs(H)/max(s,cutoff*max(s))
%                         2 -> rescaled absolute error c/(err coeff^2) * Abs(H)/(smax-smin)
%   - nbjacoby :          correction used by Hessiantype=1 (default is 1)
%   - nbsmooth :          number of metric smoothing procedure (default is 3)
%   - omega :             relaxation parameter of the smoothing procedure (default is 1.8)
%   - power :             power applied to the metric (default is 1)
%   - splitcorners :      split triangles which have 3 vertices on the outline (default is 1)
%   - verbose :           level of verbosity (default is 1)
%
%   - vertical :          is this a 2d vertical mesh (flowband, default is 0)
%   - rifts :             followed by an ARGUS file that prescribes the rifts
%   - toltip :            tolerance to move tip on an existing point of the domain outline
%   - tracks :            followed by an ARGUS file that prescribes the tracks that the mesh will stick to
%   - RequiredVertices :  mesh vertices that are required. [x,y,ref]; ref is optional
%   - tol :               if the distance between 2 points of the domain outline is less than tol, they
%                         will be merged
%
%   Examples:
%      md=bamg(md,'domain','DomainOutline.exp','hmax',3000);
%      md=bamg(md,'field',[md.inversion.vel_obs md.geometry.thickness],'hmax',20000,'hmin',1000);
%      md=bamg(md,'metric',A,'hmin',1000,'hmax',20000,'gradation',3,'anisomax',1);

%process options
options=pairoptions(varargin{:});
options=deleteduplicates(options,1);

%initialize the structures required as input of BAMG
bamg_options=struct();
bamg_geometry=bamggeom();
bamg_mesh=bamgmesh();

subdomain_ref = 1;
hole_ref = 1;
% BAMG Geometry parameters {{{
if exist(options,'domain'),

	%Check that file exists
	domainfile=getfieldvalue(options,'domain');
	if ischar(domainfile),
		if ~exist(domainfile,'file') error(['bamg error message: file ' domainfile ' not found']); end

		%read domain according to its extension: 
		[path,name,ext]=fileparts(domainfile);
		if strcmp(ext,'.exp'),
			domain=expread(domainfile);
		elseif strcmp(ext,'.shp'),
			domain=shpread(domainfile);
		else
			error(['bamg error message: file ' domainfile ' format not supported (.exp or .shp)']);
		end
	elseif isstruct(domainfile),
		domain = domainfile;
	else
		error('''domain'' type not supported yet');
	end

	holes = [];
	if exist(options,'holes'),
		holesfile=getfieldvalue(options,'holes');
		if ischar(holesfile),
			if ~exist(holesfile,'file') error(['bamg error message: file ' holesfile ' not found']); end

			%read holes according to its extension: 
			[path,name,ext]=fileparts(holesfile);
			if strcmp(ext,'.exp'),
				holes=expread(holesfile);
			elseif strcmp(ext,'.shp'),
				holes=shpread(holesfile);
			else
				error(['bamg error message: file ' holesfile ' format not supported (.shp or .exp)']);
			end
		elseif isstruct(holesfile),
			holes = holesfile;
		else
			error('''holes'' type not supported yet');
		end
	end
	subdomains = [];
	if exist(options,'subdomains'),
		subdomainsfile=getfieldvalue(options,'subdomains');
		if ischar(subdomainsfile),
			if ~exist(subdomainsfile,'file') error(['bamg error message: file ' subdomainsfile ' not found']); end

			%read subdomains according to its extension: 
			[path,name,ext]=fileparts(subdomainsfile);
			if strcmp(ext,'.exp'),
				subdomains=expread(subdomainsfile);
			elseif strcmp(ext,'.shp'),
				subdomains=shpread(subdomainsfile);
			else
				error(['bamg error message: file ' subdomainsfile ' format not supported (.shp or .exp)']);
			end
		elseif isstruct(subdomainsfile),
			subdomains = subdomainsfile;
		else
			error('''subdomains'' type not supported yet');
		end
	end

	%Build geometry 
	count=0;
	for i=1:length(domain),

		%Check that the domain is closed
		if (domain(i).x(1)~=domain(i).x(end) | domain(i).y(1)~=domain(i).y(end)),
			error('bamg error message: all contours provided in ''domain'' should be closed');
		end

		%Check that all holes are INSIDE the principle domain outline
		if i>1,
			flags=ContourToNodes(domain(i).x,domain(i).y,domain(1),0);
			if any(~flags),
				error('bamg error message: All holes should be strictly inside the principal domain');
			end
		end

		%Check orientation
		nods=domain(i).nods-1; %the domain is closed (domain[1] = domain[end])
		test = sum([(domain(i).x(2:nods+1) - domain(i).x(1:nods)).*(domain(i).y(2:nods+1) + domain(i).y(1:nods))]);
		if (i==1 && test>0) || (i>1 && test<0),
			disp('At least one contour was not correctly oriented and has been re-oriented');
			domain(i).x = flipud(domain(i).x); domain(i).y = flipud(domain(i).y);
		end

		%Flag how many edges we have so far: 
		edge_length=length(bamg_geometry.Edges);

		%Add all points to bamg_geometry
		bamg_geometry.Vertices=[bamg_geometry.Vertices; [domain(i).x(1:nods) domain(i).y(1:nods) ones(nods,1)]];
		bamg_geometry.Edges   =[bamg_geometry.Edges;    [transpose(count+1:count+nods) transpose([count+2:count+nods count+1])  1*ones(nods,1)]];

		% Flag how many edges we have now, that way we know which edges belong 
		% to the domain. Will be used later for required edges if 
		% NoBoundaryRefinement equals 1.
		new_edge_length=length(bamg_geometry.Edges);
		edges_required=(edge_length+1):new_edge_length;

		if i>1,
			bamg_geometry.SubDomains=[bamg_geometry.SubDomains; 2 count+1 1 -subdomain_ref]; 
			subdomain_ref = subdomain_ref+1;
		else
			bamg_geometry.SubDomains=[bamg_geometry.SubDomains; 2 count+1 1 0];
		end

		%update counter
		count=count+nods;
	end
	for i=1:length(holes),

		%Check that the hole is closed
		if (holes(i).x(1)~=holes(i).x(end) | holes(i).y(1)~=holes(i).y(end)),
			error('bamg error message: all contours provided in ''holes'' should be closed');
		end

		%Check that all holes are INSIDE the principal domain (principal domain should be index 0)
		flags=ContourToNodes(holes(i).x,holes(i).y,domain(1),0);
		if any(~flags), error('bamg error message: All holes should be strictly inside the principal domain'); end

		%Check that hole is correctly oriented
		nods=holes(i).nods-1; %the hole is closed (hole[1] = hole[end])
		if(sum([(holes(i).x(2:nods+1) - holes(i).x(1:nods)).*(holes(i).y(2:nods+1) + holes(i).y(1:nods))]))<0
			disp('At least one hole was not correctly oriented and has been re-oriented');
			holes(i).x = flipud(holes(i).x); holes(i).y = flipud(holes(i).y);
		end

		%Add all points to bamg_geometry
		bamg_geometry.Vertices=[bamg_geometry.Vertices; [holes(i).x(1:nods) holes(i).y(1:nods) ones(nods,1)]];
		bamg_geometry.Edges   =[bamg_geometry.Edges;    [transpose(count+1:count+nods) transpose([count+2:count+nods count+1])  1.*ones(nods,1)]];
		bamg_geometry.SubDomains=[bamg_geometry.SubDomains; 2 count+1 1 -hole_ref]; hole_ref = hole_ref+1;

		%update counter
		count=count+nods;
	end
	for i=1:length(subdomains),

		%Check that the subdomain is closed
		if (subdomains(i).x(1)~=subdomains(i).x(end) | subdomains(i).y(1)~=subdomains(i).y(end)),
			error('bamg error message: all contours provided in ''subdomains'' should be closed');
		end

		%Checks that all subdomains are INSIDE the principal domain (principal domain should be index 0)
		flags=ContourToNodes(subdomains(i).x,subdomains(i).y,domain(1),0);
		if any(~flags),
			error('bamg error message: All subdomains should be strictly inside the principal domain');
		end

		%Check that subdomain is correctly oriented
		nods=subdomains(i).nods-1; % the subdomains are closed (subdomains[1] = subdomains[end])
		if(sum([(subdomains(i).x(2:nods+1) - subdomains(i).x(1:nods)).*(subdomains(i).y(2:nods+1) + subdomains(i).y(1:nods))]))>0
			disp('At least one subdomain was not correctly oriented and has been re-oriented');
			subdomains(i).x = flipud(subdomains(i).x); 
			subdomains(i).y = flipud(subdomains(i).y);
		end

		%Add all points to bamg_geometry
		bamg_geometry.Vertices=[bamg_geometry.Vertices; [subdomains(i).x(1:nods) subdomains(i).y(1:nods) ones(nods,1)]];
		bamg_geometry.Edges   =[bamg_geometry.Edges;    [transpose(count+1:count+nods) transpose([count+2:count+nods count+1])  1.*ones(nods,1)]];
		
		bamg_geometry.SubDomains=[bamg_geometry.SubDomains; 2 count+1 1 subdomain_ref];
		subdomain_ref = subdomain_ref+1;
		
		%update counter
		count=count+nods;
	end

	if getfieldvalue(options,'vertical',0),
		if numel(getfieldvalue(options,'Markers',[]))~=size(bamg_geometry.Edges,1),
			error(['for 2d vertical mesh, ''Markers'' option is required, and should be of size ' num2str(size(bamg_geometry.Edges,1))]);
		end
	end
	if numel(getfieldvalue(options,'Markers',[]))==size(bamg_geometry.Edges,1),
		bamg_geometry.Edges(:,3)=getfieldvalue(options,'Markers');
	end

	%take care of rifts
	if exist(options,'rifts'),
		%read rift file according to its extension: 
		riftfile=getfieldvalue(options,'rifts');
		[path,name,ext]=fileparts(riftfile);
		if strcmp(ext,'.exp'),
			rift=expread(riftfile);
		elseif strcmp(ext,'.shp'),
			rift=shpread(riftfile);
		else
			error(['bamg error message: file ' riftfile ' format not supported (.exp or .shp)']);
		end

		for i=1:length(rift),

			%detect whether all points of the rift are inside the domain
			flags=ContourToNodes(rift(i).x,rift(i).y,domain(1),0);
			if ~flags,
				error('one rift has all its points outside of the domain outline'),

			elseif any(~flags),
				%We have LOTS of work to do
				disp('Rift tip outside of or on the domain has been detected and is being processed...');

				%check that only one point is outside (for now)
				if sum(~flags)~=1,
					error('bamg error message: only one point outside of the domain is supported at this time');
				end

				%Move tip outside to the first position
				if flags(1)==0,
					%OK, first point is outside (do nothing),
				elseif (flags(end)==0),
					rift(i).x=flipud(rift(i).x);
					rift(i).y=flipud(rift(i).y);
				else
					error('bamg error message: only a rift tip can be outside of the domain');
				end

				%Get coordinate of intersection point
				x1=rift(i).x(1);
				y1=rift(i).y(1);
				x2=rift(i).x(2);
				y2=rift(i).y(2);
				for j=1:length(domain(1).x)-1;
					if SegIntersect([x1 y1; x2 y2],[domain(1).x(j) domain(1).y(j); domain(1).x(j+1) domain(1).y(j+1)]),

						%Get position of the two nodes of the edge in domain
						i1=j;
						i2=j+1;

						%rift is crossing edge [i1 i2] of the domain
						%Get coordinate of intersection point (http://mathworld.wolfram.com/Line-LineIntersection.html)
						x3=domain(1).x(i1); y3=domain(1).y(i1);
						x4=domain(1).x(i2); y4=domain(1).y(i2);
						x=det([det([x1 y1; x2 y2])  x1-x2;det([x3 y3; x4 y4])  x3-x4])/det([x1-x2 y1-y2;x3-x4 y3-y4]);
						y=det([det([x1 y1; x2 y2])  y1-y2;det([x3 y3; x4 y4])  y3-y4])/det([x1-x2 y1-y2;x3-x4 y3-y4]);

						segdis= sqrt((x4-x3)^2+(y4-y3)^2);
						tipdis=[sqrt((x-x3)^2+(y-y3)^2)  sqrt((x-x4)^2+(y-y4)^2)];

						if (min(tipdis)/segdis) < getfieldvalue(options,'toltip',0),
							disp('moving tip-domain intersection point');

							%Get position of the closer point
							if tipdis(1)>tipdis(2),
								pos=i2;
							else
								pos=i1;
							end

							%This point is only in Vertices (number pos).
							%OK, now we can add our own rift
							nods=rift(i).nods-1;
							bamg_geometry.Vertices=[bamg_geometry.Vertices; [rift(i).x(2:end) rift(i).y(2:end) ones(nods,1)]];
							bamg_geometry.Edges=[bamg_geometry.Edges;...
								pos count+1  (1+i);...
								[transpose(count+1:count+nods-1) transpose(count+2:count+nods)  (1+i)*ones(nods-1,1)]];
							count=count+nods;

							break;

						else
							%Add intersection point to Vertices
							bamg_geometry.Vertices=[bamg_geometry.Vertices; x y 1];
							count=count+1;

							%Decompose the crossing edge into 2 subedges
							pos=find(bamg_geometry.Edges(:,1)==i1 & bamg_geometry.Edges(:,2)==i2);
							if isempty(pos) error('bamg error message: a problem occurred...'); end
							bamg_geometry.Edges=[bamg_geometry.Edges(1:pos-1,:);...
								bamg_geometry.Edges(pos,1) count                      bamg_geometry.Edges(pos,3);...
								count                      bamg_geometry.Edges(pos,2) bamg_geometry.Edges(pos,3);...
								bamg_geometry.Edges(pos+1:end,:)];

							%OK, now we can add our own rift
							nods=rift(i).nods-1;
							bamg_geometry.Vertices=[bamg_geometry.Vertices; [rift(i).x(2:end) rift(i).y(2:end) ones(nods,1)]];
							bamg_geometry.Edges=[bamg_geometry.Edges;...
								count  count+1  2 ;...
								[transpose(count+1:count+nods-1) transpose(count+2:count+nods)  (1+i)*ones(nods-1,1)]];
							count=count+nods;

							break;
						end
					end
				end
			else
				nods=rift(i).nods-1;
				bamg_geometry.Vertices=[bamg_geometry.Vertices; [rift(i).x(:) rift(i).y(:) ones(nods+1,1)]];
				bamg_geometry.Edges=[bamg_geometry.Edges; [transpose(count+1:count+nods) transpose(count+2:count+nods+1)  (1+i)*ones(nods,1)]];
				count=count+nods+1;
			end
		end
	end

	%Deal with tracks
	if exist(options,'tracks'),
		%read tracks
		intrack=getfieldvalue(options,'tracks');
		if all(ischar(intrack)),
			intrack=expread(intrack);
        else
            intrack=double(intrack); %for some reason, it is of class "single"
            intrack=struct('x',intrack(:,1),'y',intrack(:,2));
            intrack=repmat(intrack, [1,1]);
        end
        
        %Process multiple tracks and handle domain intersections
        track=[];
        edges=[];
        newcount=count;
        for i=1:length(intrack), 
            nods=size(intrack(i).x,1);
            newtrack=[intrack(i).x intrack(i).y];
            newedges=[transpose(newcount+1:newcount+nods-1) transpose(newcount+2:newcount+nods) 3.*ones(nods-1,1)];

            %only keep those inside
            flags=ContourToNodes(newtrack(:,1),newtrack(:,2),domainfile,0);
            edgeflags = flags(1:end-1) & flags(2:end);

            %calculate edge offset, accounting for broken tracks
            offsets=zeros(size(edgeflags,1)+1,3);
            for j=2:size(offsets,1),
                if edgeflags(j-1)==0,
                    offsets(j,:)=offsets(j-1,:)+[1 1 0];
                else
                    offsets(j,:)=offsets(j-1,:);
                end
            end
            newedges=newedges-offsets(2:end,:);

            %add track segments
            newtrack=newtrack(find(flags),:);
            newedges=newedges(find(edgeflags),:);
            nods=size(newtrack,1);
            track=[track; newtrack];
            edges=[edges; newedges];
            newcount=newcount+nods;
        end
		if(size(track,2)==2), track=[track 3.*ones(size(track,1),1)]; end
        
		%Add all points to bamg_geometry
		nods=size(track,1);
		bamg_geometry.Vertices=[bamg_geometry.Vertices; track];
		bamg_geometry.Edges=[bamg_geometry.Edges; edges];

		%update counter
		count=count+nods;
	end

	%Deal with vertices that need to be kept by mesher
	if exist(options,'RequiredVertices'),

		%recover RequiredVertices
		requiredvertices=double(getfieldvalue(options,'RequiredVertices')); %for some reason, it is of class "single"
		if(size(requiredvertices,2)==2), requiredvertices=[requiredvertices 4.*ones(size(requiredvertices,1),1)]; end

		%only keep those inside
		flags=ContourToNodes(requiredvertices(:,1),requiredvertices(:,2),domain(1),0);
		requiredvertices=requiredvertices(find(flags),:);

		%Add all points to bamg_geometry
		nods=size(requiredvertices,1);
		bamg_geometry.Vertices=[bamg_geometry.Vertices; requiredvertices];

		%update counter
		count=count+nods;

	end

	%Deal with RequiredEdges
	if getfieldvalue(options,'NoBoundaryRefinement',0)==1,
		bamg_geometry.RequiredEdges=edges_required';
	elseif getfieldvalue(options,'NoBoundaryRefinementAllBoundaries',0)==1,
		bamg_geometry.RequiredEdges=[1:size(bamg_geometry.Edges,1)]';
	end

	%process geom
	%bamg_geometry=processgeometry(bamg_geometry,getfieldvalue(options,'tol',NaN),domain(1));

elseif isstruct(md.private.bamg) & isfield(md.private.bamg,'geometry'),
	bamg_geometry=bamggeom(md.private.bamg.geometry); 
else
	%do nothing...
end
%}}}
% BAMG Mesh parameters {{{
if (~exist(options,'domain') & md.mesh.numberofvertices~=0 & strcmp(elementtype(md.mesh),'Tria')),

	if isstruct(md.private.bamg) & isfield(md.private.bamg,'mesh'),
		bamg_mesh=bamgmesh(md.private.bamg.mesh);
	else
		bamg_mesh.Vertices=[md.mesh.x md.mesh.y ones(md.mesh.numberofvertices,1)];
		bamg_mesh.Triangles=[md.mesh.elements ones(md.mesh.numberofelements,1)];
	end

	if isstruct(md.rifts.riftstruct)
		error('bamg error message: rifts not supported yet. Do meshprocessrift AFTER bamg');
	end
end
%}}}
% BAMG Options {{{
bamg_options.Crack=getfieldvalue(options,'Crack',0);
bamg_options.anisomax=getfieldvalue(options,'anisomax',1e30);
bamg_options.coeff=getfieldvalue(options,'coeff',1.);
bamg_options.cutoff=getfieldvalue(options,'cutoff',1e-5);
bamg_options.err=getfieldvalue(options,'err',0.01);
bamg_options.errg=getfieldvalue(options,'errg',0.1);
bamg_options.field=getfieldvalue(options,'field',[]);
bamg_options.gradation=getfieldvalue(options,'gradation',1.5);
bamg_options.Hessiantype=getfieldvalue(options,'Hessiantype',0);
bamg_options.hmin=getfieldvalue(options,'hmin',1e-100);
bamg_options.hmax=getfieldvalue(options,'hmax',1e100);
bamg_options.hminVertices=getfieldvalue(options,'hminVertices',[]);
bamg_options.hmaxVertices=getfieldvalue(options,'hmaxVertices',[]);
bamg_options.hVertices=getfieldvalue(options,'hVertices',[]);
bamg_options.KeepVertices=getfieldvalue(options,'KeepVertices',1);
bamg_options.maxnbv=getfieldvalue(options,'maxnbv',1e6);
bamg_options.maxsubdiv=getfieldvalue(options,'maxsubdiv',10.);
bamg_options.metric=getfieldvalue(options,'metric',[]);
bamg_options.Metrictype=getfieldvalue(options,'Metrictype',0);
bamg_options.nbjacobi=getfieldvalue(options,'nbjacobi',1);
bamg_options.nbsmooth=getfieldvalue(options,'nbsmooth',3);
bamg_options.omega=getfieldvalue(options,'omega',1.8);
bamg_options.power=getfieldvalue(options,'power',1.);
bamg_options.splitcorners=getfieldvalue(options,'splitcorners',1);
bamg_options.verbose=getfieldvalue(options,'verbose',1);
%}}}

%call BAMG
[bamgmesh_out bamggeom_out]=BamgMesher(bamg_mesh,bamg_geometry,bamg_options);

if getfieldvalue(options,'vertical',0),
	md.mesh=mesh2dvertical();
	md.mesh.x=bamgmesh_out.Vertices(:,1);
	md.mesh.y=bamgmesh_out.Vertices(:,2);
	md.mesh.elements=bamgmesh_out.Triangles(:,1:3);
	md.mesh.edges=bamgmesh_out.IssmEdges;
	md.mesh.segments=bamgmesh_out.IssmSegments(:,1:3);
	md.mesh.segmentmarkers=bamgmesh_out.IssmSegments(:,4);

	%Fill in rest of fields:
	md.mesh.numberofelements=size(md.mesh.elements,1);
	md.mesh.numberofvertices=length(md.mesh.x);
	md.mesh.numberofedges=size(md.mesh.edges,1);
	md.mesh.vertexonboundary=zeros(md.mesh.numberofvertices,1);
	md.mesh.vertexonboundary(md.mesh.segments(:,1:2))=1;

elseif getfieldvalue(options,'3dsurface',0),
	
	md.mesh=mesh3dsurface();
	md.mesh.x=bamgmesh_out.Vertices(:,1);
	md.mesh.y=bamgmesh_out.Vertices(:,2);
	md.mesh.z=md.mesh.x; md.mesh.z(:)=0;
	md.mesh.elements=bamgmesh_out.Triangles(:,1:3);
	md.mesh.edges=bamgmesh_out.IssmEdges;
	md.mesh.segments=bamgmesh_out.IssmSegments(:,1:3);
	md.mesh.segmentmarkers=bamgmesh_out.IssmSegments(:,4);

	%Fill in rest of fields:
	md.mesh.numberofelements=size(md.mesh.elements,1);
	md.mesh.numberofvertices=length(md.mesh.x);
	md.mesh.numberofedges=size(md.mesh.edges,1);
	md.mesh.vertexonboundary=zeros(md.mesh.numberofvertices,1);
	md.mesh.vertexonboundary(md.mesh.segments(:,1:2))=1;

else 
	md.mesh=mesh2d();
	md.mesh.x=bamgmesh_out.Vertices(:,1);
	md.mesh.y=bamgmesh_out.Vertices(:,2);
	md.mesh.elements=bamgmesh_out.Triangles(:,1:3);
	md.mesh.edges=bamgmesh_out.IssmEdges;
	md.mesh.segments=bamgmesh_out.IssmSegments(:,1:3);
	md.mesh.segmentmarkers=bamgmesh_out.IssmSegments(:,4);

	%Fill in rest of fields:
	md.mesh.numberofelements=size(md.mesh.elements,1);
	md.mesh.numberofvertices=length(md.mesh.x);
	md.mesh.numberofedges=size(md.mesh.edges,1);
	md.mesh.vertexonboundary=zeros(md.mesh.numberofvertices,1);
	md.mesh.vertexonboundary(md.mesh.segments(:,1:2))=1;
end

%BAMG private fields
md.private.bamg=struct();
md.private.bamg.mesh=bamgmesh(bamgmesh_out);
md.private.bamg.geometry=bamggeom(bamggeom_out);
md.mesh.elementconnectivity=md.private.bamg.mesh.ElementConnectivity;
md.mesh.elementconnectivity(find(isnan(md.mesh.elementconnectivity)))=0;

%Check for orphan
if any(~ismember(1:md.mesh.numberofvertices,sort(unique(reshape(md.mesh.elements,3*md.mesh.numberofelements,1)))))
	error('Output mesh has orphans. Check your Domain and/or RequiredVertices');
end
end 

function geom=processgeometry(geom,tol,outline) % {{{

%Deal with edges
disp('Checking Edge crossing...');
i=0;
while (i<size(geom.Edges,1)),

	%edge counter
	i=i+1;

	%Get coordinates
	x1=geom.Vertices(geom.Edges(i,1),1);
	y1=geom.Vertices(geom.Edges(i,1),2);
	x2=geom.Vertices(geom.Edges(i,2),1);
	y2=geom.Vertices(geom.Edges(i,2),2);
	color1=geom.Edges(i,3);

	j=i; %test edges located AFTER i only
	while (j<size(geom.Edges,1)),

		%edge counter
		j=j+1;

		%Skip if the two edges already have a vertex in common
		if any(ismember(geom.Edges(i,1:2),geom.Edges(j,1:2))),
			continue
		end

		%Get coordinates
		x3=geom.Vertices(geom.Edges(j,1),1);
		y3=geom.Vertices(geom.Edges(j,1),2);
		x4=geom.Vertices(geom.Edges(j,2),1);
		y4=geom.Vertices(geom.Edges(j,2),2);
		color2=geom.Edges(j,3);

		%Check if the two edges are crossing one another
		if SegIntersect([x1 y1; x2 y2],[x3 y3; x4 y4]),

			%Get coordinate of intersection point (http://mathworld.wolfram.com/Line-LineIntersection.html)
			x=det([det([x1 y1; x2 y2])  x1-x2;det([x3 y3; x4 y4])  x3-x4])/det([x1-x2 y1-y2;x3-x4 y3-y4]);
			y=det([det([x1 y1; x2 y2])  y1-y2;det([x3 y3; x4 y4])  y3-y4])/det([x1-x2 y1-y2;x3-x4 y3-y4]);

			%Add vertex to the list of vertices
			geom.Vertices(end+1,:)=[x y min(color1,color2)];
			id=size(geom.Vertices,1);

			%Update edges i and j
			edgei=geom.Edges(i,:);
			edgej=geom.Edges(j,:);
			geom.Edges(i,:)    =[edgei(1) id       edgei(3)];
			geom.Edges(end+1,:)=[id       edgei(2) edgei(3)];
			geom.Edges(j,:)    =[edgej(1) id       edgej(3)];
			geom.Edges(end+1,:)=[id       edgej(2) edgej(3)];

			%update current edge second tip
			x2=x; y2=y;
		end
	end

end

%Check point outside
disp('Checking for points outside the domain...');
i=0;
num=0;
while (i<size(geom.Vertices,1)),

	%vertex counter
	i=i+1;

	%Get coordinates
	x=geom.Vertices(i,1);
	y=geom.Vertices(i,2);
	color=geom.Vertices(i,3);

	%Check that the point is inside the domain
	if (color~=1 & ~ContourToNodes(x,y,outline(1),1)),

		%Remove points from list of Vertices
		num=num+1;
		geom.Vertices(i,:)=[];

		%update edges
		[posedges dummy]=find(geom.Edges==i);
		geom.Edges(posedges,:)=[];
		posedges=find(geom.Edges>i);
		geom.Edges(posedges)=geom.Edges(posedges)-1;

		%update counter
		i=i-1;
	end
end
if num,
	disp(['WARNING: ' num2str(num) ' points outside the domain outline have been removed']);
end

%Check point spacing
if ~isnan(tol),
	disp('Checking point spacing...');
	i=0;
	while (i<size(geom.Vertices,1)),

		%vertex counter
		i=i+1;

		%Get coordinates
		x1=geom.Vertices(i,1);
		y1=geom.Vertices(i,2);

		j=i; %test edges located AFTER i only
		while (j<size(geom.Vertices,1)),

			%vertex counter
			j=j+1;

			%Get coordinates
			x2=geom.Vertices(j,1);
			y2=geom.Vertices(j,2);

			%Check whether the two vertices are too close
			if ((x2-x1)^2+(y2-y1)^2<tol^2)

				%Remove points from list of Vertices
				geom.Vertices(j,:)=[];

				%update edges
				posedges=find(ismember(geom.Edges,j));
				geom.Edges(posedges)=i;
				posedges=find(geom.Edges>j);
				geom.Edges(posedges)=geom.Edges(posedges)-1;

				%update counter
				j=j-1;

			end
		end
	end
end
%remove empty edges
geom.Edges(find(geom.Edges(:,1)==geom.Edges(:,2)),:)=[];
end % }}}

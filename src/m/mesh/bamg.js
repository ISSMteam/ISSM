function bamg(md){
    //BAMG - mesh generation
    //
    //   Available options (for more details see ISSM website http://issm.jpl.nasa.gov/):
    //
    //   - domain :            followed by an ARGUS file that prescribes the domain outline
    //   - holes :             followed by an ARGUS file that prescribes the holes
    //   - subdomains :        followed by an ARGUS file that prescribes the list of
    //                         subdomains (that need to be inside domain)
    //
    //   - hmin :              minimum edge length (default is 10^-100)
    //   - hmax :              maximum edge length (default is 10^100)
    //   - hVertices :         imposed edge length for each vertex (geometry or mesh)
    //   - hminVertices :      minimum edge length for each vertex (mesh)
    //   - hmaxVertices :      maximum edge length for each vertex (mesh)
    //
    //   - anisomax :          maximum ratio between the smallest and largest edges (default is 10^30)
    //   - coeff :             coefficient applied to the metric (2-> twice as many elements, default is 1)
    //   - cutoff :            scalar used to compute the metric when metric type 2 or 3 are applied
    //   - err :               error used to generate the metric from a field
    //   - errg :              geometric error (default is 0.1)
    //   - field :             field of the model that will be used to compute the metric
    //                         to apply several fields, use one column per field
    //   - gradation :         maximum ratio between two adjacent edges
    //   - Hessiantype :       0 -> use double L2 projection (default)
    //                         1 -> use Green formula
    //   - KeepVertices :      try to keep initial vertices when adaptation is done on an existing mesh (default 1)
    //   - NoBoundaryRefinment: do not refine boundary, only follow contour provided (default 0)
    //   - maxnbv :            maximum number of vertices used to allocate memory (default is 10^6)
    //   - maxsubdiv :         maximum subdivision of exisiting elements (default is 10)
    //   - metric :            matrix (numberofnodes x 3) used as a metric
    //   - Metrictype :        0 -> absolute error          c/(err coeff^2) * Abs(H)        (default)
    //                         1 -> relative error          c/(err coeff^2) * Abs(H)/max(s,cutoff*max(s))
    //                         2 -> rescaled absolute error c/(err coeff^2) * Abs(H)/(smax-smin)
    //   - nbjacoby :          correction used by Hessiantype=1 (default is 1)
    //   - nbsmooth :          number of metric smoothing procedure (default is 3)
    //   - omega :             relaxation parameter of the smoothing procedure (default is 1.8)
    //   - power :             power applied to the metric (default is 1)
    //   - splitcorners :      split triangles whuch have 3 vertices on the outline (default is 1)
    //   - verbose :           level of verbosity (default is 1)
    //
    //   - vertical :          is this a 2d vertical mesh (flowband, default is 0)
    //   - rifts :             followed by an ARGUS file that prescribes the rifts
    //   - toltip :            tolerance to move tip on an existing point of the domain outline
    //   - tracks :            followed by an ARGUS file that prescribes the tracks that the mesh will stick to
    //   - RequiredVertices :  mesh vertices that are required. [x,y,ref]; ref is optional
    //   - tol :               if the distance between 2 points of the domain outline is less than tol, they
    //                         will be merged
    //
    //   Examples:
    //      md=bamg(md,'domain','DomainOutline.exp','hmax',3000);
    //      md=bamg(md,'field',[md.inversion.vel_obs md.geometry.thickness],'hmax',20000,'hmin',1000);
    //      md=bamg(md,'metric',A,'hmin',1000,'hmax',20000,'gradation',3,'anisomax',1);

    //process options
    var args = Array.prototype.slice.call(arguments);
    var options = new pairoptions(args.slice(1,args.length));
    options.deleteduplicates(1);

    //initialize the structures required as input of Bamg
    var bamg_options = {}
    var bamg_geometry = new bamggeom();
    var bamg_mesh = new bamgmesh();

    var subdomain_ref = 1;
    var hole_ref = 1;
    // Bamg Geometry parameters {{{
    if (options.exist('domain')) {

        //Check that file exists
        var domainfile=options.getfieldvalue('domain');
        if ((typeof domainfile) === 'string') {
            console.log('bamg error message: file ' + domainfile + ' loading from file path not supported - domain must be file object (.shp or .exp)');
        } else if ((typeof domainfile) === 'object') {
            domain = domainfile;
        } else {
            console.log('"domain" type not supported yet - ' + (typeof domainfile));
        }

        var holes = [];
        if (options.exist('holes')) {
            var holesfile=options.getfieldvalue('holes');
            if ((typeof holesfile) === 'string') {
                console.log('bamg error message: file ' + holesfile + ' loading from file path not supported - holes must be file object (.shp or .exp)');
            } else if ((typeof holesfile) === 'object') {
                holes = holesfile;
            } else {
                console.log('"holes" type not supported yet - ' + (typeof holesfile));
            }
        }
        var subdomains = [];
        if (options.exist('subdomains')) {
            var subdomainsfile=options.getfieldvalue('subdomains');
            if ((typeof subdomainsfile) === 'string') {
                console.log('bamg error message: file ' + subdomainsfile + ' loading from file path not supported - subdomains must be file object (.shp or .exp)');
            } else if ((typeof subdomainsfile) === 'object') {
                subdomains = subdomainsfile;
            } else {
                console.log('"subdomains" type not supported yet - ' + (typeof subdomainsfile));
            }
        }

        //Build geometry 
        var count=0;
        for (var i=0; i < domain.length; i++) {

            //Check that the domain is closed
            if (domain[i].x[0] != domain[i].x[domain[i].x.length-1] || domain[i].y[0] != domain[i].y[domain[i].y.length-1]) {
                console.log('bamg error message: all contours provided in "domain" should be closed');
            }

            //TODO: Implement ContourToNodes
            //Checks that all holes are INSIDE the principle domain outline
            //if (i>1) {
            //    flags=ContourToNodes(domain[i].x,domain[i].y,domain[0],0);
            //    if (ArrayAny(ArrayFlip(flags))) {
            //        console.log('bamg error message: All holes should be strictly inside the principal domain');
            //    }
            //}
            console.log('bamg warning message: All holes should be strictly inside the principal domain. No checks are currently implemented.');

            //Check orientation
            var nods = domain[i].nods-1; //the domain are closed 1=end;
            var test = ArraySum(ArrayMultiply(ArraySubtract(domain[i].x.slice(1,nods+1), domain[i].x.slice(0,nods)), ArrayAdd(domain[i].y.slice(1,nods+1), domain[i].y.slice(0,nods))));
            if ((i==0 && test>0) || (i>0 && test<0)) {
                console.log('At least one contour was not correctly oriented and has been re-oriented');
		domain[i].x.reverse();
		domain[i].y.reverse();
            }

            //Add all points to bamg_geometry
            for (var j = 0; j < nods; j++) {
                bamg_geometry.Vertices.push([domain[i].x[j], domain[i].y[j], 1]);
            }
            var edges1 = ArrayRange(count + 1, count + nods);
            var edges2 = ArrayConcat(ArrayRange(count + 2, count + nods), [count + 1]);
            for (var j = 0; j < nods; j++) {
                bamg_geometry.Edges.push([edges1[j], edges2[j], 1]);
            }
            if (i > 1) {
                bamg_geometry.SubDomains.push([2, count + 1, 1, -subdomain_ref]);
                subdomain_ref = subdomain_ref + 1;
            } else {
                bamg_geometry.SubDomains.push([2, count + 1, 1, 0]);
            }

            //update counter
            count=count+nods;
        }
        for (var i=0; i < holes.length; i++) {

            //Check that the subdomains is closed
            if (holes[i].x[0] != holes[i].x[holes[i].x.length-1] || holes[i].y[0] != holes[i].y[holes[i].y.length-1]) {
                console.log('bamg error message: all contours provided in "domain" should be closed');
            }
            //Checks that all holes are INSIDE the principle domain outline
            //flags=ContourToNodes(holes[i].x,holes[i].y,domain[0],0);
            //if ArrayAny(ArrayFlip(flags)) { console.log('bamg error message: All holes should be strictly inside the principal domain'); }

            //TODO: Implement ContourToNodes
            //Checks that all holes are INSIDE the principle domain outline
            //if (i>1) {
            //    flags=ContourToNodes(domain[i].x,domain[i].y,domain[0],0);
            //    if (ArrayAny(ArrayFlip(flags))) {
            //        console.log('bamg error message: All holes should be strictly inside the principal domain');
            //    }
            //}
            console.log('bamg warning message: all holes should be strictly inside the principal domain. no checks are currently implemented.');

            //Check that hole is correctly oriented
            var nods = holes[i].nods-1; //the holes are closed 1=end;
            var test = ArraySum(ArrayMultiply(ArraySubtract(holes[i].x.slice(1,nods+1), holes[i].x.slice(0,nods)), ArrayAdd(holes[i].y.slice(1,nods+1), holes[i].y.slice(0,nods))));
            if ((i==0 && test>0) || (i>0 && test<0)) {
                console.log('At least one contour was not correctly oriented and has been re-oriented');
		holes[i].x.reverse();
		holes[i].y.reverse();
            }

            //Add all points to bamg_geometry
            for (var j = 0; j < nods; j++) {
                bamg_geometry.Vertices.push([holes[i].x[j], holes[i].y[j], 1]);
            }
            var edges1 = ArrayRange(count + 1, count + nods);
            var edges2 = ArrayConcat(ArrayRange(count + 2, count + nods), [count + 1]);
            for (var j = 0; j < nods; j++) {
                bamg_geometry.Edges.push([edges1[j], edges2[j], 1]);
            }
            bamg_geometry.SubDomains.push([2, count + 1, 1, -hole_ref]);
            hole_ref = hole_ref + 1;

            //update counter
            count=count+nods;
        }
        for (var i=0; i < subdomains.length; i++) {

            //Check that the subdomains is closed
            if (subdomains[i].x[0] != subdomains[i].x[subdomains[i].x.length-1] || subdomains[i].y[0] != subdomains[i].y[subdomains[i].y.length-1]) {
                console.log('bamg error message: all contours provided in "subdomains" should be closed');
            }

            //TODO: Implement ContourToNodes
            //Checks that all holes are INSIDE the principle domain outline
            //flags=ContourToNodes(subdomains[i].x,subdomains[i].y,domain[0],0);
            //if ArrayAny(ArrayFlip(flags)) {
            //    console.log('bamg error message: All holes should be strictly inside the principal domain');
            //}
            console.log('bamg warning message: all holes should be strictly inside the principal domain. no checks are currently implemented.');

            //Check that hole is correctly oriented
            var nods=subdomains[i].nods-1; //the subdomains are closed 1=end;
            var test = ArraySum(ArrayMultiply(ArraySubtract(subdomain[i].x.slice(1,nods+1), subdomain[i].x.slice(0,nods)), ArrayAdd(subdomain[i].y.slice(1,nods+1), subdomain[i].y.slice(0,nods))));
            if ((i==0 && test>0) || (i>0 && test<0)) {
                console.log('At least one contour was not correctly oriented and has been re-oriented');
		subdomains[i].x.reverse();
		subdomains[i].y.reverse();
            }

            //Add all points to bamg_geometry
            for (var j = 0; j < nods; j++) {
                bamg_geometry.Vertices.push([subdomains[i].x[j], subdomains[i].y[j], 1]);
            }
            var edges1 = ArrayRange(count + 1, count + nods);
            var edges2 = ArrayConcat(ArrayRange(count + 2, count + nods), [count + 1]);
            for (var j = 0; j < nods; j++) {
                bamg_geometry.Edges.push([edges1[j], edges2[j], 1]);
            }
            bamg_geometry.SubDomains.push([2, count + 1, 1, subdomain_ref]);
            subdomain_ref = subdomain_ref + 1;

            //update counter
            count=count+nods;
        }
        if (options.getfieldvalue('vertical',0)) {
            if (options.getfieldvalue('Markers',[]).length != bamg_geometry.Edges.length) {
                console.log('for 2d vertical mesh, "Markers" option is required, and should be of size ' + bamg_geometry.Edges.length);
            }
        }
        if (options.getfieldvalue('Markers',[]).length == bamg_geometry.Edges.length) {
            var markers = options.getfieldvalue('Markers');
            for (var i = 0; i < markers.length; i++) {
                bamg_geometry.Edges[i][2] = markers[i];
            }
        }
        /*

        //take care of rifts
        if options.exist('rifts') {

            //Check that file exists
            riftfile=options.getfieldvalue('rifts');
            [pathr,namer,extr]=fileparts(riftfile);
            if !exist(riftfile {'file')
                console.log(['bamg error message: file ' riftfile ' not found ']);
            } else if strcmp(extr,'.exp') {
                rift=expread(riftfile);
            } else if strcmp(extr,'.shp') {
                rift=shpread(riftfile);
            }
            //read rift file according to its extension: 
            [path,name,ext]=fileparts(riftfile);
            if strcmp(ext,'.exp') {
                rift=expread(riftfile);
            } else if strcmp(ext,'.shp') {
                rift=shpread(riftfile);
            } else {
                console.log(['bamg error message: file ' riftfile ' format not supported (.shp or .exp)']);
            }

            for i=1:length(rift) {

                //detect whether all points of the rift are inside the domain
                flags=ContourToNodes(rift[i].x,rift[i].y,domain[0],0);
                if (ArrayFlip(flags)) {
                    console.log('one rift has all its points outside of the domain outline'),

                } else if (ArrayAny(ArrayFlip(flags))) {
                    //We LOTS of work to do
                    console.log('Rift tip outside of or on the domain has been detected and is being processed...');

                    //check that only one point is outside (for now)
                    if (ArraySum(ArrayFlip(flags))!=1) {
                        console.log('bamg error message: only one point outside of the domain is supported yet');
                    }

                    //Move tip outside to the first position
                    if (flags[0]==0) {
                        //OK, first point is outside (do nothing),
                    } else if (flags[flags.length-1]==0) {
                        rift[i].x=flipud(rift[i].x);
                        rift[i].y=flipud(rift[i].y);
                    } else {
                        console.log('bamg error message: only a rift tip can be outside of the domain');
                    }

                    //Get cordinate of intersection point
                    x1=rift[i].x[0]; y1=rift[i].y[0];
                    x2=rift[i].x[1]; y2=rift[i].y[1];
                    for (var j=0; j < domain[0].x)-1 j++) {
                        if SegIntersect([x1 y1; x2 y2],[domain[0].x(j) domain[0].y(j); domain[0].x(j+1) domain[0].y(j+1)]) {

                            //Get position of the two nodes of the edge in domain
                            i1=j;
                            i2=j+1;

                            //rift is crossing edge [i1 i2] of the domain
                            //Get coordinate of intersection point (http://mathworld.wolfram.com/Line-LineIntersection.html)
                            x3=domain[0].x[i1]; y3=domain[0].y[i1];
                            x4=domain[0].x[i2]; y4=domain[0].y[i2];
                            x=det([det([x1 y1; x2 y2])  x1-x2;det([x3 y3; x4 y4])  x3-x4])/det([x1-x2 y1-y2;x3-x4 y3-y4]);
                            y=det([det([x1 y1; x2 y2])  y1-y2;det([x3 y3; x4 y4])  y3-y4])/det([x1-x2 y1-y2;x3-x4 y3-y4]);

                            segdis= sqrt((x4-x3)^2+(y4-y3)^2);
                            tipdis=[sqrt((x-x3)^2+(y-y3)^2)  sqrt((x-x4)^2+(y-y4)^2)];

                            if (min(tipdis)/segdis) < options.getfieldvalue('toltip',0) {
                                disp('moving tip-domain intersection point');

                                //Get position of the closer point
                                if tipdis[0]>tipdis[1] {
                                    pos=i2;
                                } else {
                                    pos=i1;
                                }

                                //This point is only in Vertices (number pos).
                                //OK, now we can add our own rift
                                nods=rift[i].nods-1;
                                bamg_geometry.Vertices=[bamg_geometry.Vertices; [rift[i].x(2:end) rift[i].y(2:end) ones(nods,1)]];
                                bamg_geometry.Edges=[bamg_geometry.Edges;...
                                    pos count+1  (1+i);...
                                    [transpose(count+1:count+nods-1) transpose(count+2:count+nods)  (1+i)*ones(nods-1,1)]];
                                count=count+nods;

                                break;

                            } else {
                                //Add intersection point to Vertices
                                bamg_geometry.Vertices=[bamg_geometry.Vertices; x y 1];
                                count=count+1;

                                //Decompose the crossing edge into 2 subedges
                                pos=find(bamg_geometry.Edges(:,1)==i1 & bamg_geometry.Edges(:,2)==i2);
                                if isempty(pos) console.log('bamg error message: a problem occurred...'); }
                                bamg_geometry.Edges=[bamg_geometry.Edges(1:pos-1,:);...
                                    bamg_geometry.Edges(pos,1) count                      bamg_geometry.Edges(pos,3);...
                                    count                      bamg_geometry.Edges(pos,2) bamg_geometry.Edges(pos,3);...
                                    bamg_geometry.Edges(pos+1:end,:)];

                                //OK, now we can add our own rift
                                nods=rift[i].nods-1;
                                bamg_geometry.Vertices=[bamg_geometry.Vertices; [rift[i].x(2:end) rift[i].y(2:end) ones(nods,1)]];
                                bamg_geometry.Edges=[bamg_geometry.Edges;...
                                    count  count+1  2 ;...
                                    [transpose(count+1:count+nods-1) transpose(count+2:count+nods)  (1+i)*ones(nods-1,1)]];
                                count=count+nods;

                                break;
                            }
                        }
                    }
                } else {
                    nods=rift[i].nods-1;
                    bamg_geometry.Vertices=[bamg_geometry.Vertices; [rift[i].x(:) rift[i].y(:) ones(nods+1,1)]];
                    bamg_geometry.Edges=[bamg_geometry.Edges; [transpose(count+1:count+nods) transpose(count+2:count+nods+1)  (1+i)*ones(nods,1)]];
                    count=count+nods+1;
                }
            }
        }

        //Deal with tracks
        if options.exist('tracks') {

            //read tracks
            track=options.getfieldvalue('tracks');
            if all(ischar(track)) {
                A=expread(track);
                track=[];
                for i=1:length(A), 
                    track=[track; [A[i].x A[i].y]];
                }
            } else {
                track=double(track); //for some reason, it is of class "single"
            }
            if(size(track,2)==2), track=[track 3.*ones(size(track,1),1)]; }

            //only keep those inside
            flags=ContourToNodes(track(:,1),track(:,2),domainfile,0);
            track=track(find(flags),:);

            //Add all points to bamg_geometry
            nods=size(track,1);
            bamg_geometry.Vertices=[bamg_geometry.Vertices; track];
            bamg_geometry.Edges=[bamg_geometry.Edges; [transpose(count+1:count+nods-1) transpose(count+2:count+nods)  3.*ones(nods-1,1)]];

            //update counter
            count=count+nods;
        }

        //Deal with vertices that need to be kept by mesher
        if (options.exist('RequiredVertices')) {

            //recover RequiredVertices
            requiredvertices = options.getfieldvalue('RequiredVertices'); //for some reason, it is of class "single"
            if (requiredvertices[0].length == 2) {
                 requiredvertices=[requiredvertices 4.*ones(size(requiredvertices,1),1)];
            }    

            //only keep those inside
            flags=ContourToNodes(requiredvertices(:,1),requiredvertices(:,2),domain[0],0);
            requiredvertices=requiredvertices(find(flags),:);

            //Add all points to bamg_geometry
            nods=size(requiredvertices,1);
            bamg_geometry.Vertices=[bamg_geometry.Vertices; requiredvertices];

            //update counter
            count=count+nods;

        }
        */

        //Deal with RequiredEdges
        if (options.getfieldvalue('NoBoundaryRefinment', 0) == 1) {
            bamg_geometry.RequiredEdges = ArrayTranspose(ArrayRange(1, bamg_geometry.Edges.length));
        }

        //process geom
        //bamg_geometry=processgeometry(bamg_geometry,options.getfieldvalue('tol',NaN),domain[0]);

    } else if ((typeof md.priv.bamg === 'object') && ('geometry' in md.priv.bamg)) {
        bamg_geometry = new bamggeom(md.priv.bamg.geometry); 
    } else {
        //do nothing...
    }
    //}}}
    // Bamg Mesh parameters {{{
    if (!options.exist('domain') && md.mesh.numberofvertices != 0 && md.mesh.elementtype() == 'Tria') {

        if ((typeof md.priv.bamg === 'object') && ('mesh' in md.priv.bamg)) {
            bamg_mesh = new bamgmesh(md.priv.bamg.mesh);
        } else {
            for (var i = 0; i < md.mesh.numberofvertices; i++) {
                bamg_mesh.Vertices.push([md.mesh.x[i], md.mesh.y.y[i], 1]);
            }
            for (var i = 0; i < md.mesh.numberofelements; i++) {
                bamg_mesh.Triangles.push([md.mesh.elements[i][0], md.mesh.elements[i][1], md.mesh.elements[i][2], 1]);
            }
        }

        if (typeof md.rifts.riftstruct === 'object') {
            console.log('bamg error message: rifts not supported yet. Do meshprocessrift AFTER bamg');
        }
    }
    //}}}
    // Bamg Options {{{
    bamg_options.Crack=options.getfieldvalue('Crack',0);
    bamg_options.anisomax=options.getfieldvalue('anisomax',Math.pow(10,30));
    bamg_options.coeff=options.getfieldvalue('coeff',1.);
    bamg_options.cutoff=options.getfieldvalue('cutoff',Math.pow(10,-5));
    bamg_options.err=options.getfieldvalue('err',0.01);
    bamg_options.errg=options.getfieldvalue('errg',0.1);
    bamg_options.field=options.getfieldvalue('field',[]);
    bamg_options.gradation=options.getfieldvalue('gradation',1.5);
    bamg_options.Hessiantype=options.getfieldvalue('Hessiantype',0);
    bamg_options.hmin=options.getfieldvalue('hmin',Math.pow(10,-100));
    bamg_options.hmax=options.getfieldvalue('hmax',Math.pow(10,100));
    bamg_options.hminVertices=options.getfieldvalue('hminVertices',[]);
    bamg_options.hmaxVertices=options.getfieldvalue('hmaxVertices',[]);
    bamg_options.hVertices=options.getfieldvalue('hVertices',[]);
    bamg_options.KeepVertices=options.getfieldvalue('KeepVertices',1);
    bamg_options.maxnbv=options.getfieldvalue('maxnbv',Math.pow(10,6));
    bamg_options.maxsubdiv=options.getfieldvalue('maxsubdiv',10.);
    bamg_options.metric=options.getfieldvalue('metric',[]);
    bamg_options.Metrictype=options.getfieldvalue('Metrictype',0);
    bamg_options.nbjacobi=options.getfieldvalue('nbjacobi',1);
    bamg_options.nbsmooth=options.getfieldvalue('nbsmooth',3);
    bamg_options.omega=options.getfieldvalue('omega',1.8);
    bamg_options.power=options.getfieldvalue('power',1.);
    bamg_options.splitcorners=options.getfieldvalue('splitcorners',1);
    bamg_options.verbose=options.getfieldvalue('verbose',1);
    //}}}

    //call Bamg
    console.log("calling BamgMesher");
    var return_array=BamgMesher(bamg_mesh,bamg_geometry,bamg_options);
    var bamgmesh_out=return_array[0];
    var bamggeom_out=return_array[1];

    if (options.getfieldvalue('vertical', 0) != 0) {
        md.mesh                     = new mesh2dvertical();
        md.mesh.x                   = ArrayCol(bamgmesh_out.Vertices, 0);
        md.mesh.y                   = ArrayCol(bamgmesh_out.Vertices, 1);
        md.mesh.elements            = ArrayCol(bamgmesh_out.Triangles, [0, 2]);
        md.mesh.edges               = bamgmesh_out.IssmEdges;
        md.mesh.segments            = ArrayCol(bamgmesh_out.IssmSegments, [0, 2]);
        md.mesh.segmentmarkers      = ArrayCol(bamgmesh_out.IssmSegments, 3);

        //Fill in rest of fields:
        md.mesh.numberofelements    = md.mesh.elements.length;
        md.mesh.numberofvertices    = md.mesh.x.length;
        md.mesh.numberofedges       = md.mesh.edges.length;
        for (var i = 0; i < md.mesh.segments.length; i++) {
            md.mesh.vertexonboundary[md.mesh.segments[i][0]] = 1;
            md.mesh.vertexonboundary[md.mesh.segments[i][1]] = 1;
        }
    } else if (options.getfieldvalue('3dsurface', 0) != 0) {
        md.mesh                     = new mesh3dsurface();
        md.mesh.x                   = ArrayCol(bamgmesh_out.Vertices, 0);
        md.mesh.y                   = ArrayCol(bamgmesh_out.Vertices, 1);
        md.mesh.z                   = NewArrayFill(md.mesh.x.length, 0);
        md.mesh.elements            = ArrayCol(bamgmesh_out.Triangles, [0, 2]);
        md.mesh.edges               = bamgmesh_out.IssmEdges;
        md.mesh.segments            = ArrayCol(bamgmesh_out.IssmSegments, [0, 2]);
        md.mesh.segmentmarkers      = ArrayCol(bamgmesh_out.IssmSegments, 3);

        //Fill in rest of fields:
        md.mesh.numberofelements    = md.mesh.elements.length;
        md.mesh.numberofvertices    = md.mesh.x.length;
        md.mesh.numberofedges       = md.mesh.edges.length;
        for (var i = 0; i < md.mesh.segments.length; i++) {
            md.mesh.vertexonboundary[md.mesh.segments[i][0]] = 1;
            md.mesh.vertexonboundary[md.mesh.segments[i][1]] = 1;
        }
    } else { 
        md.mesh                     = new mesh2d();
        md.mesh.x                   = ArrayCol(bamgmesh_out.Vertices, 0);
        md.mesh.y                   = ArrayCol(bamgmesh_out.Vertices, 1);
        md.mesh.elements            = ArrayCol(bamgmesh_out.Triangles, [0, 2]);
        md.mesh.edges               = bamgmesh_out.IssmEdges;
        md.mesh.segments            = ArrayCol(bamgmesh_out.IssmSegments, [0, 2]);
        md.mesh.segmentmarkers      = ArrayCol(bamgmesh_out.IssmSegments, 3);

        //Fill in rest of fields:
        md.mesh.numberofelements    = md.mesh.elements.length;
        md.mesh.numberofvertices    = md.mesh.x.length;
        md.mesh.numberofedges       = md.mesh.edges.length;
        md.mesh.vertexonboundary    = NewArrayFill(md.mesh.numberofvertices, 0);
        for (var i = 0; i < md.mesh.segments.length; i++) {
            md.mesh.vertexonboundary[md.mesh.segments[i][0]] = 1;
            md.mesh.vertexonboundary[md.mesh.segments[i][1]] = 1;
        }
    }

    //Bamg private fields
    md.priv.bamg                 = [];
    md.priv.bamg.mesh            = new bamgmesh(bamgmesh_out);
    md.priv.bamg.geometry        = new bamggeom(bamggeom_out);
    md.mesh.elementconnectivity  = md.priv.bamg.mesh.ElementConnectivity;
    for (var i = 0; i < md.mesh.elementconnectivity.length; i++) {
        if (isNaN(md.mesh.elementconnectivity[i])) {
            md.mesh.elementconnectivity[i] = 0;
        }
    }

    //Check for orphan
    for (var i = 0; i < md.mesh.numberofelements; i++) {
        for (var j = 0; j < 3; j++) {
            if (md.mesh.elements[i][j] > md.mesh.numberofvertices) {
                console.log('Output mesh has orphans. Check your Domain and/or RequiredVertices');
                break;
            }
        }
    }
} 

function processgeometry(geom,tol,outline){ // {{{

//    //Deal with edges
//    disp('Checking Edge crossing...');
//    i=0;
//    while (i<size(geom.Edges,1)),
//
//        //edge counter
//        i=i+1;
//
//        //Get coordinates
//        x1=geom.Vertices(geom.Edges(i,1),1);
//        y1=geom.Vertices(geom.Edges(i,1),2);
//        x2=geom.Vertices(geom.Edges(i,2),1);
//        y2=geom.Vertices(geom.Edges(i,2),2);
//        color1=geom.Edges(i,3);
//
//        j=i; //test edges located AFTER i only
//        while (j<size(geom.Edges,1)),
//
//            //edge counter
//            j=j+1;
//
//            //Skip if the two edges already have a vertex in common
//            if ArrayAny(ismember(geom.Edges(i,1:2),geom.Edges(j,1:2))),
//                continue
//            }
//
//            //Get coordinates
//            x3=geom.Vertices(geom.Edges(j,1),1);
//            y3=geom.Vertices(geom.Edges(j,1),2);
//            x4=geom.Vertices(geom.Edges(j,2),1);
//            y4=geom.Vertices(geom.Edges(j,2),2);
//            color2=geom.Edges(j,3);
//
//            //Check if the two edges are crossing one another
//            if SegIntersect([x1 y1; x2 y2],[x3 y3; x4 y4]),
//
//                //Get coordinate of intersection point (http://mathworld.wolfram.com/Line-LineIntersection.html)
//                x=det([det([x1 y1; x2 y2])  x1-x2;det([x3 y3; x4 y4])  x3-x4])/det([x1-x2 y1-y2;x3-x4 y3-y4]);
//                y=det([det([x1 y1; x2 y2])  y1-y2;det([x3 y3; x4 y4])  y3-y4])/det([x1-x2 y1-y2;x3-x4 y3-y4]);
//
//                //Add vertex to the list of vertices
//                geom.Vertices(end+1,:)=[x y min(color1,color2)];
//                id=size(geom.Vertices,1);
//
//                //Update edges i and j
//                edgei=geom.Edges(i,:);
//                edgej=geom.Edges(j,:);
//                geom.Edges(i,:)    =[edgei[0] id       edgei(3)];
//                geom.Edges(end+1,:)=[id       edgei[1] edgei(3)];
//                geom.Edges(j,:)    =[edgej[0] id       edgej(3)];
//                geom.Edges(end+1,:)=[id       edgej[1] edgej(3)];
//
//                //update current edge second tip
//                x2=x; y2=y;
//            }
//        }
//
//    }
//
//    //Check point outside
//    disp('Checking for points outside the domain...');
//    i=0;
//    num=0;
//    while (i<size(geom.Vertices,1)),
//
//        //vertex counter
//        i=i+1;
//
//        //Get coordinates
//        x=geom.Vertices(i,1);
//        y=geom.Vertices(i,2);
//        color=geom.Vertices(i,3);
//
//        //Check that the point is inside the domain
//        if (color!=1 & !ContourToNodes(x,y,outline[0],1)),
//
//            //Remove points from list of Vertices
//            num=num+1;
//            geom.Vertices(i,:)=[];
//
//            //update edges
//            [posedges dummy]=find(geom.Edges==i);
//            geom.Edges(posedges,:)=[];
//            posedges=find(geom.Edges>i);
//            geom.Edges(posedges)=geom.Edges(posedges)-1;
//
//            //update counter
//            i=i-1;
//        }
//    }
//    if num,
//        disp(['WARNING: ' num2str(num) ' points outside the domain outline have been removed']);
//    }
//
//    //Check point spacing
//    if !isnan(tol),
//        disp('Checking point spacing...');
//        i=0;
//        while (i<size(geom.Vertices,1)),
//
//            //vertex counter
//            i=i+1;
//
//            //Get coordinates
//            x1=geom.Vertices(i,1);
//            y1=geom.Vertices(i,2);
//
//            j=i; //test edges located AFTER i only
//            while (j<size(geom.Vertices,1)),
//
//                //vertex counter
//                j=j+1;
//
//                //Get coordinates
//                x2=geom.Vertices(j,1);
//                y2=geom.Vertices(j,2);
//
//                //Check whether the two vertices are too close
//                if ((x2-x1)^2+(y2-y1)^2<tol^2)
//
//                    //Remove points from list of Vertices
//                    geom.Vertices(j,:)=[];
//
//                    //update edges
//                    posedges=find(ismember(geom.Edges,j));
//                    geom.Edges(posedges)=i;
//                    posedges=find(geom.Edges>j);
//                    geom.Edges(posedges)=geom.Edges(posedges)-1;
//
//                    //update counter
//                    j=j-1;
//
//                }
//            }
//        }
//    }
//    //remove empty edges
//    geom.Edges(find(geom.Edges(:,1)==geom.Edges(:,2)),:)=[];
} // }}}

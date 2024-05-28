function doublesToHeap(array) {
    var doubleArray = new Float64Array(array);
    var numBytes = doubleArray.length * doubleArray.BYTES_PER_ELEMENT;
    var doubleArrayPtr = Module._malloc(numBytes);
    var doubleArrayHeap = new Uint8Array(Module.HEAPU8.buffer, doubleArrayPtr, numBytes);
    doubleArrayHeap.set(new Uint8Array(doubleArray.buffer));
    return doubleArrayHeap.byteOffset;
}
function intsToHeap(array) {
    var intArray = new Int32Array(array);
    var numBytes = intArray.length * intArray.BYTES_PER_ELEMENT;
    var intArrayPtr = Module._malloc(numBytes);
    var intArrayHeap = new Uint8Array(Module.HEAPU8.buffer, intArrayPtr, numBytes);
    intArrayHeap.set(new Uint8Array(intArray.buffer));
    return intArrayHeap.byteOffset;
}
function heapToDoubles(pptr, size) {
    var ptr = Module.getValue(pptr,'i32');
    var array = Module.HEAPF64.slice(ptr / 8, ptr / 8 + size[0] * size[1]);
    return ListToMatrix(array, size[1]);
}
function heapToInts(pptr, nods) {
    var ptr = Module.getValue(pptr,'i32');
    return Module.HEAPU32.slice(ptr / 4, ptr / 4 + nods);
}
function size2d(array) {
    // If array = Array(0), size2d(array) == [0, 1]
    // If array = [1, 2, 3] size2d(array) == [3, 1]
    // If array = [[1, 2, 3], [2, 3, 4]], size2d(array) == [2, 3]
    var size0 = array.length;
    var size1 = 1;
    if (array[0] instanceof Array) {
        size1 = array[0].length;
    }
    //if (array[0] != undefined && array[0].length != undefined) {
    //    size1 = array[0].length;
    //}
    return [size0, size1];
}
function BamgMesher(bamgmesh_in, bamggeom_in, bamgopts) {
/*
       usage: var array = Triangle(domain,rifts,area);
          where: array is made of [index,x,y,segments,segmentmarkers]
          and index,x,y defines a triangulation, segments is an array made
          of exterior segments to the mesh domain outline, segmentmarkers is an array
          flagging each segment, domain a js array defining the domain outline  (sames for
          rifts) and area is the maximum area desired for any element of the resulting mesh.

          Ok, for now, we are not dealing with rifts. Also, the domain is made of only one
          profile, this to avoid passing a double** pointer to js.
*/

    //Dynamic allocations: {{{
    //Retrieve domain arrays, and allocate on Module heap:
    //For each property, calculate the size and fill with 0 if the 2nd dimension is undefined, then use the int size array to init the double array.
    //input
    var pVerticesSize_mesh_in               = intsToHeap(size2d(bamgmesh_in.Vertices));
    var pVertices_mesh_in                   = doublesToHeap(Array.prototype.concat.apply([], bamgmesh_in.Vertices));
    var pEdgesSize_mesh_in                  = intsToHeap(size2d(bamgmesh_in.Edges));
    var pEdges_mesh_in                      = doublesToHeap(Array.prototype.concat.apply([], bamgmesh_in.Edges));
    var pTrianglesSize_mesh_in              = intsToHeap(size2d(bamgmesh_in.Triangles));
    var pTriangles_mesh_in                  = doublesToHeap(Array.prototype.concat.apply([], bamgmesh_in.Triangles));
    var pCrackedEdgesSize_mesh_in           = intsToHeap(size2d(bamgmesh_in.CrackedEdges));
    var pCrackedEdges_mesh_in               = doublesToHeap(Array.prototype.concat.apply([], bamgmesh_in.CrackedEdges));
    var pVerticesOnGeomEdgeSize_mesh_in     = intsToHeap(size2d(bamgmesh_in.VerticesOnGeomEdge));
    var pVerticesOnGeomEdge_mesh_in         = doublesToHeap(Array.prototype.concat.apply([], bamgmesh_in.VerticesOnGeomEdge));
    var pVerticesOnGeomVertexSize_mesh_in   = intsToHeap(size2d(bamgmesh_in.VerticesOnGeomVertex));
    var pVerticesOnGeomVertex_mesh_in       = doublesToHeap(Array.prototype.concat.apply([], bamgmesh_in.VerticesOnGeomVertex));
    var pEdgesOnGeomEdgeSize_mesh_in        = intsToHeap(size2d(bamgmesh_in.EdgesOnGeomEdge));
    var pEdgesOnGeomEdge_mesh_in            = doublesToHeap(Array.prototype.concat.apply([], bamgmesh_in.EdgesOnGeomEdge));
    var pIssmSegmentsSize_mesh_in           = intsToHeap(size2d(bamgmesh_in.IssmSegments));
    var pIssmSegments_mesh_in               = doublesToHeap(Array.prototype.concat.apply([], bamgmesh_in.IssmSegments));

    var pVerticesSize_geom_in               = intsToHeap(size2d(bamggeom_in.Vertices));
    var pVertices_geom_in                   = doublesToHeap(Array.prototype.concat.apply([], bamggeom_in.Vertices));
    var pEdgesSize_geom_in                  = intsToHeap(size2d(bamggeom_in.Edges));
    var pEdges_geom_in                      = doublesToHeap(Array.prototype.concat.apply([], bamggeom_in.Edges));
    var pCornersSize_geom_in                = intsToHeap(size2d(bamggeom_in.Corners));
    var pCorners_geom_in                    = doublesToHeap(Array.prototype.concat.apply([], bamggeom_in.Corners));
    var pRequiredVerticesSize_geom_in       = intsToHeap(size2d(bamggeom_in.RequiredVertices));
    var pRequiredVertices_geom_in           = doublesToHeap(Array.prototype.concat.apply([], bamggeom_in.RequiredVertices));
    var pRequiredEdgesSize_geom_in          = intsToHeap(size2d(bamggeom_in.RequiredEdges));
    var pRequiredEdges_geom_in              = doublesToHeap(Array.prototype.concat.apply([], bamggeom_in.RequiredEdges));
    var pCrackedEdgesSize_geom_in           = intsToHeap(size2d(bamggeom_in.CrackedEdges));
    var pCrackedEdges_geom_in               = doublesToHeap(Array.prototype.concat.apply([], bamggeom_in.CrackedEdges));
    var pSubDomainsSize_geom_in             = intsToHeap(size2d(bamggeom_in.SubDomains));
    var pSubDomains_geom_in                 = doublesToHeap(Array.prototype.concat.apply([], bamggeom_in.SubDomains));

    var anisomax                            = bamgopts.anisomax;
    var cutoff                              = bamgopts.coeff;
    var coeff                               = bamgopts.cutoff;
    var errg                                = bamgopts.errg;
    var gradation                           = bamgopts.gradation;
    var Hessiantype                         = bamgopts.Hessiantype;
    var maxnbv                              = bamgopts.maxnbv;
    var maxsubdiv                           = bamgopts.maxsubdiv;
    var Metrictype                          = bamgopts.Metrictype;
    var nbjacobi                            = bamgopts.nbjacobi;
    var nbsmooth                            = bamgopts.nbsmooth;
    var omega                               = bamgopts.omega;
    var power                               = bamgopts.power;
    var verbose                             = bamgopts.verbose;
    var Crack                               = bamgopts.Crack;
    var KeepVertices                        = bamgopts.KeepVertices;
    var splitcorners                        = bamgopts.splitcorners;
    var hmin                                = bamgopts.hmin;
    var hmax                                = bamgopts.hmax;
    var phminVerticesSize                   = intsToHeap(size2d(bamgopts.hminVertices));
    var phminVertices                       = doublesToHeap(Array.prototype.concat.apply([], bamgopts.hminVertices));
    var phmaxVerticesSize                   = intsToHeap(size2d(bamgopts.hmaxVertices));
    var phmaxVertices                       = doublesToHeap(Array.prototype.concat.apply([], bamgopts.hmaxVertices));
    var hVerticesLength                     = bamgopts.hVerticesLength;
    var phVertices                          = doublesToHeap(Array.prototype.concat.apply([], bamgopts.hVertices));
    var pmetricSize                         = intsToHeap(size2d(bamgopts.metric));
    var pmetric                             = doublesToHeap(Array.prototype.concat.apply([], bamgopts.metric));
    var pfieldSize                          = intsToHeap(size2d(bamgopts.field));
    var pfield                              = doublesToHeap(Array.prototype.concat.apply([], bamgopts.field));
    var perrSize                            = intsToHeap(size2d([[bamgopts.err]]));
    var perr                                = doublesToHeap(Array.prototype.concat.apply([], [[bamgopts.err]]));

    //output
    var pVerticesSize_mesh_out                  = Module._malloc(4);
    var pVertices_mesh_out                      = Module._malloc(4);
    var pEdgesSize_mesh_out                     = Module._malloc(4);
    var pEdges_mesh_out                         = Module._malloc(4);
    var pTrianglesSize_mesh_out                 = Module._malloc(4);
    var pTriangles_mesh_out                     = Module._malloc(4);
    var pIssmEdgesSize_mesh_out                 = Module._malloc(4);
    var pIssmEdges_mesh_out                     = Module._malloc(4);
    var pIssmSegmentsSize_mesh_out              = Module._malloc(4);
    var pIssmSegments_mesh_out                  = Module._malloc(4);
    var pVerticesOnGeomVertexSize_mesh_out      = Module._malloc(4);
    var pVerticesOnGeomVertex_mesh_out          = Module._malloc(4);
    var pVerticesOnGeomEdgeSize_mesh_out        = Module._malloc(4);
    var pVerticesOnGeomEdge_mesh_out            = Module._malloc(4);
    var pEdgesOnGeomEdgeSize_mesh_out           = Module._malloc(4);
    var pEdgesOnGeomEdge_mesh_out               = Module._malloc(4);
    var pSubDomainsSize_mesh_out                = Module._malloc(4);
    var pSubDomains_mesh_out                    = Module._malloc(4);
    var pSubDomainsFromGeomSize_mesh_out        = Module._malloc(4);
    var pSubDomainsFromGeom_mesh_out            = Module._malloc(4);
    var pElementConnectivitySize_mesh_out       = Module._malloc(4);
    var pElementConnectivity_mesh_out           = Module._malloc(4);
    var pNodalConnectivitySize_mesh_out         = Module._malloc(4);
    var pNodalConnectivity_mesh_out             = Module._malloc(4);
    var pNodalElementConnectivitySize_mesh_out  = Module._malloc(4);
    var pNodalElementConnectivity_mesh_out      = Module._malloc(4);
    var pCrackedVerticesSize_mesh_out           = Module._malloc(4);
    var pCrackedVertices_mesh_out               = Module._malloc(4);
    var pCrackedEdgesSize_mesh_out              = Module._malloc(4);
    var pCrackedEdges_mesh_out                  = Module._malloc(4);
    var pPreviousNumberingSize_mesh_out         = Module._malloc(4);
    var pPreviousNumbering_mesh_out             = Module._malloc(4);

    var pVerticesSize_geom_out                  = Module._malloc(4);
    var pVertices_geom_out                      = Module._malloc(4);
    var pEdgesSize_geom_out                     = Module._malloc(4);
    var pEdges_geom_out                         = Module._malloc(4);
    var pCornersSize_geom_out                   = Module._malloc(4);
    var pCorners_geom_out                       = Module._malloc(4);
    var pRequiredVerticesSize_geom_out          = Module._malloc(4);
    var pRequiredVertices_geom_out              = Module._malloc(4);
    var pRequiredEdgesSize_geom_out             = Module._malloc(4);
    var pRequiredEdges_geom_out                 = Module._malloc(4);
    var pCrackedEdgesSize_geom_out              = Module._malloc(4);
    var pCrackedEdges_geom_out                  = Module._malloc(4);
    var pSubDomainsSize_geom_out                = Module._malloc(4);
    var pSubDomains_geom_out                    = Module._malloc(4);
    //}}}

    //Declare BamgMesher module:
    BamgMesherModule = Module.cwrap('BamgMesherModule', 'number',[
        'number', 'number', 'number', 'number', 'number', 'number', 'number', 'number', 'number', 'number', 'number', 'number', 'number', 'number', 'number', 'number', 'number', 'number', 'number', 'number', 'number', 'number', 'number', 'number', 'number', 'number', 'number', 'number', 'number', 'number', 'number', 'number',
        'number', 'number', 'number', 'number', 'number', 'number', 'number', 'number', 'number', 'number', 'number', 'number', 'number', 'number',
        'number', 'number', 'number', 'number', 'number', 'number', 'number', 'number', 'number', 'number', 'number', 'number', 'number', 'number', 'number', 'number',
        'number', 'number', 'number', 'number', 'number', 'number', 'number', 'number', 'number', 'number', 'number', 'number', 'number', 'number',
        'number', 'number', 'number', 'number', 'number', 'number', 'number', 'number', 'number', 'number', 'number', 'number', 'number', 'number', 'number', 'number', 'number', 'number', 'number', 'number', 'number', 'number', 'number', 'number', 'number', 'number', 'number', 'number', 'number', 'number', 'number']);

    //Call BamgMesher module:

    BamgMesherModule(
        pVerticesSize_mesh_out, pVertices_mesh_out, pEdgesSize_mesh_out, pEdges_mesh_out, pTrianglesSize_mesh_out, pTriangles_mesh_out, pIssmEdgesSize_mesh_out, pIssmEdges_mesh_out, pIssmSegmentsSize_mesh_out, pIssmSegments_mesh_out, pVerticesOnGeomVertexSize_mesh_out, pVerticesOnGeomVertex_mesh_out, pVerticesOnGeomEdgeSize_mesh_out, pVerticesOnGeomEdge_mesh_out, pEdgesOnGeomEdgeSize_mesh_out, pEdgesOnGeomEdge_mesh_out, pSubDomainsSize_mesh_out, pSubDomains_mesh_out, pSubDomainsFromGeomSize_mesh_out, pSubDomainsFromGeom_mesh_out, pElementConnectivitySize_mesh_out, pElementConnectivity_mesh_out, pNodalConnectivitySize_mesh_out, pNodalConnectivity_mesh_out, pNodalElementConnectivitySize_mesh_out, pNodalElementConnectivity_mesh_out, pCrackedVerticesSize_mesh_out, pCrackedVertices_mesh_out, pCrackedEdgesSize_mesh_out, pCrackedEdges_mesh_out, pPreviousNumberingSize_mesh_out, pPreviousNumbering_mesh_out,
        pVerticesSize_geom_out, pVertices_geom_out, pEdgesSize_geom_out, pEdges_geom_out, pCornersSize_geom_out, pCorners_geom_out, pRequiredVerticesSize_geom_out, pRequiredVertices_geom_out, pRequiredEdgesSize_geom_out, pRequiredEdges_geom_out, pCrackedEdgesSize_geom_out, pCrackedEdges_geom_out, pSubDomainsSize_geom_out, pSubDomains_geom_out,
        pVerticesSize_mesh_in, pVertices_mesh_in, pEdgesSize_mesh_in, pEdges_mesh_in, pTrianglesSize_mesh_in, pTriangles_mesh_in, pCrackedEdgesSize_mesh_in, pCrackedEdges_mesh_in, pVerticesOnGeomEdgeSize_mesh_in, pVerticesOnGeomEdge_mesh_in, pVerticesOnGeomVertexSize_mesh_in, pVerticesOnGeomVertex_mesh_in, pEdgesOnGeomEdgeSize_mesh_in, pEdgesOnGeomEdge_mesh_in, pIssmSegmentsSize_mesh_in, pIssmSegments_mesh_in,
        pVerticesSize_geom_in, pVertices_geom_in, pEdgesSize_geom_in, pEdges_geom_in, pCornersSize_geom_in, pCorners_geom_in, pRequiredVerticesSize_geom_in, pRequiredVertices_geom_in, pRequiredEdgesSize_geom_in, pRequiredEdges_geom_in, pCrackedEdgesSize_geom_in, pCrackedEdges_geom_in, pSubDomainsSize_geom_in, pSubDomains_geom_in,
        anisomax, cutoff, coeff, errg, gradation, Hessiantype, maxnbv, maxsubdiv, Metrictype, nbjacobi, nbsmooth, omega, power, verbose, Crack, KeepVertices, splitcorners, hmin, hmax, phminVerticesSize, phminVertices, phmaxVerticesSize, phmaxVertices, hVerticesLength, phVertices, pmetricSize, pmetric, pfieldSize, pfield, perrSize, perr);

    /*Dynamic copying from heap: {{{*/
    //recover mesh:
    var bamgmeshout = new bamgmesh();
    var bamggeomout = new bamggeom();

    bamgmeshout.VerticesSize                   = heapToInts(pVerticesSize_mesh_out, 2);
    bamgmeshout.Vertices                       = heapToDoubles(pVertices_mesh_out, bamgmeshout.VerticesSize);
    bamgmeshout.EdgesSize                      = heapToInts(pEdgesSize_mesh_out, 2);
    bamgmeshout.Edges                          = heapToDoubles(pEdges_mesh_out, bamgmeshout.EdgesSize);
    bamgmeshout.TrianglesSize                  = heapToInts(pTrianglesSize_mesh_out, 2);
    bamgmeshout.Triangles                      = heapToDoubles(pTriangles_mesh_out, bamgmeshout.TrianglesSize);
    bamgmeshout.IssmEdgesSize                  = heapToInts(pIssmEdgesSize_mesh_out, 2);
    bamgmeshout.IssmEdges                      = heapToDoubles(pIssmEdges_mesh_out, bamgmeshout.IssmEdgesSize);
    bamgmeshout.IssmSegmentsSize               = heapToInts(pIssmSegmentsSize_mesh_out, 2);
    bamgmeshout.IssmSegments                   = heapToDoubles(pIssmSegments_mesh_out, bamgmeshout.IssmSegmentsSize);
    bamgmeshout.VerticesOnGeomVertexSize       = heapToInts(pVerticesOnGeomVertexSize_mesh_out, 2);
    bamgmeshout.VerticesOnGeomVertex           = heapToDoubles(pVerticesOnGeomVertex_mesh_out, bamgmeshout.VerticesOnGeomVertexSize);
    bamgmeshout.VerticesOnGeomEdgeSize         = heapToInts(pVerticesOnGeomEdgeSize_mesh_out, 2);
    bamgmeshout.VerticesOnGeomEdge             = heapToDoubles(pVerticesOnGeomEdge_mesh_out, bamgmeshout.VerticesOnGeomEdgeSize);
    bamgmeshout.EdgesOnGeomEdgeSize            = heapToInts(pEdgesOnGeomEdgeSize_mesh_out, 2);
    bamgmeshout.EdgesOnGeomEdge                = heapToDoubles(pEdgesOnGeomEdge_mesh_out, bamgmeshout.EdgesOnGeomEdgeSize);
    bamgmeshout.SubDomainsSize                 = heapToInts(pSubDomainsSize_mesh_out, 2);
    bamgmeshout.SubDomains                     = heapToDoubles(pSubDomains_mesh_out, bamgmeshout.SubDomainsSize);
    bamgmeshout.SubDomainsFromGeomSize         = heapToInts(pSubDomainsFromGeomSize_mesh_out, 2);
    bamgmeshout.SubDomainsFromGeom             = heapToDoubles(pSubDomainsFromGeom_mesh_out, bamgmeshout.SubDomainsFromGeomSize);
    bamgmeshout.ElementConnectivitySize        = heapToInts(pElementConnectivitySize_mesh_out, 2);
    bamgmeshout.ElementConnectivity            = heapToDoubles(pElementConnectivity_mesh_out, bamgmeshout.ElementConnectivitySize);
    bamgmeshout.NodalConnectivitySize          = heapToInts(pNodalConnectivitySize_mesh_out, 2);
    bamgmeshout.NodalConnectivity              = heapToDoubles(pNodalConnectivity_mesh_out, bamgmeshout.NodalConnectivitySize);
    bamgmeshout.NodalElementConnectivitySize   = heapToInts(pNodalElementConnectivitySize_mesh_out, 2);
    bamgmeshout.NodalElementConnectivity       = heapToDoubles(pNodalElementConnectivity_mesh_out, bamgmeshout.NodalElementConnectivitySize);
    bamgmeshout.CrackedVerticesSize            = heapToInts(pCrackedVerticesSize_mesh_out, 2);
    bamgmeshout.CrackedVertices                = heapToDoubles(pCrackedVertices_mesh_out, bamgmeshout.CrackedVerticesSize);
    bamgmeshout.CrackedEdgesSize               = heapToInts(pCrackedEdgesSize_mesh_out, 2);
    bamgmeshout.CrackedEdges                   = heapToDoubles(pCrackedEdges_mesh_out, bamgmeshout.CrackedEdgesSize);
    bamgmeshout.PreviousNumberingSize          = heapToInts(pPreviousNumberingSize_mesh_out, 2);
    bamgmeshout.PreviousNumbering              = heapToDoubles(pPreviousNumbering_mesh_out, bamgmeshout.PreviousNumberingSize);

    bamggeomout.VerticesSize                   = heapToInts(pVerticesSize_geom_out, 2);
    bamggeomout.Vertices                       = heapToDoubles(pVertices_geom_out, bamggeomout.VerticesSize);
    bamggeomout.EdgesSize                      = heapToInts(pEdgesSize_geom_out, 2);
    bamggeomout.Edges                          = heapToDoubles(pEdges_geom_out, bamggeomout.EdgesSize);
    bamggeomout.CornersSize                    = heapToInts(pCornersSize_geom_out, 2);
    bamggeomout.Corners                        = heapToDoubles(pCorners_geom_out, bamggeomout.CornersSize);
    bamggeomout.RequiredVerticesSize           = heapToInts(pRequiredVerticesSize_geom_out, 2);
    bamggeomout.RequiredVertices               = heapToDoubles(pRequiredVertices_geom_out, bamggeomout.RequiredVerticesSize);
    bamggeomout.RequiredEdgesSize              = heapToInts(pRequiredEdgesSize_geom_out, 2);
    bamggeomout.RequiredEdges                  = heapToDoubles(pRequiredEdges_geom_out, bamggeomout.RequiredEdgesSize);
    bamggeomout.CrackedEdgesSize               = heapToInts(pCrackedEdgesSize_geom_out, 2);
    bamggeomout.CrackedEdges                   = heapToDoubles(pCrackedEdges_geom_out, bamggeomout.CrackedEdgesSize);
    bamggeomout.SubDomainsSize                 = heapToInts(pSubDomainsSize_geom_out, 2);
    bamggeomout.SubDomains                     = heapToDoubles(pSubDomains_geom_out, bamggeomout.SubDomainsSize);
    /*}}}*/

    var return_array=[bamgmeshout, bamggeomout];

    /*Free resources: */
    Module._free(pVerticesSize_mesh_out);
    Module._free(pVertices_mesh_out);
    Module._free(pEdgesSize_mesh_out);
    Module._free(pEdges_mesh_out);
    Module._free(pTrianglesSize_mesh_out);
    Module._free(pTriangles_mesh_out);
    Module._free(pIssmEdgesSize_mesh_out);
    Module._free(pIssmEdges_mesh_out);
    Module._free(pIssmSegmentsSize_mesh_out);
    Module._free(pIssmSegments_mesh_out);
    Module._free(pVerticesOnGeomVertexSize_mesh_out);
    Module._free(pVerticesOnGeomVertex_mesh_out);
    Module._free(pVerticesOnGeomEdgeSize_mesh_out);
    Module._free(pVerticesOnGeomEdge_mesh_out);
    Module._free(pEdgesOnGeomEdgeSize_mesh_out);
    Module._free(pEdgesOnGeomEdge_mesh_out);
    Module._free(pSubDomainsSize_mesh_out);
    Module._free(pSubDomains_mesh_out);
    Module._free(pSubDomainsFromGeomSize_mesh_out);
    Module._free(pSubDomainsFromGeom_mesh_out);
    Module._free(pElementConnectivitySize_mesh_out);
    Module._free(pElementConnectivity_mesh_out);
    Module._free(pNodalConnectivitySize_mesh_out);
    Module._free(pNodalConnectivity_mesh_out);
    Module._free(pNodalElementConnectivitySize_mesh_out);
    Module._free(pNodalElementConnectivity_mesh_out);
    Module._free(pCrackedVerticesSize_mesh_out);
    Module._free(pCrackedVertices_mesh_out);
    Module._free(pCrackedEdgesSize_mesh_out);
    Module._free(pCrackedEdges_mesh_out);
    Module._free(pPreviousNumberingSize_mesh_out);
    Module._free(pPreviousNumbering_mesh_out);

    Module._free(pVerticesSize_geom_out);
    Module._free(pVertices_geom_out);
    Module._free(pEdgesSize_geom_out);
    Module._free(pEdges_geom_out);
    Module._free(pCornersSize_geom_out);
    Module._free(pCorners_geom_out);
    Module._free(pRequiredVerticesSize_geom_out);
    Module._free(pRequiredVertices_geom_out);
    Module._free(pRequiredEdgesSize_geom_out);
    Module._free(pRequiredEdges_geom_out);
    Module._free(pCrackedEdgesSize_geom_out);
    Module._free(pCrackedEdges_geom_out);
    Module._free(pSubDomainsSize_geom_out);
    Module._free(pSubDomains_geom_out);

    Module._free(pVerticesSize_mesh_in);
    Module._free(pVertices_mesh_in);
    Module._free(pEdgesSize_mesh_in);
    Module._free(pEdges_mesh_in);
    Module._free(pTrianglesSize_mesh_in);
    Module._free(pTriangles_mesh_in);
    Module._free(pCrackedEdgesSize_mesh_in);
    Module._free(pCrackedEdges_mesh_in);
    Module._free(pVerticesOnGeomEdgeSize_mesh_in);
    Module._free(pVerticesOnGeomEdge_mesh_in);
    Module._free(pVerticesOnGeomVertexSize_mesh_in);
    Module._free(pVerticesOnGeomVertex_mesh_in);
    Module._free(pEdgesOnGeomEdgeSize_mesh_in);
    Module._free(pEdgesOnGeomEdge_mesh_in);
    Module._free(pIssmSegmentsSize_mesh_in);
    Module._free(pIssmSegments_mesh_in);

    Module._free(pVerticesSize_geom_in);
    Module._free(pVertices_geom_in);
    Module._free(pEdgesSize_geom_in);
    Module._free(pEdges_geom_in);
    Module._free(pCornersSize_geom_in);
    Module._free(pCorners_geom_in);
    Module._free(pRequiredVerticesSize_geom_in);
    Module._free(pRequiredVertices_geom_in);
    Module._free(pRequiredEdgesSize_geom_in);
    Module._free(pRequiredEdges_geom_in);
    Module._free(pCrackedEdgesSize_geom_in);
    Module._free(pCrackedEdges_geom_in);
    Module._free(pSubDomainsSize_geom_in);
    Module._free(pSubDomains_geom_in);

    Module._free(phminVerticesSize);
    Module._free(phminVertices);
    Module._free(phmaxVerticesSize);
    Module._free(phmaxVertices);
    Module._free(phVertices);
    Module._free(pmetricSize);
    Module._free(pmetric);
    Module._free(pfieldSize);
    Module._free(pfield);
    Module._free(perrSize);
    Module._free(perr);

    return return_array;
}

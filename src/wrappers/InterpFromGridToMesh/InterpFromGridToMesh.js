/**
 * INTERPFROMGRIDTOMESH - interpolation from a grid onto a list of points
 *
 * This function interpolates a field defined on a grid to a list of points based on a bilinear interpolation.
 *
 * Usage:
 *	var data_mesh=InterpFromGridToMesh(xIn,yIn,dataIn,xMeshIn,yMeshIn,defaultValue,interpType);\
 *
 *	xIn,yIn						: coordinates of matrix data. (x and y must be in increasing order)
 *	dataIn						: matrix holding the data to be interpolated onto the mesh
 *	xMeshIn,yMeshIn				: coordinates of the points onto which we interpolate
 *	defaultValue 				: default value if no data is found (holes)
 *	interpType (optional) 		: interpolation type
 * 	dataMesh					: array of mesh interpolated data
 */
function InterpFromGridToMesh(xIn,yIn,dataIn,dataNumColsIn,dataNumRowsIn,xMeshIn,yMeshIn,defaultValue) {
	/*
		Variables
	*/
	//{{{
	var data 			= {};
	var dataMesh 		= {};
	var dataMeshPtr 	= {};
	var ddata			= {};
	var ddataHeap 		= {};
	var ddataPtr 		= {};
	var dx				= {};
	var dxHeap 			= {};
	var dxMesh 			= {};
	var dxMeshHeap 		= {};
	var dxMeshPtr 		= {};
	var dxPtr 			= {};
	var dy				= {};
	var dyHeap 			= {};
	var dyMesh			= {};
	var dyMeshHeap 		= {};
	var dyMeshPtr 		= {};
	var dyPtr 			= {};
	var interpType 		= '';
	var meshNumRows		= 0;
	var ndata			= {};
	var nods 			= 0;
	var nx 				= {};
	var nxMesh 			= {};
	var ny 				= {};
	var nyMesh 			= {};
	var pdataMesh 		= {};
	var x 				= {};
	var xMesh 			= {};
	var y				= {};
	var yMesh 			= {};
	//}}}


	/*
		Dynamic allocations
	*/
	//{{{

	/*
		Input
	*/
	//{{{
	dx 			= new Float64Array(xIn);
	nx 			= dx.length * dx.BYTES_PER_ELEMENT;
	dxPtr 		= Module._malloc(nx);
	dxHeap 		= new Uint8Array(Module.HEAPU8.buffer, dxPtr, nx);
	dxHeap.set(new Uint8Array(dx.buffer));
	x 			= dxHeap.byteOffset;

	dy 			= new Float64Array(yIn);
	ny 			= dy.length * dy.BYTES_PER_ELEMENT;
	dyPtr 		= Module._malloc(ny);
	dyHeap 		= new Uint8Array(Module.HEAPU8.buffer, dyPtr, ny);
	dyHeap.set(new Uint8Array(dy.buffer));
	y 			= dyHeap.byteOffset;

	ddata 		= new Float64Array(dataIn);
	ndata 		= ddata.length * ddata.BYTES_PER_ELEMENT;
	ddataPtr 	= Module._malloc(ndata);
	ddataHeap 	= new Uint8Array(Module.HEAPU8.buffer, ddataPtr, ndata);
	ddataHeap.set(new Uint8Array(ddata.buffer));
	data 		= ddataHeap.byteOffset;

	dxMesh 		= new Float64Array(xMeshIn);
	nxMesh 		= dxMesh.length * dxMesh.BYTES_PER_ELEMENT;
	dxMeshPtr 	= Module._malloc(nxMesh);
	dxMeshHeap 	= new Uint8Array(Module.HEAPU8.buffer, dxMeshPtr, nxMesh);
	dxMeshHeap.set(new Uint8Array(dxMesh.buffer));
	xMesh 		= dxMeshHeap.byteOffset;

	dyMesh 		= new Float64Array(yMeshIn);
	nyMesh 		= dyMesh.length * dyMesh.BYTES_PER_ELEMENT;
	dyMeshPtr 	= Module._malloc(nyMesh);
	dyMeshHeap 	= new Uint8Array(Module.HEAPU8.buffer, dyMeshPtr, nyMesh);
	dyMeshHeap.set(new Uint8Array(dyMesh.buffer));
	yMesh 		= dyMeshHeap.byteOffset;

	nods 		= xMeshIn.length;
	meshNumRows	= xMeshIn.length;


	/*
		Retrieve interpolation type
	*/
	//{{{
	if (arguments.length === 7) {
		interpType = arguments[6];
	} else {
		interpType = 'bilinear';
	}
	//}}}

	/*
		Output
	*/
	pdataMesh = Module._malloc(4);
	//}}}

	//}}}


	/*
		Declare InterpFromGridToMesh module
	*/
	//{{{
	InterpFromGridToMeshModule = Module.cwrap(
		'InterpFromGridToMeshModule',
		'number',
		[
			'number', // output : pdataMesh
			'number', // input	: x
			'number', // input	: y
			'number', // input 	: data
			'number', // input 	: xMesh
			'number', // input	: yMesh
			'number', // input 	: defaultValue
			'number', // input	: nods
			'number', // input	: dataNumRowsIn
			'number', // input	: dataNumColsIn
			'number', // input	: meshNumRows
			'string'  // input	: interpType
		]
	);
	//}}}


	/*
		Call InterpFromGridToMesh module
	*/
	//{{{
	InterpFromGridToMeshModule(
		pdataMesh,
		x,
		y,
		data,
		xMesh,
		yMesh,
		defaultValue,
		nods,
		dataNumRowsIn,
		dataNumColsIn,
		meshNumRows,
		interpType
	);
	//}}}


	/*
		Dynamic copying from heap
	*/
	//{{{
	dataMeshPtr	= Module.getValue(pdataMesh, 'i32');
	dataMesh	= Module.HEAPF64.slice(dataMeshPtr / 8, dataMeshPtr / 8 + nods);
	//}}}


	/*
		Free resources
	*/
	//{{{
	Module._free(pdataMesh);
	//}}}


	return dataMesh;
}

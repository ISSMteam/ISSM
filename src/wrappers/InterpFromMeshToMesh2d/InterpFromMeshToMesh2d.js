function InterpFromMeshToMesh2d(indexin,xin,yin,datain,x_interpin,y_interpin){

/* INTERFROMMESHTOMESH2D - interpolation from a 2d triangular mesh onto a list of point
  
  This function interpolates a field defined on a Delaunay triangulation onto a list of points.

  Usage:
	  var data_interp=InterpFromMeshToMesh2d(index,x,y,data,x_interp,y_interp); or
	  var data_interp=InterpFromMeshToMesh2d(index,x,y,data,x_interp,y_interp,default_value);

	  index             : index of the mesh where data is defined
	  x,y               : coordinates of the nodes where data is defined
	  data              : matrix holding the data to be interpolated onto the mesh. (one column per field)
	  x_interp,y_interp : coordinates of the points onto which we interpolate.
	  default_value     : default value if point is outsite of triangulation (instead of linear interpolation)
	  data_interp       : vector of mesh interpolated data.

*/

	/*Figure out default_value: */
	if (arguments.length==7)default_value=arguments[6];
	else default_value=0;
	
	//Dynamic allocations: {{{
	//Retrieve elements and allocate on Module heap: 
	
	//input
	
	var dindex=new Int32Array(MatrixToList(indexin)); var nindex=dindex.length * dindex.BYTES_PER_ELEMENT;
	var dindexPtr= Module._malloc(nindex); var indexHeap = new Uint8Array(Module.HEAPU8.buffer,dindexPtr,nindex);
	indexHeap.set(new Uint8Array(dindex.buffer)); var index=indexHeap.byteOffset;

	var dx=new Float64Array(xin); var nx=dx.length * dx.BYTES_PER_ELEMENT;
	var dxPtr= Module._malloc(nx); var dxHeap = new Uint8Array(Module.HEAPU8.buffer,dxPtr,nx);
	dxHeap.set(new Uint8Array(dx.buffer)); var x=dxHeap.byteOffset;
	
	var dy=new Float64Array(yin); var ny=dy.length * dy.BYTES_PER_ELEMENT;
	var dyPtr= Module._malloc(nx); var dyHeap = new Uint8Array(Module.HEAPU8.buffer,dyPtr,ny);
	dyHeap.set(new Uint8Array(dy.buffer)); var y=dyHeap.byteOffset;
	
	var ddata=new Float64Array(datain); var ndata=ddata.length * ddata.BYTES_PER_ELEMENT;
	var ddataPtr= Module._malloc(ndata); var ddataHeap = new Uint8Array(Module.HEAPU8.buffer,ddataPtr,ndata);
	ddataHeap.set(new Uint8Array(ddata.buffer)); var data=ddataHeap.byteOffset;
	
	var dx_interp=new Float64Array(x_interpin); var nx_interp=dx_interp.length * dx_interp.BYTES_PER_ELEMENT;
	var dx_interpPtr= Module._malloc(nx_interp); var dx_interpHeap = new Uint8Array(Module.HEAPU8.buffer,dx_interpPtr,nx_interp);
	dx_interpHeap.set(new Uint8Array(dx_interp.buffer)); var x_interp=dx_interpHeap.byteOffset;
	
	var dy_interp=new Float64Array(y_interpin); var ny_interp=dy_interp.length * dy_interp.BYTES_PER_ELEMENT;
	var dy_interpPtr= Module._malloc(ny_interp); var dy_interpHeap = new Uint8Array(Module.HEAPU8.buffer,dy_interpPtr,ny_interp);
	dy_interpHeap.set(new Uint8Array(dy_interp.buffer)); var y_interp=dy_interpHeap.byteOffset;
	
	nel=indexin.length;
	nods=xin.length;
	nods_interp=x_interpin.length;

	//output
	var data_interp;
	var pdata_interp= Module._malloc(4); 
	//}}}

	//Declare InterpFromMeshToMesh2d module: 
	InterpFromMeshToMesh2dModule = Module.cwrap('InterpFromMeshToMesh2dModule','number',['number','number','number','number','number','number','number','number','number','number','number']);
	
	//Call InterpFromMeshToMesh2d module: 
	InterpFromMeshToMesh2dModule(pdata_interp,index,x,y,data,x_interp,y_interp,nel,nods,nods_interp,default_value);
	
	/*Dynamic copying from heap: {{{*/
	//recover mesh: 
	var data_interpptr = Module.getValue(pdata_interp,'i32');
	data_interp = Module.HEAPF64.slice(data_interpptr /8, data_interpptr/8 + nods_interp);
	/*}}}*/

	/*Free resources: */
	Module._free(pdata_interp); 

	return data_interp;
}

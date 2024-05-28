function ContourToMesh(indexin,xin,yin,contour,interptype,edgevalue){
/* CONTOURTOMESH - Flag the elements or nodes inside a contour;
	
	      Usage: ;
	         [in_nod,in_elem]=ContourToMesh(index,x,y,contourname,interptype,edgevalue);
	
	         index,x,y: mesh triangulation
	         contourname: name of .exp file containing the contours
	         interptype: string definining type of interpolation ('element', or 'node')
	         edgevalue: integer (0, 1 or 2) defining the value associated to the nodes on the edges of the polygons.
	         in_nod: vector of flags (0 or 1), of size nods if interptype is set to 'node' or 'element and node',
	            or of size 0 otherwise.
	         in_elem: vector of flags (0 or 1), of size nel if interptype is set to 'element' or 'element and node', 
	            or of size 0 otherwise.
	
	      Example: 
	         in_nod=ContourToMesh(md.elements,md.x,md.y,'Contour.exp','node',1)
	         in_elements=ContourToMesh(md.elements,md.x,md.y,'Contour.exp','element',0)
	         return_values=ContourToMesh(md.elements,md.x,md.y,'Contour.exp','element and node',0); in_nodes=return_values[0]; in_elements=return_values[1];
*/

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
	
	var dcontourx=new Float64Array(contour['x']); var nx=dcontourx.length * dcontourx.BYTES_PER_ELEMENT;
	var dcontourxPtr= Module._malloc(nx); var contourxHeap = new Uint8Array(Module.HEAPU8.buffer,dcontourxPtr,nx);
	contourxHeap.set(new Uint8Array(dcontourx.buffer)); var contourx=contourxHeap.byteOffset;

	var dcontoury=new Float64Array(contour['y']); var ny=dcontoury.length * dcontoury.BYTES_PER_ELEMENT;
	var dcontouryPtr = Module._malloc(ny); var contouryHeap = new Uint8Array(Module.HEAPU8.buffer,dcontouryPtr,ny);
	contouryHeap.set(new Uint8Array(dcontoury.buffer)); var contoury=contouryHeap.byteOffset;
	
	nel=indexin.length;
	nods=xin.length;
	contour_nods=dcontourx.length;

	//output
	var in_nod;
	var pin_nod= Module._malloc(4); 
	var in_nel;
	var pin_nel= Module._malloc(4); 
	//}}}

	//Declare ContourToMesh module: 
	ContourToMeshModule = Module.cwrap('ContourToMeshModule','number',['number','number','number','number','number','number','number','string','number','number','number']);
	
	//Call ContourToMesh module: 
	ContourToMeshModule(pin_nod,pin_nel,index,x,y,contourx,contoury,interptype,nel, nods, contour_nods, edgevalue);

	/*Dynamic copying from heap: {{{*/
	if(interptype == 'node'){
		var in_nodptr = Module.getValue(pin_nod,'i32');
		in_nod = Module.HEAPF64.slice(in_nodptr /8, in_nodptr/8 + nods);
	}
	else if (interptype == 'element'){
		var in_nelptr = Module.getValue(pin_nel,'i32');
		in_nel = Module.HEAPF64.slice(in_nelptr /8, in_nelptr/8 + nel);
	}
	else if (interptype == 'element and node'){
		var in_nodptr = Module.getValue(pin_nod,'i32');
		in_nod = Module.HEAPF64.slice(in_nodptr /8, in_nodptr/8 + nods);
		var in_nelptr = Module.getValue(pin_nel,'i32');
		in_nel = Module.HEAPF64.slice(in_nelptr /8, in_nelptr/8 + nel);
	}
	else throw Error('ContourToMeshModule error message: wrong interpolation type!');
	/*}}}*/

	/*Free resources: */
	Module._free(pin_nod); 
	Module._free(pin_nel); 
	
	if(interptype == 'node'){
		return in_nod;
	}
	else if (interptype == 'element'){
		return in_nel;
	}
	else if (interptype == 'element and node'){
		return [in_nod,in_nel];
	}
}

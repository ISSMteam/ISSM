function NodeConnectivity(elementsin,nods){
/*NodeConnectivity 
	   usage: var md.mesh.vertexconnectivity = NodeConnectivity(md.mesh.elements,md.mesh.numberofvertices);
*/

	//Dynamic allocations: {{{
	//Retrieve elements and allocate on Module heap: 
	
	//input
	var delements=new Int32Array(MatrixToList(elementsin)); var nelements=delements.length * delements.BYTES_PER_ELEMENT;
	var delementsPtr= Module._malloc(nelements); var elementsHeap = new Uint8Array(Module.HEAPU8.buffer,delementsPtr,nelements);
	elementsHeap.set(new Uint8Array(delements.buffer)); var elements=elementsHeap.byteOffset;

	//output
	var width,connectivitylinear,connectivity;
	var pwidth= Module._malloc(4); 
	var pnods= Module._malloc(4); 
	var pconnectivity= Module._malloc(4); 
	var nels=elementsin.length;
	//}}}

	//Declare NodeConnectivity module: 
	NodeConnectivityModule = Module.cwrap('NodeConnectivityModule','number',['number','number','number','number']);
	
	//Call NodeConnectivity module: 
	NodeConnectivityModule(pconnectivity,pnods,pwidth,elements,nels,nods);
	
	/*Dynamic copying from heap: {{{*/
	//recover mesh: 
	width = Module.getValue(pwidth, 'i32');
	var connectivityptr = Module.getValue(pconnectivity,'i32');
	connectivitylinear = Module.HEAPF64.slice(connectivityptr /8, connectivityptr/8 + nods*width);
	connectivity = ListToMatrix(connectivitylinear,width);
	/*}}}*/

	/*Free resources: */
	Module._free(pconnectivity); 
	Module._free(connectivitylinear); 
	Module._free(pwidth); 
	Module._free(pnods); 

	return connectivity;
}

function ElementConnectivity(elementsin,nodeconnectivityin){
/*ElementConnectivity 
	   usage: var md.mesh.elementconnectivity= ElementConnectivity(md.mesh.elements,md.mesh.vertexconnectivity);
*/

	//Dynamic allocations: {{{
	//Retrieve elements and allocate on Module heap: 
	
	//input
	nel=elementsin.length;
	nods=nodeconnectivityin.length;
	width=nodeconnectivityin[0].length;
	
	var delements=new Int32Array(MatrixToList(elementsin)); var nelements=delements.length * delements.BYTES_PER_ELEMENT;
	var delementsPtr= Module._malloc(nelements); var elementsHeap = new Uint8Array(Module.HEAPU8.buffer,delementsPtr,nelements);
	elementsHeap.set(new Uint8Array(delements.buffer)); var elements=elementsHeap.byteOffset;
	
	var dnodeconnectivity=new Int32Array(MatrixToList(nodeconnectivityin)); var nnodeconnectivity=dnodeconnectivity.length * dnodeconnectivity.BYTES_PER_ELEMENT;
	var dnodeconnectivityPtr= Module._malloc(nnodeconnectivity); var nodeconnectivityHeap = new Uint8Array(Module.HEAPU8.buffer,dnodeconnectivityPtr,nnodeconnectivity);
	nodeconnectivityHeap.set(new Uint8Array(dnodeconnectivity.buffer)); var nodeconnectivity=nodeconnectivityHeap.byteOffset;

	//output
	var connectivitylinear,connectivity;
	var pconnectivity= Module._malloc(4); 
	//}}}

	//Declare ElementConnectivity module: 
	ElementConnectivityModule = Module.cwrap('ElementConnectivityModule','number',['number','number','number','number','number','number']);
	
	//Call ElementConnectivity module: 
	ElementConnectivityModule(pconnectivity,elements, nodeconnectivity, nel, nods, width);
	
	/*Dynamic copying from heap: {{{*/
	//recover mesh: 
	var connectivityptr = Module.getValue(pconnectivity,'i32');
	connectivitylinear = Module.HEAPF64.slice(connectivityptr /8, connectivityptr/8 + nel*3);
	connectivity = ListToMatrix(connectivitylinear,3);
	/*}}}*/

	/*Free resources: */
	Module._free(pconnectivity); 
	Module._free(connectivitylinear); 

	return connectivity;
}

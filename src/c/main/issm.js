function issm(binaryin){
/*issm 
	   usage: var output = issm(input);
	      where: input is a typed array buffer created by marshall and output 
		  is a binary buffer to be read by loadresultsfromcluster.
*/

	//input
	var dbinary=new Float64Array(binaryin); var nb=dbinary.length * dbinary.BYTES_PER_ELEMENT;
	var dbinaryPtr= Module._malloc(nb); var binHeap = new Uint8Array(Module.HEAPU8.buffer,dbinaryPtr,nb);
	binHeap.set(new Uint8Array(dbinary.buffer)); var binary=binHeap.byteOffset;

	//Declare module: 
	issm= Module.cwrap('main','number',['number','number']);
	
	//Call issm:
	var output = issm(binary, 'null');
	
	return output;
}

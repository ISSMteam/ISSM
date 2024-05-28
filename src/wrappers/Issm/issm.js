function issm(fid,toolkitstring,solutionstring,modelname){
/*issm 
	   usage: var output = issm(fid,toolkitstring);
	      where: fid is a pointer to a memory buffer created by marshall, toolkitstring is a toolkits 
		  string created by ToolkitsToFile and 
		  output is a binary buffer to be read by loadresultsfromcluster.
*/
	
	/*variables: */
	var poutput,output,poutputsize,outputsize;
	var dbinaryPtr,binHeap,binary,binarybuffer,nb

	/*recover input buffer: */
	binarybuffer=fid.rawbuffer(); //binarybuffer is now an arraybuffer
	nb = fid.ptr; //size of array buffer in bytes.

	/*dyanmically allocate the raw buffer onto the Module heap: */
	dbinaryPtr= Module._malloc(nb); binHeap = new Uint8Array(Module.HEAPU8.buffer,dbinaryPtr,nb);
	binHeap.set(new Uint8Array(binarybuffer)); binary=binHeap.byteOffset;

	/*allocate output pointers: */
	poutputsize = Module._malloc(4); 
	poutput = Module._malloc(4); 

	//Declare TriMesh module: 
	issmmodule= Module.cwrap('IssmModule','number',['number','number','number','number','string','string','string']);
	
	//Call issm:
	issmmodule(poutput, poutputsize,binary, nb, toolkitstring,solutionstring,modelname);

	//recover outputs from pointers: 
	var outputsize = Module.getValue(poutputsize,'i32');
	
	var outputptr = Module.getValue(poutput,'i32');
	output = Module.HEAP8.slice(outputptr, outputptr + outputsize);
	
	return [output,outputsize];
}

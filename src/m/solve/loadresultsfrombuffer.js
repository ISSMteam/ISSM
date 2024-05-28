function loadresultsfrombuffer(md,buffer,buffersize){
//LOADRESULTSFROMBUFFER - load results of solution sequence from memory buffer 
//
//   Usage:
//      loadresultsfrombuffer(md,buffer,buffersize);

	//check number of inputs/outputs
	if (arguments.length!=3) throw Error('loadresultsfrombuffer error message: wrong usage!');

	if (!md.qmu.isdakota){

		//initialize md.results if not a structure yet
		if (MapIsEmpty(md.results)) md.results={};

		//load results: 
		structure=parseresultsfrombuffer(md,buffer,buffersize);

		//load structure onto results: 
		md.results=structure;

		return md;

	}
	else throw Error('loadresultsfrombuffer error message: qmu results not supported yet!');
}

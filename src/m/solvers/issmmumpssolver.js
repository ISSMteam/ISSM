function issmmumpssolver(){
//ISSMMUMPSSOLVER - 
//
//   Usage:
//      options=issmmumpssolver(varargin);

	//Retrieve options provided in varargin. First convert arguments to array:
	var args = Array.prototype.slice.call(arguments);

	//Then process options
	var  options = new pairoptions(args.slice(1,args.length));

	//default issmoptions options
	var issmoptions={};
	issmoptions['toolkit']='issm';
	issmoptions['mat_type']=options.getfieldvalue('mat_type','mpisparse');
	issmoptions['vec_type']=options.getfieldvalue('vec_type','mpi');
	issmoptions['solver_type']=options.getfieldvalue('solver_type','mumps');

	return issmoptions;
}

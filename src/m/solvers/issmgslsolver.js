function issmgslsolver(){
//ISSMSOLVER - 
//
//   Usage:
//      options=issmsolver(varargin);

	//Retrieve options provided in varargin. First convert arguments to array:
	var args = Array.prototype.slice.call(arguments);

	//Then process options
	var  options = new pairoptions(args.slice(1,args.length));

	//default issmoptions options
	var issmoptions={};
	issmoptions['toolkit']='issm';
	issmoptions['mat_type']=options.getfieldvalue('mat_type','dense');
	issmoptions['vec_type']=options.getfieldvalue('vec_type','seq');
	issmoptions['solver_type']=options.getfieldvalue('solver_type','gsl');

	return issmoptions;
}

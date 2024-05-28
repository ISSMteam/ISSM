function mumpsoptions(){
//MUMPSOPTIONS - 
//
//   Usage:
//      options=mumpsoptions(varargin);

	//Retrieve options provided in varargin. First convert arguments to array:
	var args = Array.prototype.slice.call(arguments);

	//Then process options
	var  options = new pairoptions(args.slice(1,args.length));

	//default issmoptions options
	var mumpsoptions={};
	mumpsoptions['toolkit']='petsc';
	mumpsoptions['mat_type']=options.getfieldvalue('mat_type','mpiaij');
	mumpsoptions['ksp_type']=options.getfieldvalue('ksp_type','preonly');
	mumpsoptions['pc_type']=options.getfieldvalue('pc_type','lu');
	mumpsoptions['pc_factor_mat_solver_type']=options.getfieldvalue('pc_factor_mat_solver_type','mumps');
	mumpsoptions['mat_mumps_icntl_14']=options.getfieldvalue('mat_mumps_icntl_14',120);
	mumpsoptions['mat_mumps_icntl_28']=options.getfieldvalue('mat_mumps_icntl_28',2);
	mumpsoptions['mat_mumps_icntl_29']=options.getfieldvalue('mat_mumps_icntl_29',2);

	return mumpsoptions;
}

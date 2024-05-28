function solverOptions=bcgslbjacobioptions(varargin)

%retrieve options provided in varargin
options=pairoptions(varargin{:});
solverOptions=struct();
solverOptions.toolkit='petsc';
solverOptions.mat_type=getfieldvalue(options, 'mat_type','mpiaij');
solverOptions.ksp_type=getfieldvalue(options, 'ksp_type','bcgsl');
solverOptions.pc_type=getfieldvalue(options, 'pc_type',  'bjacobi');
solverOptions.ksp_max_it=getfieldvalue(options,'ksp_max_it',300);
solverOptions.ksp_rtol=getfieldvalue(options,'ksp_rtol',1e-13);

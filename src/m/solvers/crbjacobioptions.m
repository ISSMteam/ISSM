function solverOptions=crbjacobioptions(varargin)

%retrieve options provided in varargin
options=pairoptions(varargin{:});
solverOptions=struct();
solverOptions.toolkit='petsc';
solverOptions.mat_type=getfieldvalue(options, 'mat_type','mpiaij');
solverOptions.ksp_type=getfieldvalue(options, 'ksp_type','cr');
solverOptions.pc_type=getfieldvalue(options, 'pc_type',  'bjacobi');


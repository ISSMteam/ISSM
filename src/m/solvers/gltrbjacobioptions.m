function solverOptions=gltrbjacobioptions(varargin)
%GLTRBJACOBIOPTIONS - PETSc solver options using the generalized Lanczos trust-region (gltr) Krylov method with block Jacobi preconditioning
%
%   Usage:
%      solverOptions=gltrbjacobioptions(varargin);

%retrieve options provided in varargin
options=pairoptions(varargin{:});
solverOptions=struct();
solverOptions.toolkit='petsc';
solverOptions.mat_type=getfieldvalue(options, 'mat_type','mpiaij');
solverOptions.ksp_type=getfieldvalue(options, 'ksp_type','gltr');
solverOptions.pc_type=getfieldvalue(options, 'pc_type',  'bjacobi');


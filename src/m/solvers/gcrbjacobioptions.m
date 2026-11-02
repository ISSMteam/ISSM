function solverOptions=gcrbjacobioptions(varargin)
%GCRBJACOBIOPTIONS - PETSc solver options using the generalized conjugate residual (gcr) Krylov method with block Jacobi preconditioning
%
%   Usage:
%      solverOptions=gcrbjacobioptions(varargin);

%retrieve options provided in varargin
options=pairoptions(varargin{:});
solverOptions=struct();
solverOptions.toolkit='petsc';
solverOptions.mat_type=getfieldvalue(options, 'mat_type','mpiaij');
solverOptions.ksp_type=getfieldvalue(options, 'ksp_type','gcr');
solverOptions.pc_type=getfieldvalue(options, 'pc_type',  'bjacobi');

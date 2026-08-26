function solverOptions=cgjacobioptions(varargin)
%CGJACOBIOPTIONS - define PETSc solver options for the Conjugate Gradient (CG) Krylov method with Jacobi preconditioning
%
%   Usage:
%      solverOptions=cgjacobioptions(varargin);

%retrieve options provided in varargin
options=pairoptions(varargin{:});
solverOptions=struct();
solverOptions.toolkit='petsc';
solverOptions.mat_type=getfieldvalue(options, 'mat_type','mpiaij');
solverOptions.ksp_type=getfieldvalue(options, 'ksp_type','cg');
solverOptions.pc_type=getfieldvalue(options, 'pc_type',  'jacobi');


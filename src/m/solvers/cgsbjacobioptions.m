function solverOptions=cgsbjacobioptions(varargin)
%CGSBJACOBIOPTIONS - define PETSc solver options for the Conjugate Gradient Squared (CGS) Krylov method with Block Jacobi preconditioning
%
%   Usage:
%      solverOptions=cgsbjacobioptions(varargin);

%retrieve options provided in varargin
options=pairoptions(varargin{:});
solverOptions=struct();
solverOptions.toolkit='petsc';
solverOptions.mat_type=getfieldvalue(options, 'mat_type','mpiaij');
solverOptions.ksp_type=getfieldvalue(options, 'ksp_type','cgs');
solverOptions.pc_type=getfieldvalue(options, 'pc_type',  'bjacobi');


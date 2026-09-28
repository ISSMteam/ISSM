function solverOptions=cgnebjacobioptions(varargin)
%CGNEBJACOBIOPTIONS - define PETSc solver options for the Conjugate Gradient on the Normal Equations (CGNE) Krylov method with Block Jacobi preconditioning
%
%   Usage:
%      solverOptions=cgnebjacobioptions(varargin);

%retrieve options provided in varargin
options=pairoptions(varargin{:});
solverOptions=struct();
solverOptions.toolkit='petsc';
solverOptions.mat_type=getfieldvalue(options, 'mat_type','mpiaij');
solverOptions.ksp_type=getfieldvalue(options, 'ksp_type','cgne');
solverOptions.pc_type=getfieldvalue(options, 'pc_type',  'bjacobi');


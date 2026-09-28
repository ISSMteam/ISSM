function solverOptions=minresbjacobioptions(varargin)
%MINRESBJACOBIOPTIONS - return PETSc options for the MINRES (minimum residual) Krylov solver with a block Jacobi preconditioner
%
%   Usage:
%      solverOptions=minresbjacobioptions();

%retrieve options provided in varargin
options=pairoptions(varargin{:});
solverOptions=struct();
solverOptions.toolkit='petsc';
solverOptions.mat_type=getfieldvalue(options, 'mat_type','mpiaij');
solverOptions.ksp_type=getfieldvalue(options, 'ksp_type','minres');
solverOptions.pc_type=getfieldvalue(options, 'pc_type',  'bjacobi');


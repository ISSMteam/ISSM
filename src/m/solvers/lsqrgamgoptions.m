function solverOptions=lsqrgamgoptions(varargin)
%LSQRGAMGOPTIONS - return PETSc options for the LSQR Krylov solver with a GAMG (geometric-algebraic multigrid) preconditioner
%
%   Usage:
%      solverOptions=lsqrgamgoptions();

%retrieve options provided in varargin
options=pairoptions(varargin{:});
solverOptions=struct();
solverOptions.toolkit='petsc';
solverOptions.mat_type=getfieldvalue(options, 'mat_type','mpiaij');
solverOptions.ksp_type=getfieldvalue(options, 'ksp_type','lsqr');
solverOptions.pc_type=getfieldvalue(options, 'pc_type',  'gamg');


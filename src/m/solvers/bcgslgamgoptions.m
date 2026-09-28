function solverOptions=bcgslgamgoptions(varargin)
%BCGSLGAMGOPTIONS - return BiCGStab(L) PETSc options with a geometric algebraic multigrid (GAMG) preconditioner
%
%   Usage:
%      options=bcgslgamgoptions;

%retrieve options provided in varargin
options=pairoptions(varargin{:});
solverOptions=struct();
solverOptions.toolkit='petsc';
solverOptions.mat_type=getfieldvalue(options, 'mat_type','mpiaij');
solverOptions.ksp_type=getfieldvalue(options, 'ksp_type','bcgsl');
solverOptions.pc_type=getfieldvalue(options, 'pc_type',  'gamg');

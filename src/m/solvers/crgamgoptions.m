function solverOptions=crgamgoptions(varargin)

%CRGAMGOPTIONS - define PETSc solver options for the Conjugate Residual (CR) Krylov method with Geometric-Algebraic Multigrid (GAMG) preconditioning
%
%   Usage:
%      solverOptions=crgamgoptions(varargin);

%retrieve options provided in varargin
options=pairoptions(varargin{:});
solverOptions=struct();
solverOptions.toolkit='petsc';
solverOptions.mat_type=getfieldvalue(options, 'mat_type','mpiaij');
solverOptions.ksp_type=getfieldvalue(options, 'ksp_type','cr');
solverOptions.pc_type=getfieldvalue(options, 'pc_type',  'gamg');


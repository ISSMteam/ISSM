function solverOptions=lcdgamgoptions(varargin)
%LCDGAMGOPTIONS - return PETSc options for the LCD (left conjugate direction) Krylov solver with a GAMG (geometric-algebraic multigrid) preconditioner
%
%   Usage:
%      solverOptions=lcdgamgoptions();

%retrieve options provided in varargin
options=pairoptions(varargin{:});
solverOptions=struct();
solverOptions.toolkit='petsc';
solverOptions.mat_type=getfieldvalue(options, 'mat_type','mpiaij');
solverOptions.ksp_type=getfieldvalue(options, 'ksp_type','lcd');
solverOptions.pc_type=getfieldvalue(options, 'pc_type',  'gamg');


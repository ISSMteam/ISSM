function solverOptions=bcgslasmoptions(varargin)
%BCGSLASMOPTIONS - return BiCGStab(L) PETSc options with an additive Schwartz method (ASM) preconditioner
%
%   Usage:
%      options=bcgslasmoptions;

%retrieve options provided in varargin
options=pairoptions(varargin{:});
solverOptions=struct();
solverOptions.toolkit='petsc';
solverOptions.mat_type=getfieldvalue(options, 'mat_type','mpiaij');
solverOptions.ksp_type=getfieldvalue(options, 'ksp_type','bcgsl');
solverOptions.pc_type=getfieldvalue(options, 'pc_type',  'asm');

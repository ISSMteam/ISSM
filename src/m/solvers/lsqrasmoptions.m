function solverOptions=lsqrasmoptions(varargin)
%LSQRASMOPTIONS - return PETSc options for the LSQR Krylov solver with an ASM (Additive Schwarz Method) preconditioner
%
%   Usage:
%      solverOptions=lsqrasmoptions();

%retrieve options provided in varargin
options=pairoptions(varargin{:});
solverOptions=struct();
solverOptions.toolkit='petsc';
solverOptions.mat_type=getfieldvalue(options, 'mat_type','mpiaij');
solverOptions.ksp_type=getfieldvalue(options, 'ksp_type','lsqr');
solverOptions.pc_type=getfieldvalue(options, 'pc_type',  'asm');


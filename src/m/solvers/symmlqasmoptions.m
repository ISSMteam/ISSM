function solverOptions=symmlqasmoptions(varargin)
%SYMMLQASMOPTIONS - return PETSc options for the SYMMLQ solver with an ASM (Additive Schwarz Method) preconditioner
%
%   Usage:
%      solverOptions=symmlqasmoptions(varargin);

%retrieve options provided in varargin
options=pairoptions(varargin{:});
solverOptions=struct();
solverOptions.toolkit='petsc';
solverOptions.mat_type=getfieldvalue(options, 'mat_type','mpiaij');
solverOptions.ksp_type=getfieldvalue(options, 'ksp_type','symmlq');
solverOptions.pc_type=getfieldvalue(options, 'pc_type',  'asm');


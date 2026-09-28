function solverOptions=tfqmrasmoptions(varargin)
%TFQMRASMOPTIONS - return PETSc options for the TFQMR solver with an ASM (Additive Schwarz Method) preconditioner
%
%   Usage:
%      solverOptions=tfqmrasmoptions(varargin);

%retrieve options provided in varargin
options=pairoptions(varargin{:});
solverOptions=struct();
solverOptions.toolkit='petsc';
solverOptions.mat_type=getfieldvalue(options, 'mat_type','mpiaij');
solverOptions.ksp_type=getfieldvalue(options, 'ksp_type','tfqmr');
solverOptions.pc_type=getfieldvalue(options, 'pc_type',  'asm');

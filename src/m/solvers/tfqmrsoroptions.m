function solverOptions=tfqmrsoroptions(varargin)
%TFQMRSOROPTIONS - return PETSc options for the TFQMR solver with a SOR (Successive Over-Relaxation) preconditioner
%
%   Usage:
%      solverOptions=tfqmrsoroptions(varargin);

%retrieve options provided in varargin
options=pairoptions(varargin{:});
solverOptions=struct();
solverOptions.toolkit='petsc';
solverOptions.mat_type=getfieldvalue(options, 'mat_type','mpiaij');
solverOptions.ksp_type=getfieldvalue(options, 'ksp_type','tfqmr');
solverOptions.pc_type=getfieldvalue(options, 'pc_type',  'sor');


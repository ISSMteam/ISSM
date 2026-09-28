function solverOptions=richardsonsoroptions(varargin)
%RICHARDSONSOROPTIONS - return PETSc options for the Richardson iteration Krylov solver with an SOR (Successive Over-Relaxation) preconditioner
%
%   Usage:
%      solverOptions=richardsonsoroptions();

%retrieve options provided in varargin
options=pairoptions(varargin{:});
solverOptions=struct();
solverOptions.toolkit='petsc';
solverOptions.mat_type=getfieldvalue(options, 'mat_type','mpiaij');
solverOptions.ksp_type=getfieldvalue(options, 'ksp_type','richardson');
solverOptions.pc_type=getfieldvalue(options, 'pc_type',  'sor');


function solverOptions=lsqrgasmoptions(varargin)
%LSQRGASMOPTIONS - return PETSc options for the LSQR Krylov solver with a GASM (Generalized Additive Schwarz Method) preconditioner
%
%   Usage:
%      solverOptions=lsqrgasmoptions();

%retrieve options provided in varargin
options=pairoptions(varargin{:});
solverOptions=struct();
solverOptions.toolkit='petsc';
solverOptions.mat_type=getfieldvalue(options, 'mat_type','mpiaij');
solverOptions.ksp_type=getfieldvalue(options, 'ksp_type','lsqr');
solverOptions.pc_type=getfieldvalue(options, 'pc_type',  'gasm');


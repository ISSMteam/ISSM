function solverOptions=richardsongasmoptions(varargin)
%RICHARDSONGASMOPTIONS - return PETSc options for the Richardson iteration Krylov solver with a GASM (Generalized Additive Schwarz Method) preconditioner
%
%   Usage:
%      solverOptions=richardsongasmoptions();

%retrieve options provided in varargin
options=pairoptions(varargin{:});
solverOptions=struct();
solverOptions.toolkit='petsc';
solverOptions.mat_type=getfieldvalue(options, 'mat_type','mpiaij');
solverOptions.ksp_type=getfieldvalue(options, 'ksp_type','richardson');
solverOptions.pc_type=getfieldvalue(options, 'pc_type',  'gasm');


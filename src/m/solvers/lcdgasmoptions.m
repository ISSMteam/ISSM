function solverOptions=lcdgasmoptions(varargin)
%LCDGASMOPTIONS - return PETSc options for the LCD (left conjugate direction) Krylov solver with a GASM (Generalized Additive Schwarz Method) preconditioner
%
%   Usage:
%      solverOptions=lcdgasmoptions();

%retrieve options provided in varargin
options=pairoptions(varargin{:});
solverOptions=struct();
solverOptions.toolkit='petsc';
solverOptions.mat_type=getfieldvalue(options, 'mat_type','mpiaij');
solverOptions.ksp_type=getfieldvalue(options, 'ksp_type','lcd');
solverOptions.pc_type=getfieldvalue(options, 'pc_type',  'gasm');


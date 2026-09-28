function solverOptions=stcggasmoptions(varargin)
%STCGGASMOPTIONS - return PETSc options for the STCG solver with a GASM (Generalized Additive Schwarz Method) preconditioner
%
%   Usage:
%      solverOptions=stcggasmoptions(varargin);

%retrieve options provided in varargin
options=pairoptions(varargin{:});
solverOptions=struct();
solverOptions.toolkit='petsc';
solverOptions.mat_type=getfieldvalue(options, 'mat_type','mpiaij');
solverOptions.ksp_type=getfieldvalue(options, 'ksp_type','stcg');
solverOptions.pc_type=getfieldvalue(options, 'pc_type',  'gasm');


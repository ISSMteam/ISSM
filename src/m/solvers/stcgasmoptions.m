function solverOptions=stcgasmoptions(varargin)
%STCGASMOPTIONS - return PETSc options for the STCG (Steihaug-Toint conjugate gradient) Krylov solver with an ASM (Additive Schwarz Method) preconditioner
%
%   Usage:
%      solverOptions=stcgasmoptions();

%retrieve options provided in varargin
options=pairoptions(varargin{:});
solverOptions=struct();
solverOptions.toolkit='petsc';
solverOptions.mat_type=getfieldvalue(options, 'mat_type','mpiaij');
solverOptions.ksp_type=getfieldvalue(options, 'ksp_type','stcg');
solverOptions.pc_type=getfieldvalue(options, 'pc_type',  'asm');


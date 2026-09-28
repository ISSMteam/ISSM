function solverOptions=stcgjacobioptions(varargin)
%STCGJACOBIOPTIONS - return PETSc options for the STCG solver with a Jacobi preconditioner
%
%   Usage:
%      solverOptions=stcgjacobioptions(varargin);

%retrieve options provided in varargin
options=pairoptions(varargin{:});
solverOptions=struct();
solverOptions.toolkit='petsc';
solverOptions.mat_type=getfieldvalue(options, 'mat_type','mpiaij');
solverOptions.ksp_type=getfieldvalue(options, 'ksp_type','stcg');
solverOptions.pc_type=getfieldvalue(options, 'pc_type',  'jacobi');


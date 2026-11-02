function solverOptions=symmlqpbjacobioptions(varargin)
%SYMMLQPBJACOBIOPTIONS - return PETSc options for the SYMMLQ solver with a point-block Jacobi preconditioner
%
%   Usage:
%      solverOptions=symmlqpbjacobioptions(varargin);

%retrieve options provided in varargin
options=pairoptions(varargin{:});
solverOptions=struct();
solverOptions.toolkit='petsc';
solverOptions.mat_type=getfieldvalue(options, 'mat_type','mpiaij');
solverOptions.ksp_type=getfieldvalue(options, 'ksp_type','symmlq');
solverOptions.pc_type=getfieldvalue(options, 'pc_type',  'pbjacobi');


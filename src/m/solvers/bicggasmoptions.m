function solverOptions=bicggasmoptions(varargin)
%BICGGASMOPTIONS - define PETSc solver options for the BiConjugate Gradient (BiCG) Krylov method with Generalized Additive Schwarz Method (GASM) preconditioning
%
%   Usage:
%      solverOptions=bicggasmoptions(varargin);

%retrieve options provided in varargin
options=pairoptions(varargin{:});
solverOptions=struct();
solverOptions.toolkit='petsc';
solverOptions.mat_type=getfieldvalue(options, 'mat_type','mpiaij');
solverOptions.ksp_type=getfieldvalue(options, 'ksp_type','bicg');
solverOptions.pc_type=getfieldvalue(options, 'pc_type',  'gasm');

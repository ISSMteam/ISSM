function solverOptions=cgsoroptions(varargin)

%CGSOROPTIONS - define PETSc solver options for the Conjugate Gradient (CG) Krylov method with Successive Over-Relaxation (SOR) preconditioning
%
%   Usage:
%      solverOptions=cgsoroptions(varargin);

%retrieve options provided in varargin
options=pairoptions(varargin{:});
solverOptions=struct();
solverOptions.toolkit='petsc';
solverOptions.mat_type=getfieldvalue(options, 'mat_type','mpiaij');
solverOptions.ksp_type=getfieldvalue(options, 'ksp_type','cg');
solverOptions.pc_type=getfieldvalue(options, 'pc_type',  'sor');


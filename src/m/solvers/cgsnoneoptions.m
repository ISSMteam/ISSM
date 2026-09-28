function solverOptions=cgsmgoptions(varargin)
%CGSNONEOPTIONS - define PETSc solver options for the Conjugate Gradient Squared (CGS) Krylov method with no preconditioning
%
%   Usage:
%      solverOptions=cgsnoneoptions(varargin);

%retrieve options provided in varargin
options=pairoptions(varargin{:});
solverOptions=struct();
solverOptions.toolkit='petsc';
solverOptions.mat_type=getfieldvalue(options, 'mat_type','mpiaij');
solverOptions.ksp_type=getfieldvalue(options, 'ksp_type','cgs');
solverOptions.pc_type=getfieldvalue(options, 'pc_type',  'none');


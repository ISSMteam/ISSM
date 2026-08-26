function solverOptions=bcgsmgoptions(varargin)
%BCGSNONEOPTIONS - define PETSc solver options for the BiCGSTAB (stabilized BiConjugate Gradient) Krylov method with no preconditioning
%
%   Usage:
%      solverOptions=bcgsnoneoptions(varargin);

%retrieve options provided in varargin
options=pairoptions(varargin{:});
solverOptions=struct();
solverOptions.toolkit='petsc';
solverOptions.mat_type=getfieldvalue(options, 'mat_type','mpiaij');
solverOptions.ksp_type=getfieldvalue(options, 'ksp_type','bcgs');
solverOptions.pc_type=getfieldvalue(options, 'pc_type',  'none');

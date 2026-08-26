function maltab=matlaboptions(varargin)
%MATLABOPTIONS - return PETSc options for the MATLAB built-in (backslash) direct solver
%
%   Usage:
%      maltab=matlaboptions();

%retrieve options provided in varargin
options=pairoptions(varargin{:});
maltab=struct();

%default matlab options
maltab.toolkit='petsc';
maltab.ksp_type='matlab';

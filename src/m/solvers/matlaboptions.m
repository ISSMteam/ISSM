function maltab=matlaboptions(varargin)
%MATLABOPTIONS - return Matlab petsc options
%
%   Usage:
%      options=matlaboptions;

%retrieve options provided in varargin
options=pairoptions(varargin{:});
maltab=struct();

%default matlab options
maltab.toolkit='petsc';
maltab.ksp_type='matlab';

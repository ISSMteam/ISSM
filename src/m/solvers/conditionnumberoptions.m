function cn=conditionnumberoptions(varargin)
%CONDITIONNUMBEROPTIONS - define PETSc GMRES solver options to estimate the condition number of the system matrix
%
%   Usage:
%      cn=conditionnumberoptions(varargin);

%retrieve options provided in varargin
options=pairoptions(varargin{:});
cn=struct();
cn.toolkit='petsc';
cn.mat_type=getfieldvalue(options,'mat_type','mpiaij');
cn.ksp_type=getfieldvalue(options,'ksp_type','gmres');
cn.pc_type=getfieldvalue(options,'pc_type','none');
cn.ksp_monitor_singular_value=getfieldvalue(options,'ksp_monitor_singular_value','');
cn.ksp_gmres_restart=getfieldvalue(options,'ksp_gmres_restart',1000);
cn.info=getfieldvalue(options,'info','');
cn.log_summary=getfieldvalue(options,'log_summary','');

function jacobiasm=jacobiasmoptions(varargin)
%ASMOPTIONS - return Additive Shwartz Method with Jacobi preconditioner petsc options
%
%   Usage:
%      options=jacobiasmoptions;

%retrieve options provided in varargin
options=pairoptions(varargin{:});
jacobiasm=struct();

%default jacobiasm options
jacobiasm.toolkit='petsc';
jacobiasm.mat_type=getfieldvalue(options,'mat_type','mpiaij');
jacobiasm.ksp_type=getfieldvalue(options,'ksp_type','gmres');
jacobiasm.pc_type=getfieldvalue(options,'pc_type','asm');
jacobiasm.sub_pc_type=getfieldvalue(options,'sub_pc_type','jacobi');
jacobiasm.pc_asm_overlap=getfieldvalue(options,'pc_asm_overlap',3);
jacobiasm.ksp_max_it=getfieldvalue(options,'ksp_max_it',100);
jacobiasm.ksp_rtol=getfieldvalue(options,'ksp_rtol',1e-15);

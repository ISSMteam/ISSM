function options=recover_qmu_options(md,varargin)
%RECOVER_SOLVE_OPTIONS - recover solution options for qmu runs.
%
%   Usage:
%      options=recover_qmu_options(md,varargin);
%
%   See also: SOLVE

%initialize options.
options=cell(0,2);

%make sure length(varargin) is even, ie options come in pairs.
if mod(length(varargin),2),
	error('recover_qmu_options error message: an even number of options is necessary');
end

%go through varargin, extract options 
for i=1:length(varargin)/2,

	options(end+1,:)={varargin{2*i-1} varargin{2*i}};

end

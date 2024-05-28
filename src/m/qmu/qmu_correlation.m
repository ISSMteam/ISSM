function factors=qmu_correlation(md,variablename,responsename)
%QMU_CORRELATION - compute correlation between qmu output and a certain input variable.
%
%   Usage:
%      factors=qmu_correlation(md,variablename,responsename)
%
%
%   Example:
%      mass_flux_drag_correlation=qmu_correlation(md,'drag','mass_flux');

if ~isfield(md.qmu.results,'dresp_dat'),
	error('qmu_correlation error message: could not find dresp_dat field in dakota results. you need to run montecarlo before computing correlations');
end

data=md.qmu.results.dresp_dat;

%go through all the rows and figure which one we are interested in.
found=0;
for i=1:numel(data),
	if strcmpi(data(i).descriptor,responsename),
		found=i;
		break;
	end
end
if found==0,
	error(['qmu_correlation error message: could not find data descriptor for response ' responsename]);
end

%get the response samples.
response_samples=data(found).sample;

%now go through variables, and compute correlation coefficient each time: 
variablenamelength=length(variablename);
index=[];
for i=1:numel(data),
	if strncmpi(variablename,data(i).descriptor,variablenamelength),
		%this observation is one we are looking for.
		index=[index;i];
	end
end

if isempty(index),
	error(['qmu_correlation error message: could not find correlation descriptor for variable ' variablename]);
end

factors=zeros(numel(index),1);
for i=1:numel(index),
	matrix=corrcoef(data(index(i)).sample,response_samples);
	factors(i)=matrix(2,1);
end

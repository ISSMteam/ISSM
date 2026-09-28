function dvar=setupriftsfriction(md,dvar,variables)
%SETUPRIFTSFRICTION - expand a rifts friction qmu variable into one entry per rift
%
%   Usage:
%      dvar=setupriftsfriction(md,dvar,variables)

%we have several rifts.

for j=1:md.rifts.numrifts
	dvar(end+1)           =variables;
	dvar(end  ).descriptor=sprintf('%s%d',variables.descriptor,j);
end

end

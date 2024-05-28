function dvar=setupriftsfriction(md,dvar,variables)

%we have several rifts.

for j=1:md.rifts.numrifts
	dvar(end+1)           =variables;
	dvar(end  ).descriptor=sprintf('%s%d',variables.descriptor,j);
end

end

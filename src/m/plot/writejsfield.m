function writejsfield(fid,name,variable,nods)
%WRITEJSFIELD - write variable to javascript file 
%
%   Usage:
%      writejsfield(fid,name,variable)
%

	%write array:
	if size(variable,2)==1,
		fprintf(fid,'<!-- %s{{{-->\n',name);
		fprintf(fid,'%s=[',name);
		for i=1:nods-1,
			fprintf(fid,'%g,',variable(i));
		end
		fprintf(fid,'%g];\n',variable(end));
		fprintf(fid,'<!--}}}-->\n');
	else
		%multi-sized array: 
		fprintf(fid,'<!-- %s{{{-->\n',name);
		fprintf(fid,'%s=[]\n',name);
		for i=1:size(variable,2),
			fprintf(fid,'%s["%i"]=[',name,i);
			for j=1:nods-1,
				fprintf(fid,'%g,',variable(j,i));
			end
			fprintf(fid,'%g];\n',variable(end,i));
		end
		fprintf(fid,'<!--}}}-->\n');
	end

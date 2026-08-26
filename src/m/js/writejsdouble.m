function writejsdouble(fid,prefix,scalar)
%WRITEJSDOUBLE - write a scalar double to a JavaScript file as an assignment statement
%
%   Usage:
%      writejsdouble(fid,prefix,scalar)
	if  isinf(scalar)
		fprintf(fid,'%s=Infinity;\n',prefix);
	else
		fprintf(fid,'%s=%g;\n',prefix,scalar);
	end
end

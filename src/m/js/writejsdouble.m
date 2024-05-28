function writejsdouble(fid,prefix,scalar)
	if  isinf(scalar),
		fprintf(fid,'%s=Infinity;\n',prefix);
	else
		fprintf(fid,'%s=%g;\n',prefix,scalar);
	end
end

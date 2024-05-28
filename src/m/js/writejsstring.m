function writejsstring(fid,prefix,string)
	fprintf(fid,'%s=''%s'';\n',prefix,string);
end

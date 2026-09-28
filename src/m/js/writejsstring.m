function writejsstring(fid,prefix,string)
%WRITEJSSTRING - write a string to a JavaScript file as an assignment statement
%
%   Usage:
%      writejsstring(fid,prefix,string)
	fprintf(fid,'%s=''%s'';\n',prefix,string);
end

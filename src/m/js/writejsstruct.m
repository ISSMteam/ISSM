function writejsstruct(fid,prefix,structure)
%WRITEJSSTRUCT - write a structure's scalar and string fields to a JavaScript file as object property assignments
%
%   Usage:
%      writejsstruct(fid,prefix,structure)

	fprintf(fid,'%s={};\n',prefix);

	fields=fieldnames(structure);
	for i=1:numel(fields)
		fieldname=fields{i};
		field=structure.(fieldname);
		if isscalar(field)
			fprintf(fid,'%s[''%s'']=%g;\n',prefix,fieldname,field);
		end
		if ischar(field)
			fprintf(fid,'%s[''%s'']=''%s'';\n',prefix,fieldname,field);
		end
	end
end

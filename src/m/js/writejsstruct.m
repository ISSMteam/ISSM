function writejsstruct(fid,prefix,structure)
	
	fprintf(fid,'%s={};\n',prefix);

	fields=fieldnames(structure);
	for i=1:numel(fields),
		fieldname=fields{i};
		field=structure.(fieldname);
		if isscalar(field),
			fprintf(fid,'%s[''%s'']=%g;\n',prefix,fieldname,field);
		end
		if ischar(field),
			fprintf(fid,'%s[''%s'']=''%s'';\n',prefix,fieldname,field);
		end
	end
end

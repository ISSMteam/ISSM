function writejscellarray(fid,prefix,cell)


	if ~iscell(cell),
		fprintf(fid,'%s=%g;\n',prefix,cell);
	else
		fprintf(fid,'%s=[',prefix);
		for i=1:length(cell),
			array=cell{i};
			fprintf(fid,'[');
			for j=1:length(array)-1,
				fprintf(fid,'%g,',array(j));
			end
			fprintf(fid,'%g]',array(end));
			if i<length(cell), fprintf(fid,','); end
		end
		fprintf(fid,'];\n');
	end


end

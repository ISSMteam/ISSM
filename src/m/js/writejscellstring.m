function writejscellstring(fid,prefix,cell)


	if ~iscell(cell),
		fprintf(fid,'%s=%g;\n',prefix,cell);
	else
		if length(cell),
			if length(cell)==1,
				fprintf(fid,'%s=[''%s''];\n',prefix,cell{1});
			else
				fprintf(fid,'%s=[''%s'',',prefix,cell{1});
				for i=2:length(cell)-1,
					fprintf(fid,'''%s'',',cell{i});
				end
				fprintf(fid,'''%s''];\n',cell{end});
			end
		else
			fprintf(fid,'%s=[];\n',prefix);
		end
	end

end

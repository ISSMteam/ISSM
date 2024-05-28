function writejs2Darray(fid,prefix,array)

	if  isscalar(array),
		fprintf(fid,'%s=%g;\n',prefix,array);
	else
		fprintf(fid,'%s=[',prefix);
		for i=1:size(array,1)-1,
			fprintf(fid,'[%g,',array(i,1));
			for j=2:size(array,2)-1,
				fprintf(fid,'%g,',array(i,j));
			end
			fprintf(fid,'%g],',array(i,end));
		end
		fprintf(fid,'[%g,',array(end,1));
		for j=2:size(array,2)-1,
			fprintf(fid,'%g,',array(end,j));
		end
		fprintf(fid,'%g]];\n',array(end,end));
	end
end

function writejs1Darray(fid,prefix,array)
%WRITEJS1DARRAY - write a 1-D array to a JavaScript file as an assignment statement
%
%   Usage:
%      writejs1Darray(fid,prefix,array)

	if isempty(array)
		fprintf(fid,'%s=[];\n',prefix);
	else if  isscalar(array),
		fprintf(fid,'%s=%g;\n',prefix,array);
	else
		fprintf(fid,'%s=[',prefix);
		for i=1:length(array)-1
			fprintf(fid,'%g,',array(i));
		end
		fprintf(fid,'%g];\n',array(end));
	end
end

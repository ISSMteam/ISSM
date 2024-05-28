function md=qmuname(md,varargin)
%INPUT function md=qmuname(md)
%Pick up the number from a file, or get it directly from the Dakota structure.  Then modify the name of this 
%model to reflect this new number.

if nargin==1,
	fid=fopen('number','r');
	number=fscanf(fid,'%i',1)
	fclose(fid);
else
	number=varargin{1};
end

%modify model name by appending number to the name
md.miscellaneous.name=[md.miscellaneous.name num2str(number)];

function archwrite(filename,varargin) % {{{
%ARCHWRITE - Write data to a field, given the file name, field name, and data.
%
%	Usage:
%		archwrite('archive101.arch','variable_name',data);
%		archwrite('archive102.arch','archive102_field1',transpose(field_value{1}),...
%					'archive_field2',transpose(field_value{2}),...
%					'archive_field3',-6.2420512521312);

	nvarargs=length(varargin);

	if nargin==1 || mod(nvarargs,2)~=0
		help archwrite
		error('Error: Invalid number of arguments provided. Please check the usage above.');
	end

	if ~exist(filename,'file'),
		fid=fopen(filename,'wb');	% create new file
	else
		fid=fopen(filename,'ab+');	% append to file
	end

	numfields=nvarargs/2;
	% Generate data to write
	for i=1:numfields,
		% write field name
		name=varargin{2*i-1};
		write_field_name(fid,name);

		% write data
		data=varargin{2*i};
		code=format_archive_code(data);
		if code==1
			fclose(fid);
			error('Error: Field data should be a vector or scalar, string was detected');
		elseif code==2 % scalar
			write_scalar(fid,data);
		elseif code==3 % vector
			write_vector(fid,data);
		end
	end
	% clean up
	fclose(fid);
end%}}}

%Helper functions{{{
function [code]=format_archive_code(format) % {{{
	%Given a format, return corresponding code (for reading and writing)
	if ischar(format)
		code=1; % string
	elseif isscalar(format)
		code=2; % scalar
	elseif (isvector(format) || ismatrix(format))
		code=3; % vector or matrix
	else
		error('Error! Please ensure arguments are strings, scalars, or vectors/matrixes.');
	end
end%}}}
function write_field_name(fid,field_name) % {{{
	% write record length for name
	rclen=4+4+length(field_name); % format code + length to write + amt of chars to write
	fwrite(fid,rclen,'int','ieee-be');

	% write field name to file
	fwrite(fid,1,'int','ieee-be');							% data will be a string
	fwrite(fid,length(field_name),'int','ieee-be');	% show how many chars need to be read
	fwrite(fid,field_name,'char','ieee-be');			% write field name
end%}}}
function write_scalar(fid,data) % {{{
	% write record length for scalar
	rclen=4+8; % format code + size of double (8 bytes)
	fwrite(fid,rclen,'int','ieee-be');

	% write scalar to file
	fwrite(fid,2,'int','ieee-be');			% data will be scalar
	fwrite(fid,data,'double','ieee-be');	% write scalar to file
end%}}}
function write_vector(fid,data) % {{{
	% write record length for vector
	row_size=size(data,1);
	col_size=size(data,2);
	rclen=4+4+4+8*row_size*col_size; % code+rowsz+colsz+(doublesz*rowamt*colamt)
	fwrite(fid,rclen,'int','ieee-be');
	
	% write vector to file
	fwrite(fid,3,'int','ieee-be');
	fwrite(fid,row_size,'int','ieee-be');
	fwrite(fid,col_size,'int','ieee-be');
	fwrite(fid,data,'double','ieee-be');
end%}}}
%}}}

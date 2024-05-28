function [archive_data]=archread(filename,varargin) % {{{
%ARCHREAD - Given a variable name list of variable names (given as separate arguments),
%	find and return the data associated with those variable names
%
%	Usage:
%		archive_data = archread('archive101.arch','archive101_field1');

	nvarargs=length(varargin);

	if nargin<2
		help archread
		error('Error: Must provide correct number of arguments.');
	end

	%open file
	fid=fopen(filename,'rb');
	if fid==-1
		error(['Error: File: ' filename ' does not exist.']);
	end

	% Get structure of file
	field_names={};
	field_data_positions=[]; % starting position (in bytes) of data
	archive_data=[];

	% pass 1: get field names and relate them to field data 
	while ~feof(fid)
		[reclen,count]=fread(fid,1,'int','ieee-be'); 
		if count == 0, % reached eof 
			break;
		end
		record_type=fread(fid,1,'int','ieee-be');
		if record_type~=1
			fclose(fid);
			error(['Error: Reading failed. Expected variable of type string, type given was ' format_archive_code_to_type(record_type)]);
		else
			field_name_length=fread(fid,1,'int','ieee-be');
			field_names{end+1}=char(fread(fid,field_name_length,'char','ieee-be')');
			% once we have a field name at index i, we will associate
			% the starting positiion of the data related to that field name
			% by using field data positions at index i
			field_data_positions(end+1)=ftell(fid);
			data_reclen=fread(fid,1,'int','ieee-be');
			fseek(fid,data_reclen,'cof'); % advance to next field name
		end
	end
	fseek(fid,0,-1); % rewind

	% Determine what fields we need (compare to varargin)
	idx=find(ismember(field_names,varargin));
	if size(idx) == 0,
		fclose(fid);
		error('Error: No matching variables found in archive.');
	end
	% Get data from field_data_positions
	for i=1:length(idx)
		fseek(fid,field_data_positions(idx(i)),0); % get starting position of data
		reclen=fread(fid,1,'int','ieee-be');
		field_type=fread(fid,1,'int','ieee-be');
		if field_type==2
			archive_data{i}=fread(fid,1,'double','ieee-be');
		elseif field_type==3
			archive_data{i}=read_vector(fid);
		else
			fclose(fid);
			error('Error: Encountered invalid field type when reading data.');
		end
	end
	fclose(fid);
end%}}}

%Helper functions{{{
function [code]=format_archive_code(format) % {{{
%Given a format, return corresponding code (for reading and writing)
	if ischar(format)
		code=1; % string
	elseif isscalar(format)
		code=2; % scalar
	elseif isvector(format) or ismatrix(format)
		code=3; % vector 
	else
		error('Error! Please ensure arguments are strings, scalars, or vectors.');
	end
end%}}}
function [str]=format_archive_code_to_type(code) % {{{
	if code == 1
		str='string';
	elseif code == 2
		str='scalar';
	elseif code == 3
		str='vector';
	else
		error(['Code ' num2str(code) ' is not associated with any known format.']);
	end
end%}}}
function [data]=read_vector(fid) % {{{
	rows=fread(fid,1,'int','ieee-be'); 
	cols=fread(fid,1,'int','ieee-be');
	data=fread(fid,[rows,cols],'double','ieee-be');
end%}}}
%}}}

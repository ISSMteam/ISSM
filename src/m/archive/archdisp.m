function archdisp(filename) % {{{
%ARCHDISP - Display the contents of a .arch file
%
%	Usage:
%		archdisp('archive101.arch');

	if nargin~=1
		help archdisplay
		error('Error: Invalid number of arguments, only one file can be displayed at a time');
	end

	if ~exist(filename,'file'),
		error('Error: File does not exist!');
	end
	
	fid=fopen(filename,'rb');
	fprintf('Source file: \n');
	fprintf('\t%s\n',filename);
	fprintf('Variables: \n');

	field_names={};
	field_sizes={};
	field_types={};
	archive_data=[];

	% Read file
	while ~feof(fid)
		[reclen,count]=fread(fid,1,'int','ieee-be');
		if count == 0, % reached eof
			break;
		end
		% Read variable name
		fread(fid,1,'int','ieee-be'); % this will always return a string
		field_name_length=fread(fid,1,'int','ieee-be');
		field_names{end+1}=char(fread(fid,field_name_length,'char','ieee-be')');
		% Read variable traits
		data_reclen=fread(fid,1,'int','ieee-be');
		field_type=fread(fid,1,'int','ieee-be');
		% return related type
		if field_type==2
			field_types{end+1}='double';
			field_sizes{end+1}='1x1'; % a scalar will always be a 1x1 cell array
			archive_data{end+1}=fread(fid,1,'double','ieee-be');
		elseif field_type==3
			field_types{end+1}='vector/matrix';
			rows=fread(fid,1,'int','ieee-be');
			cols=fread(fid,1,'int','ieee-be');
			data=fread(fid,[rows,cols],'double','ieee-be');
			field_size_cat=strcat(num2str(rows),'x',num2str(cols)); 
			field_sizes{end+1}=field_size_cat;
			archive_data{end+1}=data;
		else
			fclose(fid);
			error('Encountered invalid data type while reading archive.');
		end
	end
	fclose(fid);

	% Display contents
	for i=1:numel(field_names)
		fprintf('\t%s\n',field_names{i});
		fprintf('\t\tSize:\t\t%s\n',field_sizes{i});
		fprintf('\t\tDatatype:\t%s\n',field_types{i});	
	end
end % }}}

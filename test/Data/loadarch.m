function [s]=loadarch(archfile);
%LOADARCH - Function to read the variables from an arch-formatted file
%
%	Usage:
%		s=loadarch(archfile);
%	where:
%		archfile			name of an arch-format file
	
	if ~exist('archfile','var') || isempty(archfile)
		help loadarch
		error('loadarch usage error.');
	end

	[pathstr,name,ext]=fileparts(archfile);
	if isempty(ext)
		ext='.arch';
	end
	archfile=fullfile(pathstr,[name ext]);

	% Fields for our structure 's'
	variable_names={};
	variable_sizes={};
	variable_types={};
	variable_data={};

	% Read file and gather variable names and sizes
	fid=fopen(archfile,'rb');
	while ~feof(fid)
		[reclen,count]=fread(fid,1,'int','ieee-be');
		if count == 0, % reached eof
			break;
		end
		% Read variable traits
		fread(fid,1,'int','ieee-be');
		variable_name_len=fread(fid,1,'int','ieee-be');
		variable_names{end+1}=char(fread(fid,variable_name_len,'char','ieee-be')');
		variable_reclen=fread(fid,1,'int','ieee-be');
		variable_type=fread(fid,1,'int','ieee-be');
		% set related data type
		if variable_type==2
			variable_types{end+1}='double';
			variable_sizes{end+1}='1x1';
			variable_data{end+1}=fread(fid,1,'double','ieee-be'); % read data
		elseif variable_type==3
			variable_types{end+1}='double vector / matrix';
			rows=fread(fid,1,'int','ieee-be');
			cols=fread(fid,1,'int','ieee-be');
			variable_data{end+1}=fread(fid,[rows,cols],'double','ieee-be');
			variable_size_cat=strcat(num2str(rows),'x',num2str(cols));
			variable_sizes{end+1}=variable_size_cat;
		else
			fclose(fid);
			error('Error: Encountered invalid data type while loading arch information.');
		end
	end
	fclose(fid);
	disp(sprintf('arch-format file ''%s'' read.',archfile));	

	% Relate the variable fields to the structure
	% Access fields by doing (for example): s(1).Name -> 'VelocityXObservation'
	s=struct('FileName',archfile,'Name',variable_names,'Size',variable_sizes,'Type',variable_types,'Data',variable_data);
	% Let user know data has been successfully read
	for i=1:numel(variable_names)
		disp(sprintf('field ''%s'' of class ''%s'' and size [%s] read.',...
			variable_names{i},variable_types{i},variable_sizes{i}));
	end
end

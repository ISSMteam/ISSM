function WriteData(fid,prefix,varargin)
%WRITEDATA - write model field in binary file
%
%   Usage:
%      WriteData(fid,varargin);

%process options
options=pairoptions(varargin{:});

%Get data properties
if exist(options,'object')
	obj       = getfieldvalue(options,'object');
	fieldname = getfieldvalue(options,'fieldname');
	classname = getfieldvalue(options,'class',class(obj));
	name      = getfieldvalue(options,'name',[prefix '.' fieldname ]);
	if exist(options,'data'),
		data = getfieldvalue(options,'data');
	else
		data = obj.(fieldname);
	end
else
	data = getfieldvalue(options,'data');
	name = getfieldvalue(options,'name');
end
format  = getfieldvalue(options,'format');
mattype = getfieldvalue(options,'mattype',0);    %only required for matrices
timeserieslength = getfieldvalue(options,'timeserieslength',-1);

%Process sparse matrices
if issparse(data)
	data=full(data);
end

%Scale data if necessary
if strcmpi(format,'MatArray')
	for i=1:numel(data)
		if exist(options,'scale'),
			scale = getfieldvalue(options,'scale');
			if size(data{i},1)==timeserieslength,
				data{i}(1:end-1,:) = scale.*data{i}(1:end-1,:);
			else
				data{i} = scale.*data{i};
			end
		end
		if size(data{i},1)==timeserieslength,
			yts = getfieldvalue(options,'yts');
			data{i}(end,:) = data{i}(end,:)*yts;
		end
	end
else
	if exist(options,'scale')
		scale = getfieldvalue(options,'scale');
		if size(data,1)==timeserieslength,
			data(1:end-1,:) = scale.*data(1:end-1,:);
		else
			data  = scale.*data;
		end
	end
	if(size(data,1)==timeserieslength)
		yts = getfieldvalue(options,'yts');
		data(end,:) = data(end,:)*yts;
	end
end

%Step 1: write the name to identify this record uniquely
fwrite(fid,numel(name),'int');
fwrite(fid,name,'char');

%Step 2: write the data itself.
if     strcmpi(format,'Boolean'),% {{{
	if(numel(data)~=1), error(['field ' name ' cannot be marshalled as it has more than one element!']); end

	%first write length of record
	fwrite(fid,4+4,'int64');  %1 bool (disguised as an int)+code

	%write data code:
	fwrite(fid,FormatToCode(format),'int');

	%now write integer
	fwrite(fid,data,'int');  %send an int, not easy to send a bool
	% }}}
elseif strcmpi(format,'Integer'), % {{{
	if(numel(data)~=1), error(['field ' name ' cannot be marshalled as it has more than one element!']); end

	%first write length of record
	fwrite(fid,4+4,'int64');  %1 integer + code

	%write data code:
	fwrite(fid,FormatToCode(format),'int');

	%now write integer
	fwrite(fid,data,'int');
	% }}}
elseif strcmpi(format,'Double'), % {{{
	if(numel(data)~=1), error(['field ' name ' cannot be marshalled as it has more than one element!']); end

	%first write length of record
	fwrite(fid,4+8,'int64');  %1 double+code

	%write data code:
	fwrite(fid,FormatToCode(format),'int');

	%now write double
	fwrite(fid,data,'double');
	% }}}
elseif strcmpi(format,'String'), % {{{
	%first write length of record
	fwrite(fid,length(data)+4+4,'int64');  %string + string size + code

	%write data code:
	fwrite(fid,FormatToCode(format),'int');

	%now write string
	fwrite(fid,length(data),'int');
	fwrite(fid,data,'char');
	% }}}
elseif strcmpi(format,'BooleanMat'), % {{{

	%Get size
	s=size(data);
	%if matrix = NaN, then do not write anything
	if (s(1)==1 & s(2)==1 & isnan(data)),
		s(1)=0; s(2)=0;
	end

	%first write length of record
	fwrite(fid,4+4+8*s(1)*s(2)+4+4,'int64');  %2 integers (32 bits) + the double matrix + code + matrix type

	%write data code and matrix type:
	fwrite(fid,FormatToCode(format),'int');
	fwrite(fid,mattype,'int');

	%now write matrix
	fwrite(fid,s(1),'int');
	fwrite(fid,s(2),'int');
	if s(1)*s(2),
		fwrite(fid,data','double'); %get to the "c" convention, hence the transpose
	end
	% }}}
elseif strcmpi(format,'IntMat'), % {{{

	%Get size
	s=size(data);
	%if matrix = NaN, then do not write anything
	if (s(1)==1 & s(2)==1 & isnan(data)),
		s(1)=0; s(2)=0;
	end

	%first write length of record
	fwrite(fid,4+4+8*s(1)*s(2)+4+4,'int64');  %2 integers (32 bits) + the double matrix + code + matrix type

	%write data code and matrix type:
	fwrite(fid,FormatToCode(format),'int');
	fwrite(fid,mattype,'int');

	%now write matrix
	fwrite(fid,s(1),'int');
	fwrite(fid,s(2),'int');
	if s(1)*s(2),
		fwrite(fid,data','double'); %get to the "c" convention, hence the transpose
	end
	% }}}
elseif strcmpi(format,'DoubleMat'), % {{{

	%Get size
	s=size(data);

	if numel(s)~=2
		error('matrices that that have more than 2 dimensions are not supported');
	end

	%if matrix = NaN, then do not write anything
	if (s(1)==1 & s(2)==1 & isnan(data)),
		s(1)=0; s(2)=0;
	end

	%first write length of record
	recordlength=4+4+8*s(1)*s(2)+4+4; %2 integers (32 bits) + the double matrix + code + matrix type
	if recordlength>2^63; error(['field ' name ' cannot be marshalled because it is larger than 2^63 bytes!']); end
	fwrite(fid,recordlength,'int64');

	%write data code and matrix type:
	fwrite(fid,FormatToCode(format),'int');
	fwrite(fid,mattype,'int');

	%now write matrix
	fwrite(fid,s(1),'int');
	fwrite(fid,s(2),'int');
	if s(1)*s(2),
		fwrite(fid,data','double'); %get to the "c" convention, hence the transpose
	end
	% }}}
elseif strcmpi(format,'CompressedMat'), % {{{

	%Get size
	s=size(data);

	if (s(1)==1 & s(2)==1 & isnan(data)),
		s(1)=0; s(2)=0;
	end

	%first write length of record
	recordlength=4+4+8+8+1*(s(1)-1)*s(2)+8*s(2)+4+4; %2 integers (32 bits) + the matrix + code + matrix type
	if recordlength>2^63; error(['field ' name ' cannot be marshalled because it is larger than 2^63 bytes!']); end
	fwrite(fid,recordlength,'int64');

	%write data code and matrix type:
	fwrite(fid,FormatToCode(format),'int');
	fwrite(fid,mattype,'int');

	%write matrix size
	fwrite(fid,s(1),'int');
	fwrite(fid,s(2),'int');

	if s(1)*s(2),

		%Write offset and range
		A = data(1:end-1,:);
		offset = min(A(:));
		range = max(A(:)) - offset;
		fwrite(fid,offset,'double');
		fwrite(fid,range,'double');

		%Convert data to uint8 and write it
		A=uint8((A-offset)/range*255);
		fwrite(fid,A','uint8'); %get to the "c" convention, hence the transpose

		%Write last row as double (time)
		fwrite(fid,data(end,:),'double');
	else

		%Write empty offset and range
		fwrite(fid,0,'double');
		fwrite(fid,0,'double');
	end
	% }}}
elseif strcmpi(format,'MatArray'), % {{{

	numrecords=numel(data);

	%first get length of record
	recordlength=4+4; %number of records + code
	for i=1:numrecords,
		matrix=data{i};
		s=size(matrix);
		recordlength=recordlength+4*2+... %row and col of matrix
			s(1)*s(2)*8; %matrix of doubles
	end

	%write length of record
	fwrite(fid,recordlength,'int64');

	%write data code:
	fwrite(fid,FormatToCode(format),'int');

	%write data, first number of records
	fwrite(fid,numrecords,'int');

	%write each matrix:
	for i=1:numrecords,
		matrix=data{i};
		s=size(matrix);
		fwrite(fid,s(1),'int');
		fwrite(fid,s(2),'int');
		fwrite(fid,matrix','double');
	end
	% }}}
elseif strcmpi(format,'StringArray'), % {{{

	%first get length of string array:
	num=numel(data);
	if isnumeric(data) & num==1 & isnan(data),
		num = 0;
	end

	%now get length of record:
	recordlength=4+4; %for length of array + code
	for i=1:num,
		string=data{i};
		recordlength=recordlength+4+length(string); %for each string
	end

	%write length of record
	fwrite(fid,recordlength,'int64');

	%write data code:
	fwrite(fid,FormatToCode(format),'int');

	%now write length of string array
	fwrite(fid,num,'int');

	%now write the strings
	for i=1:num,
		string=data{i};
		fwrite(fid,length(string),'int');
		fwrite(fid,string,'char');
	end
	% }}}
else  % {{{
	error(['WriteData error message: data type: ' num2str(format) ' not supported yet! (' name ')']);
end % }}}
end

function code=FormatToCode(format) % {{{
%This routine takes the format string, and hardcodes it into an integer, which
%is passed along the record, in order to identify the nature of the dataset being
%sent.
	if     strcmpi(format,'Boolean'),
		code=1;
	elseif strcmpi(format,'Integer'),
		code=2;
	elseif strcmpi(format,'Double'),
		code=3;
	elseif strcmpi(format,'String'),
		code=4;
	elseif strcmpi(format,'BooleanMat'),
		code=5;
	elseif strcmpi(format,'IntMat'),
		code=6;
	elseif strcmpi(format,'DoubleMat'),
		code=7;
	elseif strcmpi(format,'MatArray'),
		code=8;
	elseif strcmpi(format,'StringArray'),
		code=9;
	elseif strcmpi(format,'CompressedMat'),
		code=10;
	else
		error('FormatToCode error message: data type not supported yet!');
	end
end% }}}

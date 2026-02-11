function marshall(md)
%MARSHALL - outputs a compatible binary file from @model md, for certain solution type.
%
%   The routine creates a compatible binary file from @model md
%   This binary file will be used for parallel runs in JPL-package
%
%   Usage:
%      marshall(md)

if md.verbose.solution
	disp(['marshalling file ' md.miscellaneous.name '.bin']);
end

%open file for binary writing
fid=fopen([ md.miscellaneous.name '.bin'],'wb');
if fid==-1
	error(['marshall error message: could not open ' [md.miscellaneous.name '.bin'],' file for binary writing']);
end

% Go through all model fields: check that it is a class and call checkconsistency
fields=sort(properties('model')); %sort fields so that comparison of binary files is easier
for i=1:length(fields)
	field=fields{i};

	%Some properties do not need to be marshalled
	if ismember(field,{'results' 'radaroverlay' 'toolkits' 'cluster' 'private'}),
		continue;
	end

	%Check that current field is an object
	if ~isobject(md.(field))
		error(['field ''' char(field) ''' is not an object']);
	end

	%Marshall current object
	%disp(['marshalling ' field '...']); %Uncomment for debugging
	marshall(md.(field),['md.' field],md,fid);
end

%Last, write "md.EOF" to make sure that the binary file is not corrupt
WriteData(fid,'XXX','name','md.EOF','data',true,'format','Boolean');

%close file
st=fclose(fid);

if st==-1,
	error(['marshall error message: could not close file ' [md.miscellaneous.name '.bin']]);
end

% Uncomment the following to make a copy of the binary input file for debugging 
% purposes (can be fed into scripts/BinRead.py).
% copyfile([md.miscellaneous.name '.bin'], [md.miscellaneous.name '.m.bin'])

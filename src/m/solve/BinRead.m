function md = BinReadMatlab(infile)
%BINREADMATLAB - Read ISSM input binary 
%
%   struct = BinReadMatlab(infile)

	% Open binary file for reading
	fid = fopen(infile, 'rb');
	if fid == -1
		error('Cannot open input file %s', infile);
	end

	%Prepare output
	md = model();

	while(true)
		% Step 1: read size of record name (int32)
		[recordNameSize, cnt] = fread(fid, 1, 'int32');
		if cnt ~= 1
			% EOF reached
			break;
		end

		% Read the record name 
		rawName = fread(fid, recordNameSize, '*char')';
		recordName = char(rawName);
		fprintf('Reading %s\n', recordName);

		%See if field exists in md
		parts =  strsplit(recordName, '.');
		addtomd = true;
		if numel(parts)<3
			addtomd = false;
		elseif ~isprop(md, parts{2})
			addtomd = false;
		elseif ~isprop(md.(parts{2}), parts{3})
			addtomd = false;
		end


		% Step 2: read length of record (int64)
		reclen = fread(fid, 1, 'int64');

		% Read data code (int32)
		code = fread(fid, 1, 'int32');
		formatStr = CodeToFormat(code);
		%fprintf('Format = %s (code %d)\n', formatStr, code);

		% Dispatch based on the data type
		switch code
			case FormatToCode('Boolean')
				% Boolean stored as 32‑bit integer (0/1)
				data = fread(fid, 1, 'int32');

			case FormatToCode('Integer')
				data = fread(fid, 1, 'int32');

			case FormatToCode('Double')
				data = fread(fid, 1, 'double');

			case FormatToCode('String')
				len  = fread(fid, 1, 'int');
				data = char(fread(fid, len, 'char'))';

			case {FormatToCode('BooleanMat'), FormatToCode('IntMat'), FormatToCode('DoubleMat')}

				% Read matrix type (int32) – not used further here
				mattype = fread(fid, 1, 'int32');

				% Read dimensions (two int32s)
				rows = fread(fid, 1, 'int32');
				cols = fread(fid, 1, 'int32');

				% Choose element precision according to the specific matrix type
				switch code
					case FormatToCode('BooleanMat')
						elemType = 'double';
					case FormatToCode('IntMat')
						elemType = 'int32';
					case FormatToCode('DoubleMat')
						elemType = 'double';
				end

				% Read the matrix data column‑major (Python wrote row‑major;
				% MATLAB stores column‑major, so we transpose after reading)
				data = fread(fid, rows*cols, elemType);
				data = reshape(data, [cols, rows])';


			otherwise
				% For unsupported or array types we simply skip the remaining bytes
				% (same behaviour as the Python script)
				skipBytes = reclen - 4;  % 4 bytes already consumed by the code
				fseek(fid, skipBytes, 'cof');
				fprintf('skipping %d bytes for code %d\n', skipBytes, code);
		end
		if addtomd
			md.(parts{2}).(parts{3}) = data;
		else
			disp([' -> Skipping ', recordName]);
		end

		%Change class depending on what's provided
		if strcmp(recordName, 'md.mesh.elementtype')
			switch(data)
				case 'Tria';  md.mesh = mesh2d();
				case 'Penta'; md.mesh = mesh3dprisms();
				otherwise error(['Element type ', data, ' not supported yet']);
			end
		end
	end
	fclose(fid);
end

function fmt = CodeToFormat(code) % {{{
    switch code
        case 1, fmt = 'Boolean';
        case 2, fmt = 'Integer';
        case 3, fmt = 'Double';
        case 4, fmt = 'String';
        case 5, fmt = 'BooleanMat';
        case 6, fmt = 'IntMat';
        case 7, fmt = 'DoubleMat';
        case 8, fmt = 'MatArray';
        case 9, fmt = 'StringArray';
        otherwise, error('Unsupported code: %d', code);
    end
end % }}}
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

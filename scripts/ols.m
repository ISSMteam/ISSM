function ols(varargin)

	options=pairoptions(varargin{:});

	range=getfieldvalue(options,'<',Inf);
	file=getfieldvalue(options,'file','runme.m');

	%Open runme.m file and read line by line
	fid=fopen(file,'r');

	tline = fgets(fid);
	count=1;
	strings={};
	while ischar(tline)
		tline = fgets(fid);
		if strncmpi(tline,'if perform(org,',14),
			lastchar = strfind(tline,')');
			lastchar = lastchar(end)-1;
			string=tline(17:lastchar-1);
			if strcmpi(string,'End'),
				break;
			end
			strings{end+1}= sprintf('%2i: %s',count,tline(16:lastchar));
			if count>range,
				break;
			end
			count=count+1;
		end
	end
	fclose(fid);


	%ok, we have our strings, split in two: 
	nstring=length(strings);
	split=floor(nstring/2)

	for i=1:split,
		disp(sprintf('%s     %s',pad(strings{i},60),strings{split+i}));
	end
	if nstring==(2*split),
		return;
	else
		disp(sprintf('%s     %s',pad('',60),strings{end}));
	end

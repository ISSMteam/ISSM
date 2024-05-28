function ol(varargin)

	options=pairoptions(varargin{:});

	%recover steps in calling workspace:  will be used to highlight current step.
	steps= evalin('base', 'steps');
	
	mmin=getfieldvalue(options,'>',1);
	mmax=getfieldvalue(options,'<',Inf);

	file=getfieldvalue(options,'file','runme.m');

	%Open runme.m file and read line by line
	fid=fopen(file,'r');

	tline = fgets(fid);
	count=1;
	while ischar(tline)
		tline = fgets(fid);
		if strncmpi(tline,'if perform(org,',14),
			lastchar = strfind(tline,')');
			lastchar = lastchar(end)-1;
			string=tline(17:lastchar-1);
			if strcmpi(string,'End'),
				return;
			end
			if ismember(count,steps),
				disp(sprintf('%2i: *%s',count,tline(17:lastchar-1)));
			else
				if count>=mmin & count <= mmax,
					disp(sprintf('%2i:  %s',count,tline(17:lastchar-1)));
				else
					%do nothing. 
				end
			end
			if count>mmax,
				break;
			end
			count=count+1;
		end
	end
	fclose(fid);

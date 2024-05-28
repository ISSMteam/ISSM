function results = outbinread(md,filename,steps),
%OUTBINREAD - read specific steps in outbin file
%
%   Usage:
%      results = outbinread(md,filename,steps)
%
%   Example:
%      results = outbinread(md,'dir/test.outbin',1:10)

if ~exist(filename,'file'),
	error(['Binary file ' filename ' not found']);
end

%Initialize output
results=struct();

%Open file
fid=fopen(filename,'rb');
if(fid==-1),
	error(['could not open ',filename,' for binary reading']);
end

%Read fields until the end of the file.
result  = ReadData(fid,md);
if isempty(result), error(['no results found in binary file ' filename]); end
check_nomoresteps=0;
counter = 1;
step    = result.step;
while ~isempty(result), 

	if check_nomoresteps,
		%check that the new result does not add a step, which would be an error: 
		if result.step>=1,
			error('parsing results for a steady-state core, which incorporates transient results!');
		end
	end

	%Check step, increase counter if this is a new step
	if(step~=result.step & result.step>1)
		counter = counter + 1;
		step    = result.step;
	end

	%Add result
	if(result.step==0),
		%if we have a step = 0, this is a steady state solution, don't expect more steps. 
		index = 1;
		check_nomoresteps=1;
	elseif(result.step==1),
		index = 1;
	else
		index = counter;
	end
	pos=find(result.step==steps);
	if ~isempty(pos),
		results(pos).(result.fieldname)=result.field;
		if(result.step~=-9999),
			results(pos).step=result.step;
		end
		if(result.time~=-9999),
			results(pos).time=result.time;
		end
	end
	if result.step>max(steps(:)),
		%we are done!
		return;
	end

	%read next result
	try,
		result  = ReadData(fid,md);
	catch me,
		disp('WARNING: file corrupted, trying partial recovery');
		result=[];
	end

end

fclose(fid);

function result=ReadData(fid,md) % {{{

%read field
[length,count]=fread(fid,1,'int');

if count==0,
	result=struct([]);
else
	fieldname=fread(fid,length,'char');
	fieldname=fieldname(1:end-1)';
	fieldname=char(fieldname);
	time=fread(fid,1,'double');
	step=fread(fid,1,'int');

	type=fread(fid,1,'int');
	M=fread(fid,1,'int');
	if type==1,
		field=fread(fid,M,'double');
	elseif type==2,
		field=fread(fid,M,'char');
		field=char(field(1:end-1)');
	elseif type==3,
		N=fread(fid,1,'int');
		field=transpose(fread(fid,[N M],'double'));
	else
		error(['cannot read data of type ' num2str(type) ]);
	end

	%Process units here FIXME: this should not be done here!
	yts=md.constants.yts;
	if strcmp(fieldname,'BalancethicknessThickeningRate'),
		field = field*yts;
	elseif strcmp(fieldname,'HydrologyWaterVx'),
		field = field*yts;
	elseif strcmp(fieldname,'HydrologyWaterVy'),
		field = field*yts;
	elseif strcmp(fieldname,'Vx'),
		field = field*yts;
	elseif strcmp(fieldname,'Vy'),
		field = field*yts;
	elseif strcmp(fieldname,'Vz'),
		field = field*yts;
	elseif strcmp(fieldname,'Vel'),
		field = field*yts;
	elseif strcmp(fieldname,'BasalforcingsGroundediceMeltingRate'),
		field = field*yts;
	elseif strcmp(fieldname,'BasalforcingsFloatingiceMeltingRate'),
		field = field*yts;
	elseif strcmp(fieldname,'TotalFloatingBmb'),
		field = field/10.^12*yts; %(GigaTon/year)
	elseif strcmp(fieldname,'TotalGroundedBmb'),
		field = field/10.^12*yts; %(GigaTon/year)
	elseif strcmp(fieldname,'TotalSmb'),
		field = field/10.^12*yts; %(GigaTon/year)
	elseif strcmp(fieldname,'TotalMelt'),
		field = field/10.^12*yts; %(GigaTon/year)
	elseif strcmp(fieldname,'TotalRefreeze'),
		field = field/10.^12*yts; %(GigaTon/year)
	elseif strcmp(fieldname,'SmbMassBalance'),
		field = field*yts;
	elseif strcmp(fieldname,'SmbPrecipitation'),
		field = field*yts;
	elseif strcmp(fieldname,'SmbRunoff'),
		field = field*yts;
	elseif strcmp(fieldname,'SmbCondensation'),
		field = field*yts;
	elseif strcmp(fieldname,'SmbAccumulation'),
		field = field*yts;
	elseif strcmp(fieldname,'SmbMelt'),
		field = field*yts;
	elseif strcmp(fieldname,'CalvingCalvingrate'),
		field = field*yts;
	end

	result.fieldname=fieldname;
	result.time=time;
	if result.time~=-9999,
		result.time=time/yts;
	end
	result.step=step;
	result.field=field;

end
% }}}

function results=parseresultsfromdisk(md,filename,iosplit) % {{{

if iosplit,
	results=parseresultsfromdiskiosplit(md,filename);
else
	results=parseresultsfromdiskioserial(md,filename);
	%results=parseresultsfromdiskioserialsequential(md,filename);
end
% }}}
function results=parseresultsfromdiskiosplit(md,filename) % {{{

%Open file
fid=fopen(filename,'rb');
if(fid==-1),
	error(['loadresultsfromdisk error message: could not open ',filename,' for binary reading']);
end
results=struct();

%if we have done split I/O, ie, we have results that are fragmented across patches,
%do a first pass, and figure out the structure of results
result=ReadDataDimensions(fid);
while ~isempty(result),

	%Get time and step
	results(result.step).step=result.step;
	if result.time~=-9999,
		results(result.step).time=result.time;
	end

	%Add result
	results(result.step).(result.fieldname)=NaN;

	%read next result
	result=ReadDataDimensions(fid);
end

%do a second pass, and figure out the size of the patches
fseek(fid,0,-1); %rewind
result=ReadDataDimensions(fid);
while ~isempty(result),
	%read next result
	result=ReadDataDimensions(fid);
end

%third pass, this time to read the real information
fseek(fid,0,-1); %rewind
result=ReadData(fid,md);
while ~isempty(result),

	%Get time and step
	results(result.step).step=result.step;
	if result.time~=-9999,
		results(result.step).time=result.time;
	end

	%Add result
	results(result.step).(result.fieldname)=result.field;

	%read next result
	try,
		result=ReadData(fid,md);
	catch me,
		disp('WARNING: file corrupted, results partial recovery');
		result=[];
	end

end

%close file
fclose(fid);
% }}}
function results=parseresultsfromdiskioserial(md,filename) % {{{

%Open file
fid=fopen(filename,'rb');
if(fid==-1),
	error(['loadresultsfromdisk error message: could not open ',filename,' for binary reading']);
end

%Collect all results in a cell array
allresults = {};
counter    = 1;
while(true)

	%read next result
	try,
		result = ReadData(fid,md);
	catch me,
		disp('WARNING: file corrupted, trying partial recovery');
		continue;
	end

	%Have we reached the end of the file?
	if isempty(result),
		if counter==1
			error(['no results found in binary file ' filename]);
		else
			break;
		end
	end

	%Add result to cell array
	allresults{counter} = result;
	counter = counter+1;
end
fclose(fid);

%Now, process all results and find out how many steps we have
numresults = numel(allresults);
allsteps   = zeros(numresults,1);
for i=1:numresults
	allsteps(i) = allresults{i}.step;
end
pos = find(allsteps~=-9999);
allsteps = sort(unique(allsteps(pos)));

%Ok, now construct structure
results=struct();
for i=1:numresults
	result = allresults{i};
	index = 1;
	if(result.step ~= -9999)
		index = find(result.step == allsteps);
		results(index).step=result.step;
	end
	if(result.time~=-9999)
		results(index).time=result.time;
	end
	results(index).(result.fieldname)=result.field;
end
% }}}
function results=parseresultsfromdiskioserialsequential(md,filename) % {{{

%Open file
fid=fopen(filename,'rb');
if(fid==-1),
	error(['loadresultsfromdisk error message: could not open ',filename,' for binary reading']);
end

%first pass to figure out the steps we have: 
steps=[];
while  1, 
	result  = ReadDataDimensions(fid);
	if isempty(result),
		break;
	end
	if result.step~=-9999,
		steps=[steps result.step];
	end
end

steps=unique(steps);

%create structure: 
results=struct('step',num2cell(steps));

%second pass to fill the steps we have: 
fseek(fid,0,-1); %rewind
while  1, 
	result  = ReadData(fid,md);
	if isempty(result),
		break;
	end
	if result.step==-9999,
		result.step=1;
		result.time=0;
	end
	%find where we are putting this step: 
	ind=find(steps==result.step);
	if isempty(ind),
		error('could not find a step in our pre-structure! Something is very off!');
	end

	%plug: 
	results(ind).time=result.time;
	results(ind).(result.fieldname)=result.field;
end
fclose(fid);
%}}}
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
		field=fread(fid,[N M],'double')';
	elseif type==4,
		N=fread(fid,1,'int');
		field=fread(fid,[N M],'int')';
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
	elseif strcmp(fieldname,'VxShear'),
		field = field*yts;
	elseif strcmp(fieldname,'VyShear'),
		field = field*yts;
	elseif strcmp(fieldname,'VxBase'),
		field = field*yts;
	elseif strcmp(fieldname,'VyBase'),
		field = field*yts;
	elseif strcmp(fieldname,'VxSurface'),
		field = field*yts;
	elseif strcmp(fieldname,'VySurface'),
		field = field*yts;
	elseif strcmp(fieldname,'VxAverage'),
		field = field*yts;
	elseif strcmp(fieldname,'VyAverage'),
		field = field*yts;
	elseif strcmp(fieldname,'VxDebris'),
		field = field*yts;
	elseif strcmp(fieldname,'VyDebris'),
		field = field*yts;
	elseif strcmp(fieldname,'BasalforcingsGroundediceMeltingRate'),
		field = field*yts;
	elseif strcmp(fieldname,'BasalforcingsFloatingiceMeltingRate'),
		field = field*yts;
	elseif strcmp(fieldname,'BasalforcingsSpatialDeepwaterMeltingRate'),
		field = field*yts;
	elseif strcmp(fieldname,'BasalforcingsSpatialUpperwaterMeltingRate'),
		field = field*yts;
	elseif strcmp(fieldname,'TotalFloatingBmb'),
		field = field/10.^12*yts; %(GigaTon/year)
	elseif strcmp(fieldname,'TotalFloatingBmbScaled'),
		field = field/10.^12*yts; %(GigaTon/year)
	elseif strcmp(fieldname,'TotalGroundedBmb'),
		field = field/10.^12*yts; %(GigaTon/year)
	elseif strcmp(fieldname,'TotalGroundedBmbScaled'),
		field = field/10.^12*yts; %(GigaTon/year)
	elseif strcmp(fieldname,'TotalSmb'),
		field = field/10.^12*yts; %(GigaTon/year)
	elseif strcmp(fieldname,'TotalSmbScaled'),
		field = field/10.^12*yts; %(GigaTon/year)
	elseif strcmp(fieldname,'TotalSmbMelt'),
		field = field/10.^12*yts; %(GigaTon/year)
	elseif strcmp(fieldname,'TotalSmbRefreeze'),
		field = field/10.^12*yts; %(GigaTon/year)
	elseif strcmp(fieldname,'GroundinglineMassFlux'),
		field = field/10.^12*yts; %(GigaTon/year)
	elseif strcmp(fieldname,'IcefrontMassFlux'),
		field = field/10.^12*yts; %(GigaTon/year)
	elseif strcmp(fieldname,'IcefrontMassFluxLevelset'),
		field = field/10.^12*yts; %(GigaTon/year)
	elseif strcmp(fieldname,'SmbMassBalance'),
		field = field*yts;
	elseif strcmp(fieldname,'SmbPrecipitation'),
		field = field*yts;
	elseif strcmp(fieldname,'SmbRain'),
		field = field*yts;
	elseif strcmp(fieldname,'SmbRunoff'),
		field = field*yts;
	elseif strcmp(fieldname,'SmbRunoffSubstep'),
		field = field*yts;
	elseif strcmp(fieldname,'SmbEvaporation'),
		field = field*yts;
	elseif strcmp(fieldname,'SmbRefreeze'),
		field = field*yts;
	elseif strcmp(fieldname,'SmbEC'),
		field = field*yts;
	elseif strcmp(fieldname,'SmbAccumulation'),
		field = field*yts;
	elseif strcmp(fieldname,'SmbMelt'),
		field = field*yts;
	elseif strcmp(fieldname,'SmbMAdd'),
		field = field*yts;
	elseif strcmp(fieldname,'SmbWAdd'),
		field = field*yts;
	elseif strcmp(fieldname,'CalvingCalvingrate'),
		field = field*yts;
	elseif strcmp(fieldname,'Calvingratex'),
		field = field*yts;
	elseif strcmp(fieldname,'Calvingratey'),
		field = field*yts;
	elseif strcmp(fieldname,'CalvingMeltingrate'),
		field = field*yts;
	elseif (strcmp(fieldname,'LoveKernelsReal') | strcmp(fieldname,'LoveKernelsImag')),
		nlayer = md.materials.numlayers; 
		degmax = md.love.sh_nmax; 
		nfreq  = md.love.nfreq; 
		r0 = md.love.r0; 
		g0 = md.love.g0; 
		mu0 = md.love.mu0;
		rr = md.materials.radius; 
		rho= md.materials.density;
		rho_avg = sum(rho.*diff(rr.^3)/sum(diff(rr.^3)));
		temp_field = cell(degmax+1,nfreq,nlayer+1,6);
		for ii=1:degmax+1
			for jj=1:nfreq
				for kk=1:nlayer+1
					if (kk<nlayer+1)
						ll = (ii-1)*(nlayer+1)*6 + ((kk-1)*6+1) + 3;
						temp_field{ii,jj,kk,1} = field(ll+(1-1),jj)*r0;		% mm = 4
						temp_field{ii,jj,kk,2} = field(ll+(2-1),jj)*mu0;	% mm = 5 
						temp_field{ii,jj,kk,3} = field(ll+(3-1),jj)*r0;		% mm = 6 
						temp_field{ii,jj,kk,4} = field(ll+(4-1),jj)*mu0;	% mm = 1 
						temp_field{ii,jj,kk,5} = field(ll+(5-1),jj)*r0*g0;	% mm = 2 
						temp_field{ii,jj,kk,6} = field(ll+(6-1),jj)*g0;		% mm = 3 
					else % surface  
						ll = ii*(nlayer+1)*6 - 2; 
						temp_field{ii,jj,kk,1} = field(ll+(1-1),jj)*r0;	
						temp_field{ii,jj,kk,3} = field(ll+(2-1),jj)*r0;	
						temp_field{ii,jj,kk,5} = field(ll+(3-1),jj)*r0*g0; 
						% surface BC 
						temp_field{ii,jj,kk,4} = 0; 
						if (md.love.forcing_type==9)
							temp_field{ii,jj,kk,2} = 0; 
							temp_field{ii,jj,kk,6} = (2*ii-1)/r0 - ii*field(ll+(3-1),jj)*g0;
						elseif (md.love.forcing_type==11)
							temp_field{ii,jj,kk,2} = -(2*(ii-1)+1)*rho_avg/3; 
							temp_field{ii,jj,kk,6} = (2*ii-1)/r0 - ii*field(ll+(3-1),jj)*g0; 
						end
					end
				end
			end
		end
		field=temp_field; 
	end

	if time~=-9999,
		time=time/yts;
	end

	result.fieldname=fieldname;
	result.time=time;
	result.step=step;
	result.field=field;
end
% }}}
function result=ReadDataDimensions(fid) % {{{
%READDATADIMENSIONS - read data dimensions, step and time, but not the data itself.
%
%   Usage:
%      field=ReadDataDimensions(fid)

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
	N=1; %default
	if type==1,
		fseek(fid,M*8,0);
	elseif type==2,
		fseek(fid,M,0);
	elseif type==3,
		N=fread(fid,1,'int');
		fseek(fid,N*M*8,0);
	else
		error(['cannot read data of type ' num2str(type) ]);
	end

	result.fieldname=fieldname;
	result.time=time;
	result.step=step;
	result.M=M;
	result.N=N;
end
% }}}

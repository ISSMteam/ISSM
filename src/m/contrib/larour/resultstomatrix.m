function matrix=resultstomatrix(md,resultname,field,varargin)
%RESULTSTOMATRIX - go grab in the model results structure the vector results for each time step (which is not empty), 
%                  and line them up in a matrix.  If time vector is provided, resample.
%
%   Usage:
%      matrix=resultstomatrix(model,solutioname,fieldname)
%
%   Available options:
%      - 'time'     : vector providing new time tags used to resample time
%
%   Example:
%      vel=resultstomatrix(md,'TransientSolution','Vel');
%      vel=resultstomatrix(md,'TransientSolution','Vel','time',2008:1/12:2014);
%
%   See also MODEL  resample


	options=pairoptions(varargin{:});

	results=md.results.(resultname);

	%first, figure out the size: 
	count=0;
	nods=0;
	for i=1:length(results),
		if ~isempty(results(i).(field)),
			count=count+1;
			nods=size(results(i).(field),1);
		end
	end

	if ~count, 
		error(['could not find any result ' field ' in ' resultname]);
	end

	%initialize: 
	matrix=zeros(nods+1,count);

	%fill it up: 
	count=0;
	for i=1:length(results),
		if ~isempty(results(i).(field)),
			count=count+1;
			matrix(1:end-1,count)=results(i).(field);
			matrix(end,count)=results(i).time;
		end
	end

	newtime=getfieldvalue(options,'time',[]);
	if ~isempty(newtime),
		newmatrix=zeros(nods+1,length(newtime));
		newmatrix(end,:)=newtime;
		%we are asked to reinterpolate to this new time: 

		for i=1:nods,
			warning off;
			ts=timeseries(matrix(i,:), matrix(end,:));
			ts=resample(ts,newtime);
			warning on;
			newmatrix(i,:)=ts.Data;
		end

		matrix=newmatrix;
	end

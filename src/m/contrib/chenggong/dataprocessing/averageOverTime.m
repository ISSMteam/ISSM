function averagedData = averageOverTime(data, time, startP, endP)
%averageOverTime - compute the time averaged value of data in the range [startP, endP]
%
%   data: time dependent data,
%   time: has same number of columns as data, sorted in ascending order
%   startP: start time point for averaging, can be at any time point, or
%   even not coincide with the given points in 'time'
%   endP: end time point for averaging. if this is the same as the startP,
%   or not given, then the output will be the linear interpolation of the
%   data at startP.
%
% if out of the given time range, use constant extrapolation,
% if within the range, do linear interpolation.
%
if nargin < 4
    endP = startP;
end

Nt = length(time);
if Nt ~= size(data,2)
    error('The data and time does not match!');
end

if startP > endP
    error('The given start point has to be smaller than the end point!');
elseif startP == endP
    % extrapolation if earlier than the first time point
    if startP <= time(1)
        i1 = 1;
        i2 = 2;
    % extrapolation if later than the last time point
    elseif startP >= time(Nt)
        i1 = Nt-1;
        i2 = Nt;
    else
        pos = (time > startP);
        i2 = find(pos, 1, 'first');
        i1 = i2 - 1;
    end
    % linear interpolation 
    averagedData = data(:, i1) + (data(:, i1)-data(:, i2)).*(startP-time(i1))./(time(i1)-time(i2));
else
    % find the first and last ID in the time series within the given range
    pos = ((time>=startP) & (time<=endP));
    firstId = find(pos, 1, 'first');
    lastId = find(pos, 1, 'last');
    if isempty(firstId) | isempty(lastId)
		 averagedData = nan(size(data,1), 1);
       % error('Time points out of range');
	 else
		 % computed the integral with trapzoidal rule
   	 integData = zeros(size(data,1), 1);
   	 for i = firstId:lastId-1
   	     integData = integData + (data(:,i+1) + data(:,i)) .* (time(i+1) - time(i)).*0.5;
   	 end
   	 
   	 % special treatment for the first and last time step, which are not
   	 % complete steps.
   	 if firstId == 1
   	     % extrapolation
   	     integData = integData + data(:,1) .* (time(1)-startP);
   	 else
   	     % interpolation
   	     integData = integData + 0.5.*(time(firstId)-startP)./ ...
   	         (time(firstId)-time(firstId-1)) .* ...
   	         ((time(firstId)-startP) .* data(:,firstId-1) + ...
   	         (time(firstId)+startP-2*time(firstId-1)).*data(:,firstId));
   	 end
   	 
   	 if lastId == Nt
   	     % extrapolation
   	     integData = integData + data(:,Nt) .* (endP-time(Nt));
   	 else
   	     % interpolation
   	     integData = integData + 0.5.*(time(lastId)-endP)./ ...
   	         (time(lastId+1)-time(lastId)) .* ...
   	         ((time(lastId)-endP) .* data(:,lastId+1) + ...
   	         (time(lastId)+endP-2*time(lastId+1)).*data(:,lastId));
   	 end
   	 
   	 averagedData = integData ./ (endP-startP);
	end
end


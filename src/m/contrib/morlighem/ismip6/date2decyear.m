function decyear = date2decyear(datein)
%DATE2DECYEAR - convert date to decimal year (for ISSM)
%
%   input argument must be a serial date (see datenum)
%
%   Usage:
%      year = date2decyear(datein)
%
%   Example:
%      year = date2decyear(datenum('19-May-2000'))


%Make table from date coming in
timevec = datevec(datein);

%Set everything in the date vector to 0 except for the year
%to compute the date for the beginninf of the year
timevec(:,2:end) = 0;
dateYearBegin = datenum(timevec);

%Compute date of end of year
timevec2 = timevec;
timevec2(:,1) = timevec2(:,1) + 1;
dateYearEnd = datenum(timevec2);

%Calculate the day of the year
doy = datein - dateYearBegin;

%Finally, we can create the decimal year time
decyear =  timevec(:,1) + (doy - 1) ./ (dateYearEnd - dateYearBegin);

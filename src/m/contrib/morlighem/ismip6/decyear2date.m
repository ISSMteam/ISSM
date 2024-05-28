function dateout = decyear2date(decyear)
%DECYEAR2DATE - converts decimal year to serial date
%
%   Usage:
%      dateout = decyear2date(decyear)
%
%   Example:
%      dateout = decyear2date(2011.001)

decyear = decyear(:);

%Get year
year = floor(decyear);
fraction = mod(decyear,1);

%Get date of beginning and end of year
date0 = datenum((num2str(year)),'yyyy');
date1 = datenum((num2str(year+1)),'yyyy');

%Compute number of days in year
numdays = date1 - date0;

%Compute date
dateout= date0 + fraction .* numdays;

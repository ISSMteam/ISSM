function returnvalue=isnans(array)
%ISNANS - figure out if an array is nan; wrapper to isnan from MATLAB, which does not allow this test for structures
%
%   Usage:
%      returnvalue=isnans(array);
%
%   See Also:
%      ISNAN

if isstruct(array) 
	returnvalue=0;
elseif iscell(array)
	returnvalue=0;
else
	returnvalue=isnan(array);
end

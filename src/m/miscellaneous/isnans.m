function returnvalue=isnans(array)
%ISNANS: figure out if an array is nan. wrapper to isnan from matlab which stupidly does not allow this test  for structures!
%
%  Usage:    isnans(array)
%
%  See also : ISNAN 

if isstruct(array), 
	returnvalue=0;
elseif iscell(array)
	returnvalue=0;
else
	returnvalue=isnan(array);
end

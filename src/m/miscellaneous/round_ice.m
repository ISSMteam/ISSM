function new_x=round_ice(x,numnonzeros)
%ROUND_ICE - rounds up x so that it has only numnonzeros non zero digits
%
%   numnonzeros must be an integer larger or equal to 1
%
%   Usage:
%      new_x=round_ice(x,numnonzeros)

%some checks
if (nargin ~=2 | nargout>1),
	error('round_ice usage: new_x=round_ice(x,numonzeros)');
end
if ~isnumeric(x)
	error('round_ice error message: x must be a number and numzeros an integer');
end
if round(numnonzeros)~=numnonzeros
	error('round_ice error message: numnonzeros must be an integer larger or equal to 1')
end
if any(numnonzeros<1)
	error('round_ice error message: numnonzeros must be an integer larger or equal to 1')
end
if (length(numnonzeros)~=1 & size(numnonzeros)~=size(x))
	error('round_ice error message: numnonzeros must be an integer larger or equal to 1 or a list of integers of length length(x)')
end

%figure out how long x is
lengthx=ceil(log10(abs(x)));

%if x contains 0, lengthx=-Inf
lengthx(isinf(lengthx))=1;

%get its sign
si=sign(x);

%rule out zeros
new_x=si.*round(abs(x).*10.^(-lengthx+numnonzeros)).*10.^(lengthx-numnonzeros);

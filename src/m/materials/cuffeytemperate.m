function rigidity=cuffeytemperate(temperature, waterfraction, stressexp)
%CUFFEYTEMPERATE - calculates ice rigidity as a function of temperature and waterfraction
%
%   rigidity (in s^(1/3)Pa) is the flow law parameter in the flow law sigma=B*e(1/3)
%   (Cuffey and Paterson, p75). 
%   temperature is in Kelvin degrees
%
%   Usage:
%      rigidity=cuffeytemperate(temperature, waterfraction, stressexp)

if (any(size(temperature)~=size(waterfraction))),
	error('input temperature and waterfraction should have same size!');
end
if any(temperature<0)
	error('input temperature should be in Kelvin (positive)');
end
if any(waterfraction<0 | waterfraction>1)
	error('input waterfraction should be between 0 and 1');
end

rigidity=cuffey(temperature).*(1+181.25*max(0., min(0.01, waterfraction))).^(-1/stressexp);

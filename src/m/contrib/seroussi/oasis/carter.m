function attenuation=carter(temperature)
%CARTER - attenuation as a function of temperature
%
%   TWO WAY - Attenuation (in dB/m) as a function of temperature (K)
%   From Carter at al. 2007 (Radar-based subglacial lake classification in Antarctica)
%   Figure 4
%
%   Usage:
%      attenuation=carter(temperature)

if(temperature<0)
	error('input temperature should be in Kelvin (positive)');
end
T=temperature-273.15;

Temp=[-50 -45 -40 -35 -30 -25 -20 -15 -10 -5 0]';
A=[0.0015 0.002 0.003 0.0042 0.0055 0.0083 0.012 0.0175 0.026 0.038 0.055]';

%Now, do a cubic fit between Temp and B: 
[cfun,gof,output]=fit(Temp,A,'cubicspline');
%breaks=cfun.p.breaks;
%coeff=cfun.p.coefs;

%Calculate attenuation
attenuation=cfun(T);

%Make it a 2 way attenuation
attenuation=2*attenuation;

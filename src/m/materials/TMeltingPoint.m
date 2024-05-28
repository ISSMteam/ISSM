function Tm=TMeltingPoint(reftemp, pressure)
%TMELTINGPOINT- calculate pressure melting point of ice
%
%   reftemp is the melting temperature at atmospheric pressure (initialized in md.materials.meltingpoint)   
%
%   pressure is in Pa   
%
%   Usage:
%   Tm=TMeltingPoint(md.materials.meltingpoint,pressure)

%variables
beta=7.9e-8; % K Pa^-1

Tm=reftemp-beta*pressure;


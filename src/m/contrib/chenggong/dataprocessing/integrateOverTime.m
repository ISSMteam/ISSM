function [newData, intData] = integrateOverTime(tdata, data, time, Nintdt)
% integrateOverDomain - integrating data over the whole domain
%
%	 Input:
%		tdata:		time steps of data
%		data:			time series of the data 
%		time:			time steps to be used for integration
%		Nintdt:		number of time steps in 'time' to be used for integration
%
%	 Return:
%		newData:		data projected on to 'time'
%		intData:		integrated 'newData'	

dt = abs(time(2) -time(1));
% patch the first N time steps
intTime = [[-Nintdt:-1]*dt+time(1), time];
% integration scheme
integ = ones(1, Nintdt) * dt;

% interpolate data from the time scheme 'tdata' to 'time'
newData = interp1(double(tdata), double(data), intTime, 'linear', 'extrap'); 

% integrated, only works for equidistance time steps for now
intData = conv(integ, newData);
intData = intData(Nintdt:end-Nintdt);
newData = newData(Nintdt+1: end);

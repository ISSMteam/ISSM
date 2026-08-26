function [p]=normcdf_issm(x,mu,sigma)
%NORMCDF_ISSM - wrapper for normcdf to avoid using the MATLAB statistics toolbox
%
%   Usage:
%      p=normcdf_issm(x,mu,sigma);

	p=(1.+erf((x-mu)/(sigma*sqrt(2.))))/2.;

end

%
%  wrapper for normcdf to avoid using the matlab statistics toolbox.
%
function [p]=normcdf_issm(x,mu,sigma)

	p=(1.+erf((x-mu)/(sigma*sqrt(2.))))/2.;

end

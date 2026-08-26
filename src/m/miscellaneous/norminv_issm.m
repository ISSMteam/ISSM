function [x]=norminv_issm(p,mu,sigma)
%NORMINV_ISSM - wrapper for norminv to avoid using the MATLAB statistics toolbox
%
%   Usage:
%      x=norminv_issm(p,mu,sigma);

	x=mu+sigma*sqrt(2.)*erfinv(2.*p-1.);

end

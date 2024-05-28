%
%  wrapper for norminv to avoid using the matlab statistics toolbox.
%
function [x]=norminv_issm(p,mu,sigma)

	x=mu+sigma*sqrt(2.)*erfinv(2.*p-1.);

end

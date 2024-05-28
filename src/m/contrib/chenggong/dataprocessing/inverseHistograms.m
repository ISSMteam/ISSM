function [cdfa, q_Sinv] = inverseHistograms(mu, S, Sinv)
%     Given a distribution mu compute its inverse quantile function
%     Parameters
%     ----------
%     mu     : histogram
%     S      : support of the histogram
%     Sinv   : support of the quantile function
%     Returns
%     -------
%     cdfa   : the cumulative distribution function and
%     q_Sinv : the inverse quantile function of the distribution mu

epsilon = 1e-14;
A = (mu>epsilon);
Sa = S(A);

% cummulative sum
cdf = cumtrapz(S, mu);
cdfa = cdf(A);

% set the first value to 0 and last value to 1
cdfa(1) = 0;
Sa(1) = 0;

if cdfa(end) < 1
    cdfa = [cdfa(:);1];
    Sa = [Sa(:);S(end)];
end

[~, ind] = unique(cdfa);
% linear interpolation
q_Sinv = interp1(cdfa(ind), Sa(ind), Sinv);

end

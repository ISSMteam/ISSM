%
%  calculate the same moments and confidence intervals that dakota
%  calculates for a sample.
%
%  [dresp                      ]=dakota_moments(dresp,alpha)
%  [mean,stddev,meanci,stddevci]=dakota_moments(samp ,alpha)
%
%  the required input is:
%    dresp         (structure array, responses)
%      or
%    samp          (double array, lists of samples)
%
%  and the optional input is:
%    alpha         (numeric, confidence interval of 100(1-alpha)%)
%
%  the required field of dresp is:
%    sample        (double vector, list of samples)
%
%  the required output is:
%    dresp         (structure array, responses)
%      or
%    mean          (double, mean of sample)
%    stddev        (double, standard deviation of sample)
%    meanci(2)     (double, confidence interval of mean)
%    stddevci(2)   (double, confidence interval of standard deviation)
%
%  and the output fields of dresp are:
%    mean          (double, mean of sample)
%    stddev        (double, standard deviation of sample)
%    meanci(2)     (double, confidence interval of mean)
%    stddevci(2)   (double, confidence interval of standard deviation)
%
%  for each response (or column of data) in the input array, this
%  function calculates the mean, standard deviation, and their
%  confidence intervals for a normal distribution, the same way as
%  dakota would.  if the input is a structure, the output is fields
%  in the structure; if the input is an array, the output is arrays.
%
%  dresp data would typically be contained in the dakota tabular
%  output file from a sampling analysis, read by dakota_out_parse.
%
%  "Copyright 2009, by the California Institute of Technology.
%  ALL RIGHTS RESERVED. United States Government Sponsorship
%  acknowledged. Any commercial use must be negotiated with
%  the Office of Technology Transfer at the California Institute
%  of Technology.  (NTR 47078)
%
%  This software may be subject to U.S. export control laws.
%  By accepting this  software, the user agrees to comply with
%  all applicable U.S. export laws and regulations. User has the
%  responsibility to obtain export licenses, or other export
%  authority as may be required before exporting such information
%  to foreign countries or providing access to foreign persons."
%
function [varargout]=dakota_moments(varargin)

if ~nargin
    help dakota_moments
    return
end

%%  process input data

iarg=1;
if     isstruct(varargin{iarg})
    dresp=varargin{iarg};
    iarg=iarg+1;
elseif isnumeric(varargin{iarg})
    samp=varargin{iarg};
    iarg=iarg+1;
else
    error(['Unknown data ''' inputname(1) ''' of type ''' class(varargin{1}) '''.']);
end

if iarg <= nargin && isnumeric(varargin{iarg})
    alpha=varargin{iarg};
    iarg=iarg+1;
else
%  use dakota default of 95%
    alpha=0.05;
end

%%  calculate the moments and confidence intervals by input type

if     exist('dresp','var') && ~isempty(dresp)
    for i=1:length(dresp)
        [dresp(i).mean,dresp(i).stddev,...
         dresp(i).meanci,dresp(i).stddevci]=...
            moments_calc(dresp(i).sample,alpha);
    end

    varargout{1}=dresp;

elseif exist('samp','var') && ~isempty(samp)
    mean    =zeros(1,size(samp,2));
    stddev  =zeros(1,size(samp,2));
    meanci  =zeros(2,size(samp,2));
    stddevci=zeros(2,size(samp,2));

%  could do this using vector math rather than loop, but want to allow
%  the case of differently sized samples padded by NaN's

    for i=1:size(samp,2)
        [mean(i),stddev(i),...
         meanci(:,i),stddevci(:,i)]=...
            moments_calc(samp(:,i),alpha);
    end

    varargout{1}=mean;
    varargout{2}=stddev;
    varargout{3}=meanci;
    varargout{4}=stddevci;
else
    error(['Empty data ''' inputname(1) ''' of type ''' class(varargin{1}) '''.']);
end

end

%%  function to calculate the results

function [mu,sigma,muci,sigmaci]=moments_calc(samp,alpha)

%  remove any NaN padding (should only occur at end)

    samp=samp(~isnan(samp(:)));
    nsamp=length(samp);
    prob=1-alpha/2;

%  could use Matlab normfit, but make calculations explicit
%     [mu,sigma,muci,sigmaci]=normfit(samp,alpha);

    mu   =mean(samp);
    sigma=std(samp);

	try
        muci(1,1)   =mu-tinv(prob,nsamp-1)*sigma/sqrt(nsamp);
        muci(2,1)   =mu+tinv(prob,nsamp-1)*sigma/sqrt(nsamp);
        sigmaci(1,1)=sigma*sqrt((nsamp-1)/chi2inv(prob  ,nsamp-1));
        sigmaci(2,1)=sigma*sqrt((nsamp-1)/chi2inv(1-prob,nsamp-1));
	catch me
        muci(1,1)   =mu;
        muci(2,1)   =mu;
        sigmaci(1,1)=sigma;
        sigmaci(2,1)=sigma;
	end

end

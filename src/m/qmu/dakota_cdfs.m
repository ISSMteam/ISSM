%
%  calculate the same cumulative distribution functions that dakota
%  calculates for given responses, probabilities, reliabilities, and/or
%  general reliabilities.
%
%  [dresp]=dakota_cdfs(method,dresp      ,resp,prob,rel,grel)
%  [cdf  ]=dakota_cdfs(method,samp       ,resp,prob,rel,grel)
%  [cdf  ]=dakota_cdfs(method,mean,stddev,resp,prob,rel,grel)
%
%  the required input is:
%    method        (char, 'nond_sampling' or 'nond_local_reliability')
%    dresp         (structure array, responses)
%      or
%    samp          (double array, lists of samples)
%      or
%    mean          (double vector, means)
%    stddev        (double vector, standard deviations)
%    resp          (double vector, list of responses)
%    prob          (double vector, list of probabilities)
%    rel           (double vector, list of reliabilities)
%    grel          (double vector, list of general reliabilities)
%
%  and the optional input is:
%    alpha         (numeric, confidence interval of 100(1-alpha)%)
%
%  the required field of dresp is (for nond_sampling):
%    sample        (double vector, list of samples)
%  or (for nond_local_reliability):
%    mean          (double, mean of samples)
%    stddev        (double, standard deviation of samples)
%
%  the required output is:
%    dresp         (structure array, responses)
%      or
%    cdf(:,4)      (double, array of resp/prob/rel/grel)
%
%  and the output fields of dresp are:
%    cdf(:,4)      (double, array of resp/prob/rel/grel)
%
%  for each response (or column of data) in the input array, this
%  function calculates the responses, probabilities, reliabilities
%  and general reliabilities for the cumulative distribution function,
%  the same way as dakota would.  if the input is a structure, the
%  output is a field in the structure; if the input is arrays, the
%  output is an array.
%
%  dresp data would typically be contained in the dakota tabular
%  output file from a sampling analysis, or in the dakota output file
%  from a local reliability analysis, either read by dakota_out_parse.
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
function [varargout]=dakota_cdfs(varargin)

if ~nargin
    help dakota_cdfs
    return
end

%%  process input data

iarg=1;
if ischar(varargin{iarg})
    method=varargin{iarg};
    iarg=iarg+1;
    if ~strncmpi(method,'nond_s',6) && ~strncmpi(method,'nond_l',6)
        error(['Method ''' method ''' is unrecognized.']);
    end
else
    method='';
end

if isstruct(varargin{iarg})
    dresp=varargin{iarg};
    iarg=iarg+1;
else
    if     strncmpi(method,'nond_s',6)
        samp=varargin{iarg};
        iarg=iarg+1;
    elseif strncmpi(method,'nond_l',6)
        mean  =varargin{iarg};
        iarg=iarg+1;
        stddev=varargin{iarg};
        iarg=iarg+1;
    end
end

if iarg <= nargin && isnumeric(varargin{iarg})
    resp=varargin{iarg};
    iarg=iarg+1;
else
    resp=[];
end
if iarg <= nargin && isnumeric(varargin{iarg})
    prob=varargin{iarg};
    iarg=iarg+1;
else
    prob=[];
end
if iarg <= nargin && isnumeric(varargin{iarg})
    rel =varargin{iarg};
    iarg=iarg+1;
else
    rel =[];
end
if iarg <= nargin && isnumeric(varargin{iarg})
    grel=varargin{iarg};
    iarg=iarg+1;
else
    grel=[];
end

%%  calculate the cumulative distribution functions by input type

if     exist('dresp','var') && ~isempty(dresp)
    if     strncmpi(method,'nond_s',6)
        for i=1:length(dresp)
            [dresp(i).cdf]=cdfs_samp_calc(dresp(i).sample,...
                resp,prob,rel,grel);
        end
    elseif strncmpi(method,'nond_l',6)
        for i=1:length(dresp)
            [dresp(i).cdf]=cdfs_lr_calc(dresp(i).mean,dresp(i).stddev,...
                resp,prob,rel,grel);
        end
    end

    varargout{1}=dresp;

elseif exist('samp','var') && ~isempty(samp)
    cdf=zeros(length(resp)+length(prob)+length(rel)+length(grel),...
              4,size(samp,2));

    for i=1:size(samp,2)
        [cdf(:,:,i)]=cdfs_samp_calc(samp(:,i),...
            resp,prob,rel,grel);
    end

    varargout{1}=cdf;

elseif exist('mean','var'  ) && ~isempty(mean  ) && ...
       exist('stddev','var') && ~isempty(stddev)
    cdf=zeros(length(resp)+length(prob)+length(rel)+length(grel),...
              4,length(mean));

    for i=1:length(mean)
        [cdf(:,:,i)]=cdfs_lr_calc(mean(i),stddev(i),...
            resp,prob,rel,grel);
    end

    varargout{1}=cdf;
else
    error(['Empty data ''' inputname(2) ''' of type ''' class(varargin{2}) '''.']);
end

end

%%  function to calculate the results for a sampling analysis

function [cdf]=cdfs_samp_calc(samp,resp,prob,rel,grel)

%  sort the samples and remove any NaN padding (should only occur at end)

    samp=sort(samp(~isnan(samp(:))),'ascend');
    nsamp=length(samp);

    mu   =mean(samp);
    sigma=std(samp);

    cdf=zeros(length(resp)+length(prob)+length(rel)+length(grel),4);
    cdf(:,:)=NaN;
    irow=0;

%  compute quantities, given response levels

    for i=1:length(resp)
        irow=irow+1;
        indx=bin_search_val(resp(i),samp);
        cdf(irow,1)=resp(i);
        cdf(irow,2)=indx/nsamp;
        cdf(irow,3)=(mu-resp(i))/sigma;
%        cdf(irow,4)=-sqrt(2)*erfinv((indx-nsamp/2)/(nsamp/2));
        cdf(irow,4)=sqrt(2)*erfcinv(indx/(nsamp/2));
    end

%  compute response levels, given probabilities

    for i=1:length(prob)
        irow=irow+1;
%  why not round(prob(i)*(nsamp-1)+1)?
        indx=ceil(prob(i)*(nsamp));
        if     (indx < 1)
            indx=1;
        elseif (indx > nsamp)
            indx=nsamp;
        end
        cdf(irow,1)=samp(indx);
        cdf(irow,2)=prob(i);
    end

%  compute response levels, given reliabilities

    for i=1:length(rel)
        irow=irow+1;
        cdf(irow,1)=mu-sigma*rel(i);
        cdf(irow,3)=rel(i);
    end

%  compute response levels, given general reliabilities

    for i=1:length(grel)
        irow=irow+1;
%         indx=ceil(nsamp/2+nsamp/2*erf(-grel(i)/sqrt(2)));
        indx=ceil((nsamp/2)*erfc(grel(i)/sqrt(2)));
        if     (indx < 1)
            indx=1;
        elseif (indx > nsamp)
            indx=nsamp;
        end
        cdf(irow,1)=samp(indx);
        cdf(irow,4)=grel(i);
    end

end

%%  function to calculate the results for a local reliability analysis

function [cdf]=cdfs_lr_calc(mu,sigma,resp,prob,rel,grel)

    cdf=zeros(length(resp)+length(prob)+length(rel)+length(grel),4);
    irow=0;

%  compute quantities, given response levels

    for i=1:length(resp)
        irow=irow+1;
        cdf(irow,1)=resp(i);
        cdf(irow,2)=normcdf_issm(resp(i),mu,sigma);
        cdf(irow,3)=(mu-resp(i))/sigma;
        cdf(irow,4)=(mu-resp(i))/sigma;
    end

%  compute quantities, given probabilities

    for i=1:length(prob)
        irow=irow+1;
        cdf(irow,1)=norminv_issm(prob(i),mu,sigma);
        cdf(irow,2)=prob(i);
        cdf(irow,3)=-norminv_issm(prob(i),0,1);
        cdf(irow,4)=-norminv_issm(prob(i),0,1);
    end

%  compute quantities, given reliabilities

    for i=1:length(rel)
        irow=irow+1;
        cdf(irow,1)=mu-sigma*rel(i);
        cdf(irow,2)=normcdf_issm(-rel(i),0,1);
        cdf(irow,3)=rel(i);
        cdf(irow,4)=rel(i);
    end

%  compute quantities, given general reliabilities

    for i=1:length(grel)
        irow=irow+1;
        cdf(irow,1)=mu-sigma*grel(i);
        cdf(irow,2)=normcdf_issm(-grel(i),0,1);
        cdf(irow,3)=grel(i);
        cdf(irow,4)=grel(i);
    end

end
%%
%  function to perform a recursive binary search for a matrix of values
%  in an ordered vector (loop separately outside of recursion for
%  efficiency purposes)
%
%  function [ind]=bin_search(val,vect)
%
function [ind]=bin_search(val,vect)

ind=zeros(size(val));

for i=1:numel(val)
    ind(i)=bin_search_val(val(i),vect);
end

end
%%
%  function to perform a recursive binary search in an ordered vector,
%  returning low if value does not exist (more efficient than find or
%  ismember, which must use linear searches and/or sort)
%
%  function [ind]=bin_search_val(val,vect)
%
function [ind]=bin_search_val(val,vect)

imid=floor((1+length(vect))/2);

if (val == vect(imid))
    ind=imid;
elseif (val < vect(imid))
    if (imid > 1)
        ind=     bin_search(val,vect(1:imid-1));
    else
        ind=0;
    end
elseif (val > vect(imid))
    if (imid < length(vect))
        ind=imid+bin_search(val,vect(imid+1:length(vect)));
    else
        ind=length(vect);
    end
else
    ind=NaN;
end

end

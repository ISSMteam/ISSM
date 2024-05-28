%
%  plot a stacked bar chart of the response levels in the cdf
%  for the sample and confidence intervals.
%
%  []=plot_rlev_bars_ci(dresp      ,params)
%  []=plot_rlev_bars_ci(dresp,descr,params)
%  []=plot_rlev_bars_ci(sampr,descr,params)
%
%  where the required input is:
%    dresp         (structure array, responses)
%      or
%    dresp         (structure array, responses)
%    descr         (cell array, list of response descriptions desired)
%      or
%    sampr         (double array, lists of response samples)
%    descr         (cell array, list of response descriptions)
%
%  the required fields of dresp are:
%    descriptor    (char, description)
%    cdf(:,4)      (double matrix, CDF table)
%
%  and the optional fields of dresp are:
%    mean          (double, mean of sample)
%    stddev        (double, standard deviation of sample)
%    meanci(2)     (double, confidence interval of mean)
%    stddevci(2)   (double, confidence interval of standard deviation)
%
%  the optional input is:
%    params        (string/numeric, parameter names and values)
%
%  where the optional parameters are:
%    ymin          (numeric, minimum of y-axis)
%    ymax          (numeric, maximum of y-axis)
%    xtlrot        (numeric, rotation in degrees of x-tick labels)
%    lstr          (cell array, legend labels)
%
%  for each response in the input array, this function plots
%  a stacked bar plot of the responses, where the bars are
%  stacked by the response levels corresponding to the given
%  probabilities in the CDF, and annotates it with the
%  description.  the response levels for the normal distribution
%  and the confidence intervals are also plotted.  the legend
%  labels can be given or constructed from the probabilities.
%
%  dresp data would typically be contained in the dakota tabular
%  output file from a sampling analysis, read by dakota_out_parse.
%
%  "Copyright 2009, by the California Institute of Technology.
%  ALL RIGHTS RESERVED. United States Government Sponsorship
%  acknowledged. Any commercial use must be negotiated with
%  the Office of Technology Transfer at the California Institute
%  of Technology.  (J. Schiermeier, NTR 47078)
%
%  This software may be subject to U.S. export control laws.
%  By accepting this  software, the user agrees to comply with
%  all applicable U.S. export laws and regulations. User has the
%  responsibility to obtain export licenses, or other export
%  authority as may be required before exporting such information
%  to foreign countries or providing access to foreign persons."
%
function []=plot_rlev_bars_ci(varargin)

if ~nargin
    help plot_rlev_bars_ci
    return
end

%%  process input data and assemble into dresp as needed

%  responses

iarg=1;
if isstruct(varargin{iarg})
    dresp=varargin{iarg};
    iarg=iarg+1;

%     if iarg <= nargin && (iscell(varargin{iarg}) || ischar(varargin{iarg}))
    if iarg <= nargin && iscell(varargin{iarg})
        dresp=struc_desc(dresp,varargin{iarg});
        iarg=iarg+1;
    end
else
    sampr=varargin{iarg};
    iarg=iarg+1;

    if     iarg <= nargin && iscell(varargin{iarg})
        descr=varargin{iarg};
        iarg=iarg+1;
%     elseif iarg <= nargin && ischar(varargin{iarg})
%         descr=cellstr(varargin{iarg});
%         iarg=iarg+1;
    else
        descr=cell(1:size(sampr,2));
    end

    dresp=struct([]);
    for i=1:size(sampr,2)
        dresp(end+1).sample=sampr(:,i);
        if ~isempty(descr)
            dresp(i).descriptor=descr{i};
        else
            dresp(i).descriptor=['dresp_' num2str(i)];
        end
    end
end

%  parameters

while (iarg <= nargin-1)
    if ischar(varargin{iarg})
        if ~isempty(strmatch(varargin{iarg},...
                {'ymin','ymax','xtlrot','lstr'},...
                'exact'))
            eval([varargin{iarg} '=varargin{iarg+1};']);
            disp([varargin{iarg} '=' any2str(varargin{iarg+1}) ';']);
        else
            warning([varargin{iarg} '=' any2str(varargin{iarg+1}) ' is not recognized.']);
        end
    else
        error(['''' any2str(varargin{iarg}) ''' is not a parameter name.']);
    end
    iarg=iarg+2;
end

%%  calculate any missing information (noting that dresp is local)

for i=1:length(dresp)
    if ~isfield(dresp(i),'mean') || isempty(dresp(i).mean) || ...
       ~isfield(dresp(i),'stddev') || isempty(dresp(i).stddev) || ...
       ~isfield(dresp(i),'meanci') || isempty(dresp(i).meanci) || ...
       ~isfield(dresp(i),'stddevci') || isempty(dresp(i).stddevci)
%  calculate 95% confidence intervals (same as dakota)
        [dresp(i).mean,dresp(i).stddev,...
         dresp(i).meanci,dresp(i).stddevci]=...
            normfit_issm(sampr(:,i),0.05);
        display('Using calculated normal fits from sample data.')
    end

    if ~isfield(dresp(i),'cdf') || isempty(dresp(i).cdf)
%  use minus/plus integer standard deviations
        sdvect=[-4 -3 -2 -1 0 1 2 3 4];
        dresp(i).cdf(:,2)=normcdf_issm(sdvect,0,1);
        dresp(i).cdf(:,1)=norminv_issm(dresp(i).cdf(:,2),...
                                       dresp(i).mean,dresp(i).stddev);
        display('Using integer standard deviations for percentages.')

        if ~exist('lstr','var') || isempty(lstr)
            lstr=cell(1,size(dresp(i).cdf,1));
            for j=1:size(dresp(i).cdf,1)
                if sdvect(j)
                    lstr{j}=sprintf('mu %+d sigma',sdvect(j));
                else
                    lstr{j}='mu';
                end
            end
        end
    end
end

%%  assemble the data into a matrix and calculate the increments

descr=cell (1,0);
lcdfr=zeros(1,length(dresp));
for i=1:length(dresp)
    lcdfr(i)=size(dresp(i).cdf,1);
end
cdfr=zeros(0,max(lcdfr));

%  fill in the cdf data

for i=1:length(dresp)
    if ~isempty(dresp(i).cdf)
        descr(end+1)=cellstr([dresp(i).descriptor]);
        cdfr(end+1,:)=dresp(i).cdf(:,1);
        if isfield(dresp(i),'mean'  ) && ~isempty(dresp(i).mean  ) && ...
           isfield(dresp(i),'stddev') && ~isempty(dresp(i).stddev)
            descr(end+1)=cellstr([dresp(i).descriptor ' norm']);
            cdfr(end+1,:)=norminv_issm(dresp(i).cdf(:,2),dresp(i).mean,dresp(i).stddev);
        end
        if isfield(dresp(i),'meanci'  ) && ~isempty(dresp(i).meanci  ) && ...
           isfield(dresp(i),'stddevci') && ~isempty(dresp(i).stddevci)
            descr(end+1)=cellstr([dresp(i).descriptor ' norm-+']);
            descr(end+1)=cellstr([dresp(i).descriptor ' norm--']);
            descr(end+1)=cellstr([dresp(i).descriptor ' norm+-']);
            descr(end+1)=cellstr([dresp(i).descriptor ' norm++']);
            cdfr(end+1,:)=norminv_issm(dresp(i).cdf(:,2),dresp(i).meanci(1),dresp(i).stddevci(2));
            cdfr(end+1,:)=norminv_issm(dresp(i).cdf(:,2),dresp(i).meanci(1),dresp(i).stddevci(1));
            cdfr(end+1,:)=norminv_issm(dresp(i).cdf(:,2),dresp(i).meanci(2),dresp(i).stddevci(1));
            cdfr(end+1,:)=norminv_issm(dresp(i).cdf(:,2),dresp(i).meanci(2),dresp(i).stddevci(2));
        end
    end
end

%  calculate the increments

for i=1:size(cdfr,1)
    for j=find(cdfr(i,:),1,'last'):-1:2
        cdfr(i,j)=cdfr(i,j)-cdfr(i,j-1);
    end
end

%%  draw the stacked bar plot

%  if there's only one row, Matlab 7.5 interprets it as a column,
%  so add an extra row, then reduce xlim

if length(descr) == 1
    cdfr=[cdfr; cdfr];
end

figure
hl1=bar(cdfr,'stacked');
%  set barseries properties for lowest value
whitebg('white')
set(hl1(1),'FaceColor','white')
set(hl1(1),'Visible','off')

ax1=gca;
if length(descr) == 1
    set(ax1,'xlim',[0.5 1.5])
end

ylim('auto')
[ylims]=ylim;
if exist('ymin','var')
    ylims(1)=ymin;
end
if exist('ymax','var')
    ylims(2)=ymax;
end
ylim(ylims)

set(ax1,'xtick',1:length(descr))
set(ax1,'xticklabel',descr)
if exist('xtlrot','var')
    htl=rotateticklabel(ax1,xtlrot);
    tlext=zeros(length(htl),4);
    for i=1:length(htl)
        tlext(i,:)=get(htl(i),'Extent');
    end
end

%  add the annotation

title('Response Levels for Specified Probabilities (PMA)');
xlabel('Response');
if exist('xtlrot','var')
    xlext=get(get(ax1,'xlabel'),'Extent');
    nskip=ceil(max(tlext(:,4))/xlext(4));
    xlabel(cellstr([repmat('        ',nskip,1);'Response']));
    clear nskip xlext tlext
end
ylabel('Response Level');

if ~exist('lstr','var') || isempty(lstr)
    lstr=cell(1,max(lcdfr));
    for i=1:max(lcdfr)
        lstr(i)=cellstr(sprintf('%g%%',...
            100*dresp(find(lcdfr == max(lcdfr),1,'first')).cdf(i,2)));
    end
    if ~isempty(find(lcdfr < max(lcdfr),1,'first'))
        warning('Variable number of probabilities for responses.');
    end
end

hleg1=legend(ax1,lstr,'Location','EastOutside',...
             'Orientation','vertical','Interpreter','none');

end

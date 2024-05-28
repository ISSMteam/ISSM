%
%  plot a stacked bar chart of the response levels in the cdf.
%
%  []=plot_rlev_bars(dresp      ,params)
%  []=plot_rlev_bars(dresp,descr,params)
%
%  where the required input is:
%    dresp         (structure array, responses)
%      or
%    dresp         (structure array, responses)
%    descr         (cell array, list of response descriptions desired)
%
%  the required fields of dresp are:
%    descriptor    (char, description)
%    cdf(:,4)      (double matrix, CDF table)
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
%  description.  the legend labels can be given or constructed
%  from the probabilities.
%
%  this data would typically be contained in the dakota output
%  file and read by dakota_out_parse.
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
function []=plot_rlev_bars(varargin)

if ~nargin
    help plot_rlev_bars
    return
end

%%  process input data and assemble into matrices and increments

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

    descr=cell (1,length(dresp));
    lcdfr=zeros(1,length(dresp));
    for i=1:length(dresp)
        lcdfr(i)=size(dresp(i).cdf,1);
    end
    cdfr=zeros(length(dresp),max(lcdfr));

    for i=1:length(dresp)
        descr(i)=cellstr(dresp(i).descriptor);
        if ~isempty(dresp(i).cdf)
            cdfr(i,1)=dresp(i).cdf(1,1);
            for j=2:size(dresp(i).cdf,1)
                if (dresp(i).cdf(j,1) > dresp(i).cdf(j-1,1))
                    cdfr(i,j)=dresp(i).cdf(j,1)-dresp(i).cdf(j-1,1);
                end
            end
        end
    end
else
    error(['''' inputname(iarg) ''' is not a structure.']);
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

%%  draw the stacked bar plot

%  if there's only one row, Matlab 7.5 interprets it as a column,
%  so add an extra row, then reduce xlim

if length(dresp) == 1
    cdfr=[cdfr; cdfr];
end

figure
hl1=bar(cdfr,'stacked');
%  set barseries properties for lowest value
whitebg('white')
set(hl1(1),'FaceColor','white')
set(hl1(1),'Visible','off')

ax1=gca;
if length(dresp) == 1
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

title('Response Levels for Specified Probabilities (PMA)')
xlabel('Response');
if exist('xtlrot','var')
    xlext=get(get(ax1,'xlabel'),'Extent');
    nskip=ceil(max(tlext(:,4))/xlext(4));
    xlabel(cellstr([repmat('        ',nskip,1);'Response']));
    clear nskip xlext tlext
end
ylabel('Response Level')

if ~exist('lstr','var') || isempty(lstr)
    lstr=cell(1,max(lcdfr));
    for i=1:max(lcdfr)
        lstr(i)=cellstr(sprintf('%g%%',...
            100*dresp(find(lcdfr == max(lcdfr),1,'first')).cdf(i,2)));
    end
    if ~isempty(find(lcdfr < max(lcdfr),1))
        warning('Variable number of probabilities for responses.');
    end
end

hleg1=legend(ax1,lstr,'Location','EastOutside',...
             'Orientation','vertical','Interpreter','none');

end

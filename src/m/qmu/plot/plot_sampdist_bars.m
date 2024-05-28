%
%  plot a stacked bar chart of the sample distributions.
%
%  []=plot_sampdist_bars(dresp      ,params)
%  []=plot_sampdist_bars(dresp,descr,params)
%  []=plot_sampdist_bars(sampr,descr,params)
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
%    sample        (double vector, list of samples)
%
%  and the optional fields of dresp are:
%    min           (double, minimum of sample)
%    quart1        (double, first quartile of sample)
%    median        (double, median of sample)
%    quart3        (double, third quartile of sample)
%    max           (double, maximum of sample)
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
%  a stacked bar plot of the list of samples, where the bars
%  are stacked by the four quartiles, and annotates it with
%  the description.  the quartiles will be calculated from the
%  samples if they do not already exist.
%
%  this data would typically be contained in the dakota tabular
%  output file and read by dakota_out_parse.
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
function []=plot_sampdist_bars(varargin)

if ~nargin
    help plot_sampdist_bars
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

if ~isfield(dresp,'min')    || ~isfield(dresp,'quart1') || ...
   ~isfield(dresp,'median') || ~isfield(dresp,'quart3') || ...
   ~isfield(dresp,'max')
    for i=1:length(dresp)
        dresp(i).min   =min         (dresp(i).sample);
        dresp(i).quart1=prctile_issm(dresp(i).sample,25);
        dresp(i).median=median      (dresp(i).sample);
        dresp(i).quart3=prctile_issm(dresp(i).sample,75);
        dresp(i).max   =max         (dresp(i).sample);
    end
end

%%  assemble the data into a matrix and calculate the increments

descr=cell (1,length(dresp));
data =zeros(length(dresp),5);

for i=1:length(dresp)
    descr(i)=cellstr(dresp(i).descriptor);
    data(i,1)=dresp(i).min;
    data(i,2)=dresp(i).quart1-dresp(i).min;
    data(i,3)=dresp(i).median-dresp(i).quart1;
    data(i,4)=dresp(i).quart3-dresp(i).median;
    data(i,5)=dresp(i).max   -dresp(i).quart3;
end

%%  draw the stacked bar plot

%  if there's only one row, Matlab 7.5 interprets it as a column,
%  so add an extra row, then reduce xlim

if length(dresp) == 1
    data=[data; data];
end

figure
hl1=bar(data,'stacked');
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

title('Sample Distributions of Responses')
xlabel('Response')
if exist('xtlrot','var')
    xlext=get(get(ax1,'xlabel'),'Extent');
    nskip=ceil(max(tlext(:,4))/xlext(4));
    xlabel(cellstr([repmat('        ',nskip,1);'Response']));
    clear nskip xlext tlext
end
ylabel('Value')

if ~exist('lstr','var') || isempty(lstr)
    lstr={'minimum' 'quartile 1' 'median' 'quartile 3' 'maximum'};
end

hleg1=legend(ax1,lstr,'Location','EastOutside',...
             'Orientation','vertical','Interpreter','none');

end

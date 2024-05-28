%
%  plot a stacked bar chart of the importance factors.
%
%  []=plot_if_bars(dresp      ,params)
%  []=plot_if_bars(dresp,descr,params)
%
%  where the required input is:
%    dresp         (structure array, responses)
%      or
%    dresp         (structure array, responses)
%    descr         (cell array, list of response descriptions desired)
%
%  the required fields of dresp are:
%    descriptor    (char, description)
%    var           (cell array, variables)
%    impfac        (double array, importance factors)
%
%  the optional input is:
%    params        (string/numeric, parameter names and values)
%
%  where the optional parameters are:
%    ymin          (numeric, minimum of y-axis)
%    ymax          (numeric, maximum of y-axis)
%    ifmin         (double, minimum importance factor)
%    isort         (numeric, sort flag:  0, no sorting;
%                                        1, highest at bottom;
%                                       -1, lowest at bottom)
%    xtlrot        (numeric, rotation in degrees of x-tick labels)
%
%  for each response in the input array, this function plots
%  a stacked bar plot of the responses, where the bars are
%  stacked by the importance factors, and annotates it with the
%  description.  the legend labels are constructed from the
%  variable list.
%
%  this data would typically be contained in the dakota output
%  file from a local reliability analysis and read by
%  dakota_out_parse.
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
function []=plot_if_bars(varargin)

if ~nargin
    help plot_if_bars
    return
end

%%  process input data and assemble into matrices

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
    lifr =zeros(1,length(dresp));
    for i=1:length(dresp)
        lifr(i)=length(dresp(i).impfac);
    end
    ifr =zeros(length(dresp),max(lifr));
    dvar=dresp(find(lifr == max(lifr),1,'first')).var;

    for i=1:length(dresp)
        descr(i)=cellstr(dresp(i).descriptor);
        ifr(i,1:lifr(i))=dresp(i).impfac;
    end
else
    error(['''' inputname(iarg) ''' is not a structure.']);
end

%  parameters

while (iarg <= nargin-1)
    if ischar(varargin{iarg})
        if ~isempty(strmatch(varargin{iarg},...
                {'ymin','ymax','ifmin','isort','xtlrot'},...
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

if ~exist('ifmin','var') || isempty(ifmin)
    ifmin=0;
end

if ~exist('isort','var') || isempty(isort)
    isort=0;
end

%%  sort the data, if necessary

if (isort)
    ifmean=mean(ifr,1);
    if (isort > 0)
        [ifmean,index]=sort(ifmean,'descend');
    else
        [ifmean,index]=sort(ifmean,'ascend' );
    end
    clear ifmean

    dvar=dvar(index);
    ifr =ifr (:,index);
end

%%  filter the data, if necessary

if (ifmin > 0)
    nif=length(dvar);
    dvar(nif+1,1)=cellstr(sprintf('others < %f',ifmin));
    ifr (:,nif+1)=0.;

    nif2=0;
    dvar2=cell (size(dvar));
    ifr2 =zeros(size(ifr ));

%  sum filtered rows and copy unfiltered rows

    for i=1:nif
        if (max(ifr(:,i)) < ifmin)
            ifr(:,nif+1)=ifr(:,nif+1)+ifr(:,i);
        else
            nif2=nif2+1;
            dvar2(nif2)  =dvar(i);
            ifr2 (:,nif2)=ifr (:,i);
        end
    end

%  copy sums

    dvar2(nif2+1)  =dvar(nif+1);
    ifr2 (:,nif2+1)=ifr (:,nif+1);

%  copy back and truncate filtered rows

    clear dvar ifr
    if (isort >= 0)
        dvar(1:nif2+1)  =dvar2(1:nif2+1);
        ifr (:,1:nif2+1)=ifr2 (:,1:nif2+1);
    else
        dvar(1       )  =dvar2(  nif2+1);
        dvar(2:nif2+1)  =dvar2(1:nif2  );
        ifr (:,1       )=ifr2 (:,  nif2+1);
        ifr (:,2:nif2+1)=ifr2 (:,1:nif2  );
    end
    clear nif nif2 dvar2 ifr2
end

%%  draw the stacked bar plot

%  if there's only one row, Matlab 7.5 interprets it as a column,
%  so add an extra row, then reduce xlim

if length(dresp) == 1
    ifr=[ifr; ifr];
end

figure
hl1=bar(ifr,'stacked');

ax1=gca;
if length(dresp) == 1
    set(ax1,'xlim',[0.5 1.5])
end

% set(ax1,'ylim',[0 1.2])
% ylim('auto')
% [ylims]=ylim;
[ylims]=[0 1.2];
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

title('Importance Factors')
xlabel('Response')
if exist('xtlrot','var')
    xlext=get(get(ax1,'xlabel'),'Extent');
    nskip=ceil(max(tlext(:,4))/xlext(4));
    xlabel(cellstr([repmat('        ',nskip,1);'Response']));
    clear nskip xlext tlext
end
ylabel('Importance Factor Value')

hleg1=legend(ax1,dvar,'Location','EastOutside',...
             'Orientation','vertical','Interpreter','none');

end

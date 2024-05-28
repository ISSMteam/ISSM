%
%  plot a line plot or bar chart of the spectrum of importance factors.
%
%  []=plot_if_spec(dresp      ,params)
%  []=plot_if_spec(dresp,descr,params)
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
%    cplot         (char, 'l'/'b'/'g' to plot line/bar)
%    cline         (char, 'off'/'no'/'-'/'--'/':'/etc. to change lines)
%    lwidth        (numeric, line width in points, default 0.5)
%    cmark         (char, 'on'/'yes'/'+'/'o'/'*'/etc. to change markers)
%    msize         (numeric, marker size in points, default 6)
%    ymin          (numeric, minimum of y-axis)
%    ymax          (numeric, maximum of y-axis)
%    ygrid         (char, 'on' to turn on y-grid lines)
%    ylog          (char, 'yes' to use log y-axis)
%    ifmin         (double, minimum importance factor)
%    isort         (numeric, sort flag:  0, no sorting;
%                                        1, highest at right;
%                                       -1, lowest at right)
%    xtlrot        (numeric, rotation in degrees of x-tick labels)
%
%  for each response in the input array, this function plots
%  a bar plot of the importance factors and annotates it with the
%  description.  the legend labels are constructed from the
%  response list.
%
%  this data would typically be contained in the dakota output
%  file from a local reliability analysis and read by
%  dakota_out_parse.
%
%  "Copyright 2010, by the California Institute of Technology.
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
function []=plot_if_spec(varargin)

if ~nargin
    help plot_if_spec
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
                {'cplot','cline','lwidth','cmark','msize',...
                 'ymin','ymax','ygrid','ylog',...
                 'ifmin','isort','xtlrot'},...
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

if ~exist('cplot','var') || isempty(cplot)
    cplot='line';
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
        [ifmean,index]=sort(ifmean,'ascend' );
    else
        [ifmean,index]=sort(ifmean,'descend');
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

%%  draw the line or bar plot

%  if there's only one row, Matlab 7.5 interprets it as a column,
%  so add an extra row, then reduce xlim

if length(dvar) == 1
    ifr=[ifr; ifr];
end

%newplot();
figure
ax1=axes;
if strncmpi(cplot,'l',1)
    if ~exist('cline','var') || strncmpi(cline,'on' ,2) || strncmpi(cline,'y',1)
        cline='-';
    elseif strncmpi(cline,'off',3) || strncmpi(cline,'n',1)
        cline='none';
    end
    if ~exist('lwidth','var')
        lwidth=0.5;
    end
    if ~exist('cmark','var') || strncmpi(cmark,'off',3) || strncmpi(cmark,'n',1)
        cmark='none';
    elseif strncmpi(cmark,'on' ,2) || strncmpi(cmark,'y',1) || ...
           (length(cmark) > 1)
        cmark='+';
    end
    if ~exist('msize','var')
        msize=6;
    end

    hl1=plot(ax1,ifr','LineStyle',cline,'LineWidth',lwidth,...
             'Marker',cmark,'MarkerSize',msize);
else
    hl1=bar(ax1,ifr');
end

ax1=gca;
if length(dvar) == 1
    set(ax1,'xlim',[0.5 1.5])
else
    set(ax1,'xlim',[0.5 length(dvar)+0.5])
end

if exist('ylog','var') && strncmpi(ylog,'y',1)
    set(ax1,'yscale','log')
    ylim('auto')
    [ylims]=ylim;
else
    [ylims]=[0 1.2];
end
if exist('ymin','var')
    ylims(1)=ymin;
end
if exist('ymax','var')
    ylims(2)=ymax;
end
ylim(ylims)
if strncmpi(cplot,'b',1) && exist('ylog','var') && strncmpi(ylog,'y',1)
    set(hl1,'basevalue',ylims(1))
end

set(ax1,'xtick',1:length(dvar))
set(ax1,'xticklabel',dvar)
if exist('xtlrot','var')
    htl=rotateticklabel(ax1,xtlrot);
    tlext=zeros(length(htl),4);
    for i=1:length(htl)
        tlext(i,:)=get(htl(i),'Extent');
    end
end

%  add the annotation

if length(descr) > 1
    title('Importance Factors of Responses')
    hleg1=legend(ax1,descr,'Location','EastOutside',...
                 'Orientation','vertical','Interpreter','none');
else
    title(['Importance Factors of ' descr{1}],'Interpreter','none')
end
xlabel('Variable')
if exist('xtlrot','var')
    xlext=get(get(ax1,'xlabel'),'Extent');
    nskip=ceil(max(tlext(:,4))/xlext(4));
    xlabel(cellstr([repmat('        ',nskip,1);'Variable']));
    clear nskip xlext tlext
end
ylabel('Importance Factor Value')

end

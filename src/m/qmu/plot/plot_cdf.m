%
%  plot the CDF of the responses.
%
%  []=plot_cdf(dresp      ,params)
%  []=plot_cdf(dresp,descr,params)
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
%    cplot         (char, 'p'/'r'/'g' to plot prob/reli/genrel)
%    xmin          (numeric, minimum of x-axis)
%    xmax          (numeric, maximum of x-axis)
%    xgrid         (char, 'on' to turn on x-grid lines)
%    ymin1         (numeric, minimum of y-axis)
%    ymax1         (numeric, maximum of y-axis)
%    ygrid1        (char, 'on' to turn on y-grid lines)
%    ynorm         (char, 'yes' to use normal probability y-axis)
%    yprob         (double vector, list of probabilities for y-axis)
%    cline1        (char, 'off'/'no'/'-'/'--'/':'/etc. to change lines)
%    lwidth1       (numeric, line width in points, default 0.5)
%    cmark1        (char, 'on'/'yes'/'+'/'o'/'*'/etc. to change markers)
%    msize1        (numeric, marker size in points, default 6)
%    ymin2         (numeric, minimum of y-axis)
%    ymax2         (numeric, maximum of y-axis)
%    ygrid2        (char, 'on' to turn on y-grid lines)
%    cline2        (char, 'off'/'no'/'-'/'--'/':'/etc. to change lines)
%    lwidth2       (numeric, line width in points, default 0.5)
%    cmark2        (char, 'on'/'yes'/'+'/'o'/'*'/etc. to change markers)
%    msize2        (numeric, marker size in points, default 6)
%    pdfplt        (char, 'bar'/'line'/'off' for pdf plots)
%    pdfleg        (char, 'off' to turn off pdf legends)
%    cmap          (char or numeric, colormap definition)
%
%  for each response in the input array, this function plots
%  a line plot of the CDF and annotates it with the description.
%
%  this data would typically be contained in the dakota output
%  file and read by dakota_out_parse.
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
function []=plot_cdf(varargin)

if ~nargin
    help plot_cdf
    return
end

%%  process input data and assemble into matrices as needed

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
    error(['''' inputname(iarg) ''' is not a structure.']);
end

%  parameters

while (iarg <= nargin-1)
    if ischar(varargin{iarg})
        if ~isempty(strmatch(varargin{iarg},...
                {'cplot','xmin','xmax','xgrid',...
                 'ymin1','ymax1','ygrid1','ynorm','yprob',...
                 'cline1','lwidth1','cmark1','msize1',...
                 'ymin2','ymax2','ygrid2',...
                 'cline2','lwidth2','cmark2','msize2',...
                 'pdfplt','pdfleg','cmap'},...
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

if ~exist('cplot','var') || strncmpi(cplot,'p',1)
    cplot='Probabilit';
    iplot=2;
elseif strncmpi(cplot,'r',1)
    cplot='Reliabilit';
    iplot=3;
elseif strncmpi(cplot,'g',1)
    cplot='General Reliabilit';
    iplot=4;
end

%  assemble data into matrices, based on parameters

descr=cell (1,length(dresp));
lcdfr=zeros(1,length(dresp));
for i=1:length(dresp)
    lcdfr(i)=size(dresp(i).cdf,1);
end
xdata=zeros(max(lcdfr),length(dresp));
xdata(:,:)=NaN;
ydata=zeros(max(lcdfr),length(dresp));
ydata(:,:)=NaN;

for i=1:length(dresp)
    descr(i)=cellstr(dresp(i).descriptor);
    if ~isempty(dresp(i).cdf)
        xdata(1:lcdfr(i),i)=dresp(i).cdf(:,1);
        if strncmpi(cplot,'p',1) && ...
           exist('ynorm','var') && strncmpi(ynorm,'y',1)
             ydata(1:lcdfr(i),i)=norminv_issm(dresp(i).cdf(:,iplot),0,1);
        else
             ydata(1:lcdfr(i),i)=             dresp(i).cdf(:,iplot);
        end
    end
end

%%  draw the line plot

%newplot();
figure
ax1=axes;

if ~exist('cline1','var') || strncmpi(cline1,'on' ,2) || strncmpi(cline1,'y',1)
    cline1='-';
elseif strncmpi(cline1,'off',3) || strncmpi(cline1,'n',1)
    cline1='none';
end
if ~exist('lwidth1','var')
    lwidth1=0.5;
end
if ~exist('cmark1','var') || strncmpi(cmark1,'off',3) || strncmpi(cmark1,'n',1)
    cmark1='none';
elseif strncmpi(cmark1,'on' ,2) || strncmpi(cmark1,'y',1) || ...
       (length(cmark1) > 1)
    cmark1='+';
end
if ~exist('msize1','var')
    msize1=6;
end

hl1=line(xdata,ydata,'Parent',ax1,'LineStyle',cline1,'LineWidth',lwidth1,...
         'Marker',cmark1,'MarkerSize',msize1);

%  set lines to have a continuous colormap

if exist('cmap','var')
    colormap(cmap)
end

cmap=colormap;
for i=1:length(hl1)
    if (length(hl1) > 1)
        imap=round((i-1)/(length(hl1)-1)*(size(cmap,1)-1))+1;
    else
        imap=1;
    end
    set(hl1(i),'Color',cmap(imap,:))
end

xlim('auto')
[xlims]=xlim;
if exist('xmin','var')
    xlims(1)=xmin;
end
if exist('xmax','var')
    xlims(2)=xmax;
end
xlim(xlims)
if exist('xgrid','var') && strncmpi(xgrid,'on',2)
    set(ax1,'XGrid','on')
end
if strncmpi(cplot,'p',1)
    [ylims]=[0 1];
else
    ylim('auto')
    [ylims]=ylim;
end
if exist('ymin1','var')
    ylims(1)=ymin1;
end
if exist('ymax1','var')
    ylims(2)=ymax1;
end
ylim(ylims)
if exist('ygrid1','var') && strncmpi(ygrid1,'on',2)
    set(ax1,'YGrid','on')
end

%  add the annotation

if strncmpi(cplot,'p',1) && ...
   exist('ynorm','var') && strncmpi(ynorm,'y',1)
%  copied and adapted from matlab normplot function
    if ~exist('yprob','var')
        yprob = [0.001 0.003 0.01 0.02 0.05 0.10 0.25 0.5...
                 0.75 0.90 0.95 0.98 0.99 0.997 0.999];
    end

%  matlab docs say str2mat is obsolete, so replace (and eliminate redundancy)
%     label1= str2mat('0.001','0.003', '0.01','0.02','0.05','0.10','0.25','0.50');
%     label2= str2mat('0.75','0.90','0.95','0.98','0.99','0.997', '0.999');
%     label = [label1;label2];
    label=cell(1,length(yprob));
    for i=1:length(yprob)
        label(i)=cellstr(num2str(yprob(i)));
    end

    tick  = norminv_issm(yprob,0,1);
    set(ax1,'YTick',tick,'YTickLabel',label);
    ylim([tick(1) tick(end)])
end

title([cplot 'ies for Response Levels'])
xlabel('Response Level')
ylabel([cplot 'y'])

hleg1=legend(ax1,descr(1:length(hl1)),'Location','NorthWest',...
             'Color','none','Interpreter','none');

%%  generate the probability distribution functions

if ~exist('pdfplt','var') || ~strcmpi(pdfplt,'off')
    xpdf =zeros(max(lcdfr)-1,length(dresp));
    xpdf (:,:)=NaN;
    ypdf =zeros(max(lcdfr)-1,length(dresp));
    ypdf (:,:)=NaN;
    xplot=zeros(2*max(lcdfr),length(dresp));
    xplot(:,:)=NaN;
    yplot=zeros(2*max(lcdfr),length(dresp));
    yplot(:,:)=NaN;

    for i=1:length(dresp)
        descr(i)=cellstr([dresp(i).descriptor ' pdf']);
        if ~isempty(dresp(i).cdf) && (lcdfr(i) > 1)
            xpdf (1:lcdfr(i)-1,i)=(dresp(i).cdf(1:lcdfr(i)-1,1)+...
                                   dresp(i).cdf(2:lcdfr(i)  ,1))/2;
            ypdf (1:lcdfr(i)-1,i)=(dresp(i).cdf(2:lcdfr(i)  ,2)-...
                                   dresp(i).cdf(1:lcdfr(i)-1,2))./...
                                  (dresp(i).cdf(2:lcdfr(i)  ,1)-...
                                   dresp(i).cdf(1:lcdfr(i)-1,1));
            for j=1:lcdfr(i)-1
                xplot(2*(j-1)+1,i)=dresp(i).cdf(j  ,1);
                yplot(2*(j-1)+1,i)=ypdf (j,i);
                xplot(2*(j-1)+2,i)=dresp(i).cdf(j+1,1);
                yplot(2*(j-1)+2,i)=ypdf (j,i);
            end
        end
    end

    if strcmpi(pdfplt,'line')
        xplot=xpdf;
        yplot=ypdf;
    end

%  draw the line plot

%  (see "Using Multiple X- and Y-Axes" and "Overlaying Other
%  Plots on Bar Graphs", or search on "YAxisLocation right")

%     hold all
%     hold on
%     plot(edges,cdf)
%     plotyy([],[],edges,cdf)

%  ticks from the bar plot will show through on the right side,
%  so make equal number of ticks for the line plot on right side

    nytick=length(get(ax1,'YTick'));
%     ylim('auto')
%     [ylims]=ylim;
    [ylims]=[0 ceil(max(max(yplot))/0.1-0.1)*0.1];
    if exist('ymin2','var')
        ylims(1)=ymin2;
    end
    if exist('ymax2','var')
        ylims(2)=ymax2;
    else
        ylims(2)=ylims(1)+(nytick-1)/(nytick-1-1)*(ylims(2)-ylims(1));
    end
%     ylim(ylims)
    ytinc =(ylims(2)-ylims(1))/(nytick-1);

    ax2=axes('Position',get(ax1,'Position'),...
             'XLim',get(ax1,'XLim'),...
             'XTick',get(ax1,'XTick'),...
             'YLim',ylims,...
             'YTick',[ylims(1):ytinc:ylims(2)],...
             'XAxisLocation','bottom','YAxisLocation','right',...
             'Color','none','Layer','top');

    if ~exist('cline2','var') || strncmpi(cline2,'on' ,2) || strncmpi(cline2,'y',1)
        cline2='-';
    elseif strncmpi(cline2,'off',3) || strncmpi(cline2,'n',1)
        cline2='none';
    end
    if ~exist('lwidth2','var')
        lwidth2=0.5;
    end
    if ~exist('cmark2','var') || strncmpi(cmark2,'off',3) || strncmpi(cmark2,'n',1)
        cmark2='none';
    elseif strncmpi(cmark2,'on' ,2) || strncmpi(cmark2,'y',1) || ...
           (length(cmark2) > 1)
        cmark2='+';
    end
    if ~exist('msize2','var')
        msize2=6;
    end

    hl2=line(xplot,yplot,'Parent',ax2,'LineStyle',cline2,'LineWidth',lwidth2,...
             'Marker',cmark2,'MarkerSize',msize2);

%  set line property colors line property

    cmap=colormap;
    for i=1:length(hl2)
        if (length(hl2) > 1)
            imap=round((i-1)/(length(hl2)-1)*(size(cmap,1)-1))+1;
        else
            imap=1;
        end
        set(hl2(i),'Color',get(hl1(i),'Color'))
     end

%  add the annotation

    ylabel('Probability Distribution Function');

    if ~exist('cdfleg','var') || ~strcmpi(cdfleg,'off')
% legend doesn't combine with bar chart above
        hleg2=legend(ax2,descr(1:length(hl2)),'Location','NorthEast',...
                     'Color','none','Interpreter','none');
%         set(hleg2,'Color','white')
    end

    set(gcf,'PaperPositionMode','auto')
%     hold off
end

end

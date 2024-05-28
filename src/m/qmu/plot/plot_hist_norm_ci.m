%
%  plot a relative histogram and cdf optionally along with a
%  normal distribution for the sample and confidence intervals.
%
%  []=plot_hist_norm_ci(dresp      ,params)
%  []=plot_hist_norm_ci(dresp,descr,params)
%  []=plot_hist_norm_ci(sampr,descr,params)
%
%  where the required input is:
%    dresp         (structure array, responses)
%      or
%    dresp         (structure array, responses)
%    descr         (cell array, list of response descriptions desired)
%      or
%    sampr         (double array, lists of samples)
%    descr         (cell array, list of descriptions)
%
%  the required fields of dresp are:
%    descriptor    (char, description)
%    sample        (double vector, list of samples)
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
%  and the optional input is:
%    hmin          (numeric, minimum for histogram)
%    hmax          (numeric, maximum for histogram)
%    hnint         (numeric, number of intervals for histogram)
%    ymin1         (numeric, minimum of histogram y-axis)
%    ymax1         (numeric, maximum of histogram y-axis)
%    ymin2         (numeric, minimum of cdf y-axis)
%    ymax2         (numeric, maximum of cdf y-axis)
%    nrmplt        (char, 'line' or 'off' to change nrm plots from 'bar')
%    ciplt         (char, 'line' or 'off' to change ci plots from 'bar')
%    cdfplt        (char, 'off' to turn off cdf line plots)
%    cdfleg        (char, 'off' to turn off cdf legends)
%    cmap          (char or numeric, colormap definition)
%
%  for each response in the input array, this function
%  calculates and plots a relative histogram and CDF of the list
%  of samples, and annotates it with the description.  in
%  addition, the normal distribution and CDF are plotted, and
%  four CDF's are plotted for the confidence intervals.
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
function []=plot_hist_norm_ci(varargin)

if ~nargin
    help plot_hist_norm_ci
    return
end

%%  process input data and assemble into matrices as needed

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
    lsamp=zeros(1,length(dresp));
    for i=1:length(dresp)
        lsamp(i)=length(dresp(i).sample);
    end
    sampr=zeros(max(lsamp),length(dresp));
    sampr(:,:)=NaN;

    mu     =zeros(1,length(dresp));
    sigma  =zeros(1,length(dresp));
    muci   =zeros(2,length(dresp));
    sigmaci=zeros(2,length(dresp));

    for i=1:length(dresp)
        descr(i)=cellstr(dresp(i).descriptor);
        sampr(1:lsamp(i),i)=dresp(i).sample;
        mu     (i)  =dresp(i).mean;
        sigma  (i)  =dresp(i).stddev;
        muci   (:,i)=dresp(i).meanci;
        sigmaci(:,i)=dresp(i).stddevci;
    end
else
    sampr=varargin{iarg};
    iarg=iarg+1;

    lsamp(1:size(sampr,2))=size(sampr,1);

    if     iarg <= nargin && iscell(varargin{iarg})
        descr=varargin{iarg};
        iarg=iarg+1;
%     elseif iarg <= nargin && ischar(varargin{iarg})
%         descr=cellstr(varargin{iarg});
%         iarg=iarg+1;
    else
        descr=cell(1,size(sampr,2));
    end

    mu     =zeros(1,size(sampr,2));
    sigma  =zeros(1,size(sampr,2));
    muci   =zeros(2,size(sampr,2));
    sigmaci=zeros(2,size(sampr,2));
    for i=1:size(sampr,2)
%  calculate 95% confidence intervals (same as dakota)
        [mu(i),sigma(i),muci(:,i),sigmaci(:,i)]=...
            normfit_issm(sampr(:,i),0.05);
    end
    display('Using calculated normal fits.')
end

for i=1:length(descr)
    if isempty(descr{i})
        descr(i)={['resp_' num2str(i)]};
    end
end

%  parameters

while (iarg <= nargin-1)
    if ischar(varargin{iarg})
        if ~isempty(strmatch(varargin{iarg},...
                {'hmin','hmax','hnint',...
                 'ymin1','ymax1','ymin2','ymax2',...
                 'nrmplt','ciplt','cdfplt','cdfleg','cmap'},...
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

%%  generate the intervals

if ~exist('hmin','var')
    hmin=min(min(sampr));
end
if ~exist('hmax','var')
    hmax=max(max(sampr));
end
if ~exist('hnint','var')
    hnint=50;
end
edges=hmin:(hmax-hmin)/hnint:hmax;

%%  generate the histogram counts and make them relative

%  note that for the histc function:
%  n(k) counts the value x(i) if edges(k) <= x(i) < edges(k+1).
%  The last bin counts any values of x that match edges(end).
%  Values outside the values in edges are not counted.
%  Use -inf and inf in edges to include all non-NaN values.

dhistc=histc(sampr,edges);
for i=1:size(sampr,2)
    dbelow(i)  =length(find(sampr(:,i)<edges(  1)))/lsamp(i);
    dhistc(:,i)=dhistc(:,i)                        /lsamp(i);
    dabove(i)  =length(find(sampr(:,i)>edges(end)))/lsamp(i);
end

ncol=size(sampr,2);
if exist('mu','var') && exist('sigma','var')
    for i=1:ncol
        dbelow(end+1)=normcdf_issm(edges(  1),mu(i),sigma(i));
        dhistc(1:size(dhistc,1)-1,end+1)=...
            normcdf_issm(edges(2:end  ),mu(i),sigma(i))-...
            normcdf_issm(edges(1:end-1),mu(i),sigma(i));
        dabove(end+1)=norminv_issm(edges(end),mu(i),sigma(i));
        if exist('descr','var')
            descr(end+1)={[descr{i} ' norm']};
        end
    end
end

if exist('muci','var') && exist('sigmaci','var')
    for i=1:ncol
        dbelow(end+1)=normcdf_issm(edges(  1),muci(1,i),sigmaci(2,i));
        dhistc(1:size(dhistc,1)-1,end+1)=...
            normcdf_issm(edges(2:end  ),muci(1,i),sigmaci(2,i))-...
            normcdf_issm(edges(1:end-1),muci(1,i),sigmaci(2,i));
        dabove(end+1)=norminv_issm(edges(end),muci(1,i),sigmaci(2,i));
        if exist('descr','var')
            descr(end+1)={[descr{i} ' norm-+']};
        end
    end
    for i=1:ncol
        dbelow(end+1)=normcdf_issm(edges(  1),muci(1,i),sigmaci(1,i));
        dhistc(1:size(dhistc,1)-1,end+1)=...
            normcdf_issm(edges(2:end  ),muci(1,i),sigmaci(1,i))-...
            normcdf_issm(edges(1:end-1),muci(1,i),sigmaci(1,i));
        dabove(end+1)=norminv_issm(edges(end),muci(1,i),sigmaci(1,i));
        if exist('descr','var')
            descr(end+1)={[descr{i} ' norm--']};
        end
    end
    for i=1:ncol
        dbelow(end+1)=normcdf_issm(edges(  1),muci(2,i),sigmaci(1,i));
        dhistc(1:size(dhistc,1)-1,end+1)=...
            normcdf_issm(edges(2:end  ),muci(2,i),sigmaci(1,i))-...
            normcdf_issm(edges(1:end-1),muci(2,i),sigmaci(1,i));
        dabove(end+1)=norminv_issm(edges(end),muci(2,i),sigmaci(1,i));
        if exist('descr','var')
            descr(end+1)={[descr{i} ' norm+-']};
        end
    end
    for i=1:ncol
        dbelow(end+1)=normcdf_issm(edges(  1),muci(2,i),sigmaci(2,i));
        dhistc(1:size(dhistc,1)-1,end+1)=...
            normcdf_issm(edges(2:end  ),muci(2,i),sigmaci(2,i))-...
            normcdf_issm(edges(1:end-1),muci(2,i),sigmaci(2,i));
        dabove(end+1)=norminv_issm(edges(end),muci(2,i),sigmaci(2,i));
        if exist('descr','var')
            descr(end+1)={[descr{i} ' norm++']};
        end
    end
end

%  draw the bar plot

figure
if ~exist('nrmplt','var')
    nrmplt='bar';
end
if ~exist('ciplt','var')
    ciplt='bar';
end

hold all
% hl1=bar(edges(1:end-1),dhistc(1:end-1,1:6*ncol));
% hl1=line(edges(1:end-1),dhistc(1:end-1,1:6*ncol));
if     strncmpi(nrmplt,'b',1)
    if     strncmpi(ciplt,'b',1)
        hl1=bar (edges(1:end-1),dhistc(1:end-1,1:6*ncol));
    elseif strncmpi(ciplt,'l',1)
        hl1=bar (edges(1:end-1),dhistc(1:end-1,1:2*ncol));
        hl1(2*ncol+1:6*ncol)=...
            line(edges(1:end-1),dhistc(1:end-1,2*ncol+1:6*ncol),...
                 'LineWidth',2);
    elseif strncmpi(ciplt,'off',3) || strncmpi(ciplt,'n',1)
        hl1=bar (edges(1:end-1),dhistc(1:end-1,1:2*ncol));
    end
elseif strncmpi(nrmplt,'l',1)
    if     strncmpi(ciplt,'b',1) || strncmpi(ciplt,'l',1)
        hl1=bar (edges(1:end-1),dhistc(1:end-1,1:1*ncol));
        hl1(1*ncol+1:2*ncol)=...
            line(edges(1:end-1),dhistc(1:end-1,1*ncol+1:2*ncol),...
                 'LineWidth',2);
        hl1(2*ncol+1:6*ncol)=...
            line(edges(1:end-1),dhistc(1:end-1,2*ncol+1:6*ncol),...
                 'LineWidth',1);
    elseif strncmpi(ciplt,'off',3) || strncmpi(ciplt,'n',1)
        hl1=bar (edges(1:end-1),dhistc(1:end-1,1:1*ncol));
        hl1(1*ncol+1:2*ncol)=...
            line(edges(1:end-1),dhistc(1:end-1,1*ncol+1:2*ncol),...
                 'LineWidth',2);
    end
elseif strncmpi(nrmplt,'off',3) || strncmpi(nrmplt,'n',1)
    hl1=bar (edges(1:end-1),dhistc(1:end-1,1:1*ncol));
end
ax1=gca;
hold off

%  set barseries properties for clarity

if (ncol > 1) || strncmpi(nrmplt,'b',1)
    for i=1:length(hl1)
        if strcmpi(get(hl1(i),'Type'),'hggroup')
            set(hl1(i),'BarWidth',1,'EdgeColor','none');
        end
    end
end

%  set bars and lines to have a continuous colormap
%  (if barseries is "flat", must interpolate and round to colormap)

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
    if     strcmpi(get(hl1(i),'Type'),'hggroup')
        if ischar(get(hl1(i),'FaceColor')) && ...
           strcmpi(get(hl1(i),'FaceColor'),'flat')
            set(hl1(i),'FaceColor',cmap(imap,:))
        else
            set(hl1(i),'FaceColor',get(hl1(i),'FaceColor'))
        end
    elseif strcmpi(get(hl1(i),'Type'),'line')
        set(hl1(i),'Color',cmap(imap,:))
    end
end

xlim('auto')
[xlims]=xlim;
if exist('hmin','var')
    xlims(1)=edges(1);
end
if exist('hmax','var')
    xlims(2)=edges(end-1);
end
xlim(xlims)

ylim('auto')
[ylims]=ylim;
if exist('ymin1','var')
    ylims(1)=ymin1;
end
if exist('ymax1','var')
    ylims(2)=ymax1;
end
ylim(ylims)

%  add the annotation

if exist('cdfplt','var') && ...
   (strncmpi(cdfplt,'off',3) || strncmpi(cdfplt,'n',1))
    title('Relative Frequency Histogram')
else
    title('Relative Frequency Histogram with CDF')
end
xlabel('Interval Edge Value');
ylabel('Relative Frequency');

if exist('descr','var')
    hleg1=legend(ax1,descr(1:length(hl1)),'Location','NorthWest',...
                 'Color','none','Interpreter','none');
else
    hleg1=legend(ax1);
end

%%  generate the cumulative distribution functions

if ~exist('cdfplt','var') || ...
   (~strncmpi(cdfplt,'off',3) && ~strncmpi(cdfplt,'n',1))
%     cdf=zeros(size(dhistc));
%     cdf(1,:)=dhistc(1,:);
%     for i=2:size(dhistc,1)
%         cdf(i,:)=cdf(i-1,:)+dhistc(i,:);
%     end
    cdf=cumsum(dhistc);
    for i=1:size(dhistc,2)
        cdf(:,i)=dbelow(i)+cdf(:,i);
    end
    if exist('descr','var')
        ncol=length(descr);
        for i=1:ncol
            cdescr(i)={[descr{i} ' cdf']};
        end
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
    [ylims]=[0 ceil(max(max(cdf))/0.1-0.1)*0.1];
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
    hl2=line(edges(1:end-1),cdf(1:end-1,1:length(hl1)),'Parent',ax2);

%  set line property colors to match barseries or line property
%  (if barseries is "flat", must interpolate and round to colormap)

    cmap=colormap;
    for i=1:length(hl2)
        if (length(hl2) > 1)
            imap=round((i-1)/(length(hl2)-1)*(size(cmap,1)-1))+1;
        else
            imap=1;
        end
        if     strcmpi(get(hl1(i),'Type'),'hggroup')
            if ischar(get(hl1(i),'FaceColor')) && ...
               strcmpi(get(hl1(i),'FaceColor'),'flat')
                set(hl2(i),'Color',cmap(imap,:))
            else
                set(hl2(i),'Color',get(hl1(i),'FaceColor'))
            end
        elseif strcmpi(get(hl1(i),'Type'),'line')
            set(hl2(i),'Color',get(hl1(i),'Color'))
        end
    end

%  add the annotation

    ylabel('Cumulative Percent');

    if ~exist('cdfleg','var') || ...
       (~strncmpi(cdfleg,'off',3) && ~strncmpi(cdfleg,'n',1))
% legend doesn't combine with bar chart above
        if exist('cdescr','var')
            hleg2=legend(ax2,cdescr(1:length(hl2)),'Location','NorthEast',...
                         'Color','none','Interpreter','none');
%             set(hleg2,'Color','white')
        else
            hleg2=legend(ax2);
%             set(hleg2,'Color','white')
        end
    end

    set(gcf,'PaperPositionMode','auto')
%     hold off
end

end

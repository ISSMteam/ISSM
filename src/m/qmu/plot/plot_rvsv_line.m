%
%  plot line plots of responses vs. variables.
%
%  []=plot_rvsv_line(dvar       ,dresp      ,params)
%  []=plot_rvsv_line(dvar ,descv,dresp,descr,params)
%  []=plot_rvsv_line(sampv,descv,sampr,descr,params)
%
%  where the required input is:
%    dvar          (structure array, variables)
%      or
%    dvar          (structure array, variables)
%    descv         (cell array, list of variable descriptions desired)
%      or
%    sampv         (double array, lists of variable samples)
%    descv         (cell array, list of variable descriptions)
%
%    dresp         (structure array, responses)
%      or
%    dresp         (structure array, responses)
%    descr         (cell array, list of response descriptions desired)
%      or
%    sampr         (double array, lists of response samples)
%    descr         (cell array, list of response descriptions)
%
%  the required fields of dvar and dresp are:
%    descriptor    (char, description)
%    sample        (double vector, list of samples)
%
%  the optional input is:
%    params        (string/numeric, parameter names and values)
%
%  where the optional parameters are:
%    nplotr        (numeric, number of plot rows)
%    nplotc        (numeric, number of plot columns)
%    ymin          (numeric, minimum of y-axis)
%    ymax          (numeric, maximum of y-axis)
%    yscat         (char, 'off' to turn off y-axis scattergram)
%
%  for each variable/response combination in the input array, this
%  function plots a line plot.  all of the variables and responses
%  are plotted on the same axes, if nplotr and nplotc are not
%  specified, so some scaling might otherwise be desired.
%
%  dvar and dresp data would typically be contained in the dakota
%  tabular output file from a sampling or parametric analysis, and
%  read by dakota_out_parse.
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
function []=plot_rvsv_line(varargin)

if ~nargin
    help plot_rvsv_line
    return
end

%%  process input data and assemble into matrices as needed

%  variables

iarg=1;
if isstruct(varargin{iarg})
    dvar=varargin{iarg};
    iarg=iarg+1;

%     if iarg <= nargin && (iscell(varargin{iarg}) || ischar(varargin{iarg}))
    if iarg <= nargin && iscell(varargin{iarg})
        dvar=struc_desc(dvar,varargin{iarg});
        iarg=iarg+1;
    end

    descv=cell (1,length(dvar));
    lsamp=zeros(1,length(dvar));
    for i=1:length(dvar)
        lsamp(i)=length(dvar(i).sample);
    end
    sampv=zeros(max(lsamp),length(dvar));
    sampv(:,:)=NaN;

    for i=1:length(dvar)
        descv(i)=cellstr(dvar(i).descriptor);
        sampv(1:lsamp(i),i)=dvar(i).sample;
    end
else
    sampv=varargin{iarg};
    iarg=iarg+1;

    if     iarg <= nargin && iscell(varargin{iarg})
        descv=varargin{iarg};
        iarg=iarg+1;
%     elseif iarg <= nargin && ischar(varargin{iarg})
%         descv=cellstr(varargin{iarg});
%         iarg=iarg+1;
    else
        descv=cell(1,size(sampv,2));
    end
end

for i=1:length(descv)
    if isempty(descv{i})
        descv(i)={['var_' i]};
    end
end

%  responses

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

    for i=1:length(dresp)
        descr(i)=cellstr(dresp(i).descriptor);
        sampr(1:lsamp(i),i)=dresp(i).sample;
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
        descr=cell(1,size(sampr,2));
    end
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
                {'nplotr','nplotc',...
                 'ymin','ymax','yscat'},...
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

if     ~exist('nplotr','var') && ~exist('nplotc','var')
    nplotr=1;
    nplotc=1;
elseif ~exist('nplotr','var')
    nplotr=ceil(size(sampr,2)*size(sampv,2)/nplotc);
elseif ~exist('nplotc','var')
    nplotc=ceil(size(sampr,2)*size(sampv,2)/nplotr);
end

%%  filter, sort, and plot the data

%  while it would be preferable for the outer loop to be responses,
%  it is more efficient for the outer loop to be variables

figure
haxes=[];
hplot=[];
cdesc={};
hscat=[];

iplot=0;

for ivar=1:size(sampv,2)
    [vval,indxv,indxvi]=unique(sampv(:,ivar),'first');
    indxv2=setdiff(1:size(sampv,1),indxv);

    for iresp=1:size(sampr,2)

%  initialize the subplot

        if (ivar*iresp == 1) || ...
           (nplotr*nplotc > 1)
            iplot=iplot+1;
            haxes(end+1)=subplot(nplotr,nplotc,iplot);
            hold all
        end

        hplot(end+1)=plot   (sampv(indxv ,ivar),sampr(indxv ,iresp),'-x');
        cdesc(end+1)={[descr{iresp} ' wrt ' descv{ivar}]};
        if ~exist('yscat','var') || ...
           (~strncmpi(yscat,'off',3) && ~strncmpi(yscat,'n',1))
            hscat(end+1)=scatter(sampv(indxv2,ivar),sampr(indxv2,iresp),'+k');
%  see "controlling legends" in Matlab on-line docs
%         cdesc(end+1)={['constant ' descv{ivar}]};
            set(get(get(hscat(end),'Annotation'),'LegendInformation'),...
                'IconDisplayStyle','off'); % Exclude line from legend
        end

%  add the annotation

        if (ivar*iresp == size(sampv,2)*size(sampr,2)) || ...
           (nplotr*nplotc > 1)
            hold off

            ylim('auto')
            [ylims]=ylim;
            if exist('ymin','var')
                ylims(1)=ymin;
            end
            if exist('ymax','var')
                ylims(2)=ymax;
            end
            ylim(ylims)

            if (size(sampv,2) == 1) || (nplotr*nplotc > 1)
                xlabc=descv{ivar};
            else
                xlabc='Variables';
            end
            if (size(sampr,2) == 1) || (nplotr*nplotc > 1)
                ylabc=descr{iresp};
            else
                ylabc='Responses';
            end
            title([ylabc ' vs. ' xlabc],'Interpreter','none');
            xlabel(xlabc,'Interpreter','none');
            ylabel(ylabc,'Interpreter','none');

            if (nplotr*nplotc == 1) && (size(sampv,2)*size(sampr,2) > 1)
                legend(cdesc,'Location','EastOutside','Interpreter','none');
            end
        end
    end
end

end

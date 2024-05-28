%
%  plot surface plots of variables vs. responses.
%
%  []=plot_rvsv_surf(dvar       ,dresp      ,params)
%  []=plot_rvsv_surf(dvar ,descv,dresp,descr,params)
%  []=plot_rvsv_surf(sampv,descv,sampr,descr,params)
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
%    zmin          (numeric, minimum of z-axis)
%    zmax          (numeric, maximum of z-axis)
%    cmin          (numeric, minimum of colorbar)
%    cmax          (numeric, maximum of colorbar)
%
%  for each response in the input array, this function plots a
%  surface plot.  there should be two and only two variables.
%  each response will be in a separate surface plot; hence the
%  need for nplotr and nplotc.
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
function []=plot_rvsv_surf(varargin)

if ~nargin
    help plot_rvsv_surf
    return
end

%%  process input data and assemble into matrices as needed

iarg=1;

%  variables

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

if (size(sampv,2) ~= 2)
    error('Two and only two variables required for surface plot.');
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
                 'zmin','zmax',...
                 'cmin','cmax'},...
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
    nplotr=ceil(sqrt(size(sampr,2)));
    nplotc=ceil(size(sampr,2)/nplotr);
elseif ~exist('nplotr','var')
    nplotr=ceil(size(sampr,2)/nplotc);
elseif ~exist('nplotc','var')
    nplotc=ceil(size(sampr,2)/nplotr);
end

%%  filter, sort, and plot the data

figure
haxes=[];
hsurf=[];

[x,ix,ixi]=unique(sampv(:,1),'first');
[y,iy,iyi]=unique(sampv(:,2),'first');

for iresp=1:size(sampr,2)
    z=zeros(length(x),length(y));
    for i=1:size(sampr,1)
        z(ixi(i),iyi(i))=sampr(i,iresp);
    end

%  initialize the subplot

    haxes(iresp)=subplot(nplotr,nplotc,iresp);
%     hsurf(iresp)=surfc(x,y,z,'FaceColor','interp');
    surfc(x,y,z,'FaceColor','interp');

    zlim('auto')
    [zlims]=zlim;
    if exist('zmin','var')
        zlims(1)=zmin;
    end
    if exist('zmax','var')
        zlims(2)=zmax;
    end
    zlim(zlims)

%  add the annotation

    title([descr{iresp} ' wrt ' descv{1} ' and ' descv{2}],...
          'Interpreter','none');
    xlabel(descv{1},'Interpreter','none');
    ylabel(descv{2},'Interpreter','none');
    zlabel(descr{iresp},'Interpreter','none');

    caxis('auto')
    [cmini,cmaxi]=caxis;
    if exist('cmin','var')
        cmini=cmin;
    end
    if exist('cmax','var')
        cmaxi=cmax;
    end
    caxis([cmini cmaxi])

    colorbar
end

end

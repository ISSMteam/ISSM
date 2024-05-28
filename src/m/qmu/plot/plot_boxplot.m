%
%  plot a box plot of the responses.
%
%  []=plot_boxplot(dresp      ,params)
%  []=plot_boxplot(dresp,descr,params)
%  []=plot_boxplot(sampr,descr,params)
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
%  the optional input is:
%    params        (string/numeric, parameter names and values)
%
%  where the optional parameters are:
%    ymin          (numeric, minimum of y-axis)
%    ymax          (numeric, maximum of y-axis)
%
%  for each response in the input array, this function plots a
%  matlab box plot of the list of samples and annotates it
%  with the description.  the lists of samples need not all be
%  the same length.
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
function []=plot_boxplot(varargin)

if ~nargin
    help plot_boxplot
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
                {'ymin','ymax'},...
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

%%  draw the plot

%  draw box plot

figure
boxplot(sampr,'labels',descr,'notch','on')
ax1=gca;

%  add the annotation

ylim('auto')
[ylims]=ylim;
if exist('ymin','var')
    ylims(1)=ymin;
end
if exist('ymax','var')
    ylims(2)=ymax;
end
ylim(ylims)

if (size(sampr,2) == 1)
    tlabc=descr{1};
else
    tlabc='Responses';
end
title(['Box Plot of ' tlabc],'Interpreter','none');
xlabel('Response','Interpreter','none');
ylabel('Value'   ,'Interpreter','none');

end

function plotmodel(md,varargin)
%At command prompt, type plotdoc for help on documentation.

%First process options
options=plotoptions(varargin{:});

%Get figure number and number of plots
figurenumber=options.figurenumber;
numberofplots=options.numberofplots;

%get number of subplots
subplotwidth=ceil(sqrt(numberofplots));

%if nlines and ncols specified, then bypass.
if ~exist(options.list{1},'nlines') & ~exist(options.list{1},'ncols')
	 nlines=ceil(numberofplots/subplotwidth);
	 ncols=subplotwidth;
else
	nlines=getfieldvalue(options.list{1},'nlines',NaN);
	ncols=getfieldvalue(options.list{1},'ncols',NaN);
	if isnan(nlines), nlines = ceil(numberofplots/ncols);  end
	if isnan(ncols),  ncols  = ceil(numberofplots/nlines); end
end

%check that nlines and ncols were given at the same time!
if ((exist(options.list{1},'ncols') & ~exist(options.list{1},'ncols')) | (~exist(options.list{1},'ncols') & exist(options.list{1},'ncols')))
	error('plotmodel error message: nlines and ncols  need to be specified together, or not at all');
end

% Add option for subplot, only support one subplot per one plotmodel call
subindex = 0;
if (exist(options.list{1},'subplot'))
    subops = getfieldvalue(options.list{1},'subplot',NaN);
    nlines = subops(1);
    ncols = subops(2);
    subindex = subops(3);
    % change the width
    subplotwidth = ncols;
    % Turn off new figure option for subplot
    options.list{1}=addfield(options.list{1},'figurestatement','off');
    % Turn off clf option for subplot
    options.list{1}=addfield(options.list{1},'clf','off');
end

%go through subplots
if numberofplots,

	%Create figure 
	if (strcmpi(getfieldvalue(options.list{1},'figurestatement','on'),'on'))
		f=figure(figurenumber); 
	else
		f=gcf;
	end
	if strcmpi(getfieldvalue(options.list{1},'clf','on'),'on'),
		clf;
	end
	if strcmpi(getfieldvalue(options.list{1},'visible','on'),'off'),
		set(f,'Visible','Off');
	end

	if exist(options.list{1},'figposition'), % {{{
		figposition=getfieldvalue(options.list{1},'figposition');
		if ischar(figposition),
			if strcmpi(figposition,'larour'),
				set(gcf,'Position',[1604 4 1594 1177]);
			elseif strcmpi(figposition,'larour2'),
				set(gcf,'Position',[756    62   827   504]);
			elseif strcmpi(figposition,'mathieu'),
				set(gcf,'Position',[300 1 1580 1150]);
			elseif strcmpi(figposition,'byron'),
				set(gcf,'Position',[40 1580 560*1.25 420*1.25]);  %left, bottom, width, height; W=560 H=420 is default
			elseif strcmpi(figposition,'fullscreen'),
				set(gcf,'Position',get(0,'ScreenSize'));
			elseif strcmpi(figposition,'halfright'),
				screen=get(0,'ScreenSize');
				left=screen(1); bott=screen(2); widt=screen(3); heig=screen(4)-25;
				set(gcf,'Position',fix([left+widt/2 bott widt/2 heig]));
			elseif strcmpi(figposition,'halfleft'),
				screen=get(0,'ScreenSize');
				left=screen(1); bott=screen(2); widt=screen(3); heig=screen(4)-25;
				set(gcf,'Position',fix([left bott widt/2 heig]));
			elseif strcmpi(figposition,'square'),
				screen=get(0,'ScreenSize');
				left=screen(1); bott=screen(2); widt=min(screen(3)-25,screen(4)-25);
				set(gcf,'Position',fix([left+(screen(3)-widt) bott widt widt]));
			elseif strcmpi(figposition,'portrait'),
				%reformat with letter paper size (8.5" x 11")
				screen=get(0,'ScreenSize');
				left=screen(1); bott=screen(2); widt=screen(3); heig=screen(4)-25;
				portrait=fix([left+widt-(heig*8.5/11) bott heig*8.5/11 heig]);
				set(gcf,'Position',portrait)
			elseif strcmpi(figposition,'landscape'),
				%reformat with letter paper size (8.5" x 11")
				screen=get(0,'ScreenSize');
				left=screen(1); bott=screen(2); widt=screen(3); heig=screen(4)-25;
				landscape=fix([left+widt-(heig*11/8.5) bott heig*11/8.5 heig]);
				set(gcf,'Position',landscape)
			else
				disp('''figposition'' string not supported yet');
			end
		else
			set(gcf,'Position',figposition);
		end
	elseif strcmpi(getenv('USER'),'inwoo')
		set(gcf,'Position',[910 242 560 420]);
	end % }}}

	%Use zbuffer renderer (smoother colors) and white background
	set(f,'Renderer','zbuffer','color',getfieldvalue(options.list{1},'figurebackgroundcolor','w'));

	%Go through all data plottable and close window if an error occurs
	try,
		for i=1:numberofplots,
            if subindex
                plot_manager(getfieldvalue(options.list{i},'model',md),options.list{i},subplotwidth,nlines,ncols,subindex);
            else
                plot_manager(getfieldvalue(options.list{i},'model',md),options.list{i},subplotwidth,nlines,ncols,i);
            end
            %List all unused options
			if getfieldvalue(options.list{i},'displayunused',1),
				displayunused(options.list{i})
			end
		end
	catch me,
		%figure(figurenumber),close;
		rethrow(me);
	end
else
	error('plotmodel error message: no output data found.');
end

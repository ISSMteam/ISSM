function scaleruler(options)
%SCALERULER - overlay a scale ruler on current plot
%
%   Usage:
%      scaleruler(options)

%get options
structure  = getfieldvalue(options,'scaleruler');
fontcolor  = getfieldvalue(options,'fontcolor','k');
fontweight = getfieldvalue(options,'fontweight','n');
fontsize   = getfieldvalue(options,'scaleruler_fontsize',16);
unitscale  = getfieldvalue(options,'unit',1.);
%fontcolor  = 'w';

%Go through structure and fill missing arguments
if length(structure)~=5
	error('plotmodel error message: bad number of input arguments for scaleruler: [x0 y0 length thickness numberofticks]');
end

%retrieve scale parameters
x0            = double(structure(1))/unitscale;
y0            = double(structure(2))/unitscale;
lengthscale   = double(structure(3))/unitscale;
widthscale    = double(structure(4))/unitscale;
numberofticks = double(structure(5));

%If only one tick, just draw a rectangle
if numberofticks==1,
	t=text(x0+lengthscale/2,y0+2*widthscale,2,[num2str(lengthscale*unitscale/1000) ' km'],...
		'FontSize',fontsize,'FontWeight',fontweight,'Color',fontcolor,'HorizontalAlignment','center','VerticalAlignment','bottom');
	if ~verLessThan('matlab', '8.4')
		set(t,'Layer','front');
	end
	p=patch([x0 x0+lengthscale x0+lengthscale x0],[y0 y0 y0+widthscale y0+widthscale],2*ones(1,4),fontcolor,'Edgecolor',fontcolor);
else
	%initialize some coordinates
	unitlength=lengthscale/(numberofticks -1);
	flag=-1;

	Bd=[x0 y0];
	Bu=[x0 y0+widthscale];
	Tick=0;

	%Text
	xt=Bu(1);
	yt=Bu(2)+widthscale;
	text(xt,yt,2,num2str(Tick),'FontSize',fontsize,'FontWeight',fontweight,'Color',fontcolor,'HorizontalAlignment','center','VerticalAlignment','bottom');

	%loope over the patches
	for i=1:numberofticks-1,
		Au=Bu;
		Ad=Bd;
		Bu=[Au(1)+unitlength Ad(2)+widthscale];
		Bd=[Ad(1)+unitlength Ad(2)];
		Tick=(Tick+unitlength)*unitscale;

		%pathes
		if flag==-1
			p=patch([Ad(1) Bd(1) Bu(1) Au(1)],[Ad(2) Bd(2) Bu(2) Au(2)],2*ones(1,4),'Black');
		else
			p=patch([Ad(1) Bd(1) Bu(1) Au(1)],[Ad(2) Bd(2) Bu(2) Au(2)],2*ones(1,4),'White');
		end

		%flip flag
		flag=-flag;

		%Text
		xt=Bu(1);
		yt=Bu(2)+widthscale;
		if i~=numberofticks-1,
			text(xt,yt,2,num2str(round_ice(Tick/1000,3)),'FontSize',fontsize,'FontWeight',fontweight,'Color',fontcolor,'HorizontalAlignment','center','VerticalAlignment','bottom');
		end
	end
	text(xt,yt,2,num2str(round_ice(Tick/1000,3)),'FontSize',fontsize,'FontWeight',fontweight,'Color',fontcolor,'HorizontalAlignment','center','VerticalAlignment','bottom');
	% add leading spaces depending on length of label string
	str=' km';
	for i=1:numel(num2str(round_ice(Tick/1000,3))),
		str=[' ' str];
	end
	text(xt,yt,2,str,'FontSize',fontsize,'FontWeight',fontweight,'Color',fontcolor,'HorizontalAlignment','left','VerticalAlignment','bottom');
end

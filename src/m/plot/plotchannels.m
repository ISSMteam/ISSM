function plotchannels(md,channelarea,varargin)
%PLOTCHANNELS - plot GlaDS channel area
%
%   Usage:
%      plotchannels(md,channelarea,options)
%      plotchannels(md,channeldischarge,options)
%
%   List of supported options
%      - 'min' minimum channel area displayed (default is max(channelarea))
%      - 'max' maximum channel area displayed (default is min(channelarea))
%      - 'colormap' colormap used (default is 'gray')
%      - 'linewidth' linewidth (default is 2)
%      - 'quiver' only used for discharge(default is 2)
%
%   Example:
%      plotchannels(md,md.results.TransientSolution(end).ChannelArea,'min',15,'max',36);
%      plotchannels(md,md.results.TransientSolution(end).ChannelDischarge,'min',15,'max',36,'quiver',1);

%Process options
options = pairoptions(varargin{:});

%is quiver?
isquiver = (getfieldvalue(options,'quiver',0)==1);
linewidth = getfieldvalue(options,'linewidth',2);

%define level
if isquiver
	level = abs(channelarea);
else
	level = channelarea;
end

%Some processing
Min = getfieldvalue(options,'min',min(level));
Max = getfieldvalue(options,'max',max(level));

%Create colormap
options=addfielddefault(options,'colormap',flipud(gray()));
palette = getcolormap(options);
numcolors=size(palette,1);
levels=round_ice(linspace(Min,Max,numcolors+1),2);

colorind=zeros(length(level),1);
for i=1:numcolors
	pos=find((level>=levels(i)) & (level<=levels(i+1)) );
	colorind(pos)=i;
end
colorind(find(level>levels(end)))=numcolors;

%Reconstruct edges
% {{{
tic
%Maximum number of edges
maxnbf = 3*md.mesh.numberofelements;
%Initialize intermediaries
edges = zeros(maxnbf,3);
exchange = zeros(maxnbf,1);
%Chaining algorithm
head_minv = -1*ones(md.mesh.numberofvertices,1);
next_face = zeros(maxnbf,1);
%Initialize number of faces
nbf = 0;
for i=1:md.mesh.numberofelements
	for j=1:3
		%Get the two indices of the edge number j of the ith triangle
		v1 = md.mesh.elements(i,j);
		if(j==3)
			v2 = md.mesh.elements(i,1);
		else
			v2 = md.mesh.elements(i,j+1);
		end
		%sort
		if(v2<v1)
			v3=v2; v2=v1; v1=v3;
		end
		%This edge a priori has not been processed yet
		exists = false;
		%Go through all processed faces connected to v1 and check whether we have seen this edge yet
		e=head_minv(v1);
		while(e~=-1)
			if(edges(e,2)==v2)
				exists = true;
				break;
			end
			e=next_face(e);
		end
		%If this edge is new, add it to the lists
		if(~exists)
			%Update edges
			edges(nbf+1,1) = v1; %vertex 1
			edges(nbf+1,2) = v2; %vertex 2
			edges(nbf+1,3) = i;  %element 1 (ignore second one)
			if(v1~=md.mesh.elements(i,j)) exchange(nbf+1)=1; end
			%Update chain
			next_face(nbf+1) = head_minv(v1);
			head_minv(v1)    = nbf+1;
			%Increase number of faces
			nbf=nbf+1;
		end
	end
end

edges = edges(1:nbf,:);
pos = find(exchange);
v3 = edges(pos,1);
edges(pos,1) = edges(pos,2);
edges(pos,2) = v3;
toc
% }}}

%Change edges formatting so that plot looks ok
myedges = [edges(:,1:2) edges(:,1)]';

%Loop over all levels and plot
for i=1:numcolors
	pos=find(colorind==i);
	hprime=plot(md.mesh.x(myedges(:,pos)),md.mesh.y(myedges(:,pos)),'-','Color',palette(i,:),'LineWidth',linewidth);
	if i==1; hold on; end

	if isquiver
		x1 = md.mesh.x(edges(pos,1)); x2 = md.mesh.x(edges(pos,2));
		y1 = md.mesh.y(edges(pos,1)); y2 = md.mesh.y(edges(pos,2));
		len = sqrt( (x2-x1).^2 + (y2-y1).^2);

		xq = mean([x1 x2],2);
		yq = mean([y1 y2],2);
		tx = sign(channelarea(pos)).*(x2-x1)./len .* len/10;
		ty = sign(channelarea(pos)).*(y2-y1)./len .* len/10;
		px = -ty;
		py = tx;

		%hprime=quiver(xq,yq,tx,ty,'color',palette(i,:),'showarrowhead','on','autoscale','off');
		num = length(pos);
		patch('Faces',reshape(1:3*num,num,3),'Vertices',[[xq;xq+(-tx+px);xq+(-tx-px)],[yq;yq+(-ty+py);yq+(-ty-py)]],'FaceColor',palette(i,:),'EdgeColor','none');
	end
end

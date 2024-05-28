function plot_qmustddev(md,options,nlines,ncols,i)
%PLOT_QMUMEAN - plot stddev of a scaled response 
%
%   Usage:
%      plot_qmustddev(md,options,nlines,ncols,i);
%
%   See also: PLOTMODEL

%plot mesh
subplot(nlines,ncols,i); 

%edgecolor
edgecolor=getfieldvalue(options,'edgecolor','none');

%process data and model
[x y z elements is2d isplanet]=processmesh(md,[],options);

%find response function
if exist(options,'qmudata'), 
	descriptor=getfieldvalue(options,'qmudata'); 
	if ~ischar(descriptor),
		error('plot_qmustddev error message:  descriptor should be a string');
	end
else 
	error('plot_qmustddev error message:  provide descriptor of response function in ''qmudata'' option');
end

%go pick up the response: 
allresponses=md.qmu.results.dresp_out;
responses=zeros(md.qmu.numberofpartitions,1);

count=1;
for i=1:length(allresponses),
	d=allresponses(i).descriptor;
	if strncmpi(d,'scaled_',7),
		d=d(8:end);
		if strncmpi(d,descriptor,length(descriptor)),
			responses(count)=allresponses(i).stddev/allresponses(i).mean*100;
			count=count+1;
		end
	end
end

%log?
if exist(options,'log'),
	responses=log(responses)/log(getfieldvalue(options,'log'));
end

%now, project onto vertices
responses_on_node=responses(md.qmu.partition+1);

%plot
A=elements(:,1); B=elements(:,2); C=elements(:,3); 
patch( 'Faces', [A B C], 'Vertices', [x y z],'FaceVertexCData', responses_on_node,'FaceColor','interp','EdgeColor',edgecolor);

%apply options
options=addfielddefault(options,'title',['Stddev  distribution of ' descriptor ' in %']);
options=addfielddefault(options,'colorbar',0);
applyoptions(md,[],options);

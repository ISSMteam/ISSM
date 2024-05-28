function plot_importancefactors(md,options,width,ii)
%PLOT_IMPORTANCEFACTORS - plot importance factors
%
%   Usage:
%      plot_importancefactors(md,options,width,i);
%
%   See also: PLOTMODEL

%first recover design variable descriptor
if exist(options,'designvariable'),
	descriptor=getfieldvalue(options,'designvariable');
else
	error('plot_importancefactors error message: Need to supply design variable descriptor');
end
descriptorlength=length(descriptor);

%then recover responsfunction name
if exist(options,'responsefunction'),
	responsefunctiondescriptor=getfieldvalue(options,'responsefunction');
else
	error('plot_importancefactors error message: Need to supply response function descriptor');
end

%go through all response functions and find the one corresponding to the correct responsefunctiondescriptor
responsefunctions=md.qmu.results{2};
found=0;
for i=1:length(responsefunctions),
	if strcmpi(responsefunctions(i).descriptor,responsefunctiondescriptor),
		found=i;
		break;
	end
end
if ~found,
	error('plot_importancefactors error message: could not find correct response function');
end
responsefunctions=responsefunctions(found);
nfun=size(responsefunctions.desvar,1);

%Now recover response to the correct desgin variable
importancefactors=zeros(md.qmu.numberofpartitions,1);
count=0;
for i=1:nfun,
	desvar=responsefunctions.desvar{i};
	if strncmpi(desvar,descriptor,descriptorlength),
		count=count+1;
		importancefactors(count)=responsefunctions.impfac(i);
	end
end
if count==0,
	error('plot_importancefactors error message: could not find to response functions with corresponding design variable');
end

%log?
if exist(options,'log'),
	logvalue=getfieldvalue(options,'log');
	importancefactors=log(importancefactors)/log(logvalue);
end

%Ok, get partitioning.
[epart npart]=MeshPartition(md,md.qmu.numberofpartitions);

%distribute importance factor
nodeimportance=importancefactors(npart);

%process data and model
[x y z elements is2d isplanet]=processmesh(md,[],options);

%edgecolor
edgecolor=getfieldvalue(options,'edgecolor','none');

%standard plot:
subplot(width,width,ii);

%ok, plot nodeimportance now.
if is2d,
	A=elements(:,1); B=elements(:,2); C=elements(:,3); 
	patch( 'Faces', [A B C], 'Vertices', [x y z],'FaceVertexCData', nodeimportance,'FaceColor','interp','EdgeColor',edgecolor);
else
	error('plot_importancefactors error message: 3d meshes not supported yet');
end

%apply options
applyoptions(md,[],options);

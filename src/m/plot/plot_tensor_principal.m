function plot_tensor_principal(md,options,width,i,tensor,type,plot_options)
%PLOT_TENSOR_PRINCIPAL - plot principal values
%
%   Usage:
%      plot_tensor_principal(md,options,width,i,tensor,type,plot_options);
%
%   See also: PLOTMODEL, PLOT_UNIT, PLOT_MANAGER

%Compute the indexes of the components plots
upperplots=fix((i-1)/width);
if upperplots==0, leftplots=i-1; else leftplots=i-width*upperplots-1; end
if (dimension(md.mesh)==2)%3 components -> 3 indexes
	index1=4*width*upperplots+2*leftplots+1;
	index2=index1+1;
	index3=index1+width*2;
	index4=index3+1;
	newwidth=2*width;
elseif dimension(md.mesh)==3%6 components -> 6 indexes
	index1=3*3*width*upperplots+3*leftplots+1;
	index2=index1+1;
	index3=index1+2;
	index4=index1+width*3;
	index5=index4+1;
	index6=index4+2;
	newwidth=3*width;
end

%plot principal axis
type1=[type 'axis1'];
plot_tensor_principalaxis(md,options,newwidth,index1,tensor,type1,plot_options);
type2=[type 'axis2'];
plot_tensor_principalaxis(md,options,newwidth,index2,tensor,type2,plot_options);
if  dimension(md.mesh)==3
	type3=[type 'axis3'];
	plot_tensor_principalaxis(md,options,newwidth,index3,tensor,type3,plot_options);
end

%now plot principal values
[x y z elements is2d isplanet]=processmesh(md,[],options);
[tensor.principalvalue1 datatype]=processdata(md,tensor.principalvalue1,options);
[tensor.principalvalue2 datatype]=processdata(md,tensor.principalvalue2,options);
if  dimension(md.mesh)==3
	[tensor.principalvalue3 datatype]=processdata(md,tensor.principalvalue3,options);
end

if dimension(md.mesh)==2,
	subplot(2*width,2*width,index3)
	plot_unit(x,y,z,elements,tensor.principalvalue1,is2d,isplanet,datatype,options)
	Apply_options_tensor(md,options,type,'principal value 1')
	subplot(2*width,2*width,index4)
	plot_unit(x,y,z,elements,tensor.principalvalue2,is2d,isplanet,datatype,options)
	Apply_options_tensor(md,options,type,'principal value 2')
else
	subplot(3*width,3*width,index4)
	plot_unit(x,y,z,elements,tensor.principalvalue1,is2d,isplanet,datatype,options)
	Apply_options_tensor(md,options,type,'principal value 1')
	subplot(3*width,3*width,index5)
	plot_unit(x,y,z,elements,tensor.principalvalue2,is2d,isplanet,datatype,options)
	Apply_options_tensor(md,options,type,'principal value 2')
	subplot(3*width,3*width,index6)
	plot_unit(x,y,z,elements,tensor.principalvalue3,is2d,isplanet,datatype,options)
	Apply_options_tensor(md,options,type,'principal value 3')
end
end

function Apply_options_tensor(md,options,type,component)
%apply options
if ismember('_',type) %user plotet stress_tensor
	strings=strsplit_strict(type,'_');
	string=strings{1};
else %default plot: user requested stress
	string=type;
end
options=changefieldvalue(options,'title',[upper(string(1)) string(2:end) ' ' component]);
options=changefieldvalue(options,'colorbar',2);
applyoptions(md,[],options);
end

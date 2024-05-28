function plot_tensor_components(md,options,width,i,tensor,type,plot_options)
%PLOT_TENSOR_COMPONENT - plot component of a tensor
%
%   Usage:
%      plot_tensor_components(md,options,width,i,tensor,type,plot_option);
%
%   See also: PLOTMODEL, PLOT_UNIT, PLOT_MANAGER

%Compute the indexes of the components plots
upperplots=fix((i-1)/width);
if upperplots==0, leftplots=i-1; else leftplots=i-width*upperplots-1; end
if dimension(md.mesh)==2 %3 components -> 3 indexes
	index1=4*width*upperplots+2*leftplots+1;
	index2=index1+1;
	index3=index1+width*2;
elseif dimension(md.mesh)==3%6 components -> 6 indexes
	index1=3*3*width*upperplots+3*leftplots+1;
	index2=index1+1;
	index3=index1+2;
	index4=index1+width*3;
	index5=index4+1;
	index6=index4+2;
end

%process data and model
[x y z elements is2d isplanet]=processmesh(md,[],options);
[tensor.xx datatype]=processdata(md,tensor.xx,options);
[tensor.yy datatype]=processdata(md,tensor.yy,options);
[tensor.xy datatype]=processdata(md,tensor.xy,options);
if  dimension(md.mesh)==3
	[tensor.xz datatype]=processdata(md,tensor.xz,options);
	[tensor.yz datatype]=processdata(md,tensor.yz,options);
	[tensor.zz datatype]=processdata(md,tensor.zz,options);
end

if dimension(md.mesh)==2,
	subplot(2*width,2*width,index1),
	plot_unit(x,y,z,elements,tensor.xx,is2d,isplanet,datatype,options)
	Apply_options_tensor(md,options,type,'xx')
	subplot(2*width,2*width,index2),
	plot_unit(x,y,z,elements,tensor.yy,is2d,isplanet,datatype,options)
	Apply_options_tensor(md,options,type,'yy')
	subplot(2*width,2*width,index3),
	plot_unit(x,y,z,elements,tensor.xy,is2d,isplanet,datatype,options)
	Apply_options_tensor(md,options,type,'xy')
else
	subplot(3*width,3*width,index1),
	plot_unit(x,y,z,elements,tensor.xx,is2d,isplanet,datatype,options)
	Apply_options_tensor(md,options,type,'xx')
	subplot(3*width,3*width,index2),
	plot_unit(x,y,z,elements,tensor.yy,is2d,isplanet,datatype,options)
	Apply_options_tensor(md,options,type,'yy')
	subplot(3*width,3*width,index3),
	plot_unit(x,y,z,elements,tensor.zz,is2d,isplanet,datatype,options)
	Apply_options_tensor(md,options,type,'zz')
	subplot(3*width,3*width,index4),
	plot_unit(x,y,z,elements,tensor.xy,is2d,isplanet,datatype,options)
	Apply_options_tensor(md,options,type,'xy')
	subplot(3*width,3*width,index5),
	plot_unit(x,y,z,elements,tensor.xz,is2d,isplanet,datatype,options)
	Apply_options_tensor(md,options,type,'xz')
	subplot(3*width,3*width,index6),
	plot_unit(x,y,z,elements,tensor.yz,is2d,isplanet,datatype,options)
	Apply_options_tensor(md,options,type,'yz')
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
	applyoptions(md,[],options);
end

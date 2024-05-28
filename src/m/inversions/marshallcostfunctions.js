function  marshallcostfunctions(cost_functions){
	for(var i=0;i<cost_functions.length;i++){
		if(cost_functions[i]==101) data[i]='SurfaceAbsVelMisfit';
		if(cost_functions[i]==102) data[i]='SurfaceRelVelMisfit';
		if(cost_functions[i]==103) data[i]='SurfaceLogVelMisfit';
		if(cost_functions[i]==104) data[i]='SurfaceLogVxVyMisfit';
		if(cost_functions[i]==105) data[i]='SurfaceAverageVelMisfit';
		if(cost_functions[i]==201) data[i]='ThicknessAbsMisfit';
		if(cost_functions[i]==501) data[i]='DragCoefficientAbsGradient';
		if(cost_functions[i]==502) data[i]='RheologyBbarAbsGradient';
		if(cost_functions[i]==503) data[i]='ThicknessAbsGradient';
		if(cost_functions[i]==504) data[i]='ThicknessAlongGradient';
		if(cost_functions[i]==505) data[i]='ThicknessAcrossGradient';
		if(cost_functions[i]==506) data[i]='BalancethicknessMisfit';
		if(cost_functions[i]==507) data[i]='RheologyBAbsGradient';
		if(cost_functions[i]==601) data[i]='SurfaceAbsMisfit';
	}
	return data;
}

function ismodelselfconsistent(md){
//ISMODELSELFCONSISTENT - check that model forms a closed form solvable problem.
//
//   Usage:
//      ismodelselfconsistent(md),

	//initialize consistency as true
	md.priv.isconsistent=true;

	//Get solution and associated analyses
	solution=md.priv.solution;
	if(typeof solution !== 'string')throw Error('ismodelselfconsistent: did not provide correct solution type in the private class!');
	
	var analyses = AnalysisConfiguration(solution);

	//Go through a model field, check that it is a class, and call checkconsistency
	for(field in md){

		//Some properties do not need to be checked
		if (field == 'results' | field == 'debug' | field == 'radaroverlay'){
			continue;
		}

		//Check that current field is a class
		if(typeof md[field] == 'function'){
			continue;
		}

		//Check consistency of the class
		md[field].checkconsistency(md,solution,analyses);
	}

	//error message if mode is not consistent
	if (md.priv.isconsistent==false){
		throw Error('Model not consistent, see messages above');
	}
}

function AnalysisConfiguration(solutiontype){ // {{{
	//ANALYSISCONFIGURATION - return type of analyses, number of analyses 
	//
	//   Usage:
	//      [analyses]=AnalysisConfiguration(solutiontype);

	var analyses=[];
		
	if(solutiontype === 'StressbalanceSolution'){
		analyses=['StressbalanceAnalysis','StressbalanceVerticalAnalysis','StressbalanceSIAAnalysis','L2ProjectionBaseAnalysis'];
		
	}else if(solutiontype ==='SteadystateSolution'){
		analyses=['StressbalanceAnalysis','StressbalanceVerticalAnalysis','StressbalanceSIAAnalysis','L2ProjectionBaseAnalysis','ThermalAnalysis','MeltingAnalysis'];
		
	}else if(solutiontype ==='ThermalSolution'){
		analyses=['EnthalpyAnalysis','ThermalAnalysis','MeltingAnalysis'];
		
	}else if(solutiontype ==='MasstransportSolution'){
		analyses=['MasstransportAnalysis'];
		
	}else if(solutiontype ==='BalancethicknessSolution'){
		analyses=['BalancethicknessAnalysis'];
		
	}else if(solutiontype ==='Balancethickness2Solution'){
		analyses=['Balancethickness2Analysis'];
		
	}else if(solutiontype ==='BalancethicknessSoftSolution'){
		analyses=['BalancethicknessAnalysis'];
		
	}else if(solutiontype ==='BalancevelocitySolution'){
		analyses=['BalancevelocityAnalysis'];
		
	}else if(solutiontype ==='SurfaceSlopeSolution'){
		analyses=['L2ProjectionBaseAnalysis'];
		
	}else if(solutiontype ==='BedSlopeSolution'){
		analyses=['L2ProjectionBaseAnalysis'];
		
	}else if(solutiontype ==='GiaSolution'){
		analyses=['GiaIvinsAnalysis'];
		
	}else if(solutiontype ==='TransientSolution'){
		analyses=['StressbalanceAnalysis','StressbalanceVerticalAnalysis','StressbalanceSIAAnalysis','L2ProjectionBaseAnalysis','ThermalAnalysis','MeltingAnalysis','EnthalpyAnalysis','MasstransportAnalysis','HydrologyShaktiAnalysis','HydrologyGladsAnalysis'];
		
	}else if(solutiontype ==='SealevelriseSolution'){
		analyses=['SealevelriseAnalysis'];
		
	}else if(solutiontype ==='HydrologySolution'){
		analyses=['L2ProjectionBaseAnalysis','HydrologyShreveAnalysis','HydrologyDCInefficientAnalysis','HydrologyDCEfficientAnalysis'];
		
	}else if(solutiontype ==='DamageEvolutionSolution'){
		analyses=['DamageEvolutionAnalysis'];
	}else{
		throw Error(sprintf("%s%s%s\n",' solution type: ',solutiontype,' not supported yet!'));
	}
	return analyses;
} // }}}

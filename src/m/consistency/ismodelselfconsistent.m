function ismodelselfconsistent(md)
%ISMODELSELFCONSISTENT - check that model forms a closed form solvable problem.
%
%   Usage:
%      ismodelselfconsistent(md);

%initialize consistency as true
md.private.isconsistent=true;

%Get solution and associated analyses
solution=md.private.solution;
[analyses]=AnalysisConfiguration(solution);

%Go through a model field, check that it is a class, and call checkconsistency
fields=properties('model');
for i=1:length(fields)
	field=fields{i};

	%Some properties do not need to be checked
	if ismember(field,{'results' 'debug' 'radaroverlay'}),
		continue;
	end

	%Check that current field is an object
	if ~isobject(md.(field))
		md=checkmessage(md,['field ''' char(field) ''' is not an object']);
		continue;
	end

	%Check consistency of the object
	md=checkconsistency(md.(field),md,solution,analyses);
end

%error message if mode is not consistent
if md.private.isconsistent==false
	error('Model not consistent, see messages above');
end
end

function [analyses]=AnalysisConfiguration(solutiontype), % {{{
%ANALYSISCONFIGURATION - return type of analyses, number of analyses
%
%   Usage:
%      [analyses]=AnalysisConfiguration(solutiontype);

	if strcmp(solutiontype,'StressbalanceSolution')
		analyses={'StressbalanceAnalysis','StressbalanceVerticalAnalysis','StressbalanceSIAAnalysis','L2ProjectionBaseAnalysis'};
	elseif strcmp(solutiontype,'SteadystateSolution')
		analyses={'StressbalanceAnalysis','StressbalanceVerticalAnalysis','StressbalanceSIAAnalysis','L2ProjectionBaseAnalysis','ThermalAnalysis','MeltingAnalysis','EnthalpyAnalysis','AgeAnalysis'};
	elseif strcmp(solutiontype,'ThermalSolution')
		analyses={'EnthalpyAnalysis','ThermalAnalysis','MeltingAnalysis'};
	elseif strcmp(solutiontype,'MasstransportSolution')
		analyses={'MasstransportAnalysis'};
	elseif strcmp(solutiontype,'OceantransportSolution')
		analyses={'OceantransportAnalysis'};
	elseif strcmp(solutiontype,'BalancethicknessSolution')
		analyses={'BalancethicknessAnalysis'};
	elseif strcmp(solutiontype,'Balancethickness2Solution')
		analyses={'Balancethickness2Analysis'};
	elseif strcmp(solutiontype,'BalancethicknessSoftSolution')
		analyses={'BalancethicknessAnalysis'};
	elseif strcmp(solutiontype,'BalancevelocitySolution')
		analyses={'BalancevelocityAnalysis'};
	elseif strcmp(solutiontype,'SurfaceSlopeSolution')
		analyses={'L2ProjectionBaseAnalysis'};
	elseif strcmp(solutiontype,'BedSlopeSolution')
		analyses={'L2ProjectionBaseAnalysis'};
	elseif strcmp(solutiontype,'GiaSolution')
		analyses={'GiaIvinsAnalysis'};
	elseif strcmp(solutiontype,'LoveSolution')
		analyses={'LoveAnalysis'};
	elseif strcmp(solutiontype,'EsaSolution')
		analyses={'EsaAnalysis'};
	elseif strcmp(solutiontype,'TransientSolution')
		analyses={'StressbalanceAnalysis','StressbalanceVerticalAnalysis','StressbalanceSIAAnalysis','L2ProjectionBaseAnalysis','ThermalAnalysis','MeltingAnalysis','EnthalpyAnalysis','MasstransportAnalysis','OceantransportAnalysis','HydrologyShaktiAnalysis','HydrologyGladsAnalysis','HydrologyShreveAnalysis','HydrologyTwsAnalysis','HydrologyDCInefficientAnalysis','HydrologyDCEfficientAnalysis','SealevelchangeAnalysis','AgeAnalysis','HydrologyArmapwAnalysis','AgeAnalysis','DebrisAnalysis'};
	elseif strcmp(solutiontype,'SealevelchangeSolution')
		analyses={'SealevelchangeAnalysis'};
	elseif strcmp(solutiontype,'HydrologySolution')
		analyses={'L2ProjectionBaseAnalysis','HydrologyShreveAnalysis','HydrologyDCInefficientAnalysis','HydrologyDCEfficientAnalysis','HydrologyGladsAnalysis','HydrologyShaktiAnalysis','HydrologyTwsAnalysis','HydrologyArmapwAnalysis'};
	elseif strcmp(solutiontype,'DamageEvolutionSolution')
		analyses={'DamageEvolutionAnalysis'};
    elseif strcmp(solutiontype,'SamplingSolution')
		analyses={'SamplingAnalysis'};    
	else
		error(' solution type: %s' , solutiontype, ' not supported yet!');
	end
end % }}}

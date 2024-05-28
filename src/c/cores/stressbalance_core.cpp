/*!\file: stressbalance_core.cpp
 * \brief: core of the stressbalance solution
 */

#include "./cores.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../analyses/analyses.h"
#include "../shared/shared.h"
#include "../modules/modules.h"
#include "../solutionsequences/solutionsequences.h"

void stressbalance_core(FemModel* femmodel){

	/*Start profiler*/
	femmodel->profiler->Start(STRESSBALANCECORE);

	/*parameters: */
	bool       dakota_analysis,control_analysis;
	int        domaintype;
	bool       isSIA,isSSA,isL1L2,isMOLHO,isHO,isFS,isNitsche;
	bool       save_results;
	int        solution_type;
	int        numoutputs        = 0;
	char     **requested_outputs = NULL;
	Analysis  *analysis          = NULL;

	/* recover parameters:*/
	femmodel->parameters->FindParam(&domaintype,DomainTypeEnum);
	femmodel->parameters->FindParam(&isSIA,FlowequationIsSIAEnum);
	femmodel->parameters->FindParam(&isSSA,FlowequationIsSSAEnum);
	femmodel->parameters->FindParam(&isL1L2,FlowequationIsL1L2Enum);
	femmodel->parameters->FindParam(&isMOLHO,FlowequationIsMOLHOEnum);
	femmodel->parameters->FindParam(&isHO,FlowequationIsHOEnum);
	femmodel->parameters->FindParam(&isFS,FlowequationIsFSEnum);
	femmodel->parameters->FindParam(&isNitsche,FlowequationIsNitscheEnum);
	femmodel->parameters->FindParam(&save_results,SaveResultsEnum);
	femmodel->parameters->FindParam(&solution_type,SolutionTypeEnum);
	femmodel->parameters->FindParam(&numoutputs,StressbalanceNumRequestedOutputsEnum);
	femmodel->parameters->FindParam(&control_analysis,InversionIscontrolEnum);
	if(numoutputs) femmodel->parameters->FindParam(&requested_outputs,&numoutputs,StressbalanceRequestedOutputsEnum);

	if(VerboseSolution()) _printf0_("   computing new velocity\n");
	/*Compute slopes if necessary */
	if(isSIA || (isFS && domaintype==Domain2DverticalEnum)) surfaceslope_core(femmodel);
	if(isFS){
		/*We need bed slopoes for the non-penetration condition*/
		bedslope_core(femmodel);
		femmodel->SetCurrentConfiguration(StressbalanceAnalysisEnum);
		ResetFSBasalBoundaryConditionx(femmodel->elements,femmodel->nodes,femmodel->vertices,femmodel->loads,femmodel->materials,femmodel->parameters);

		/*We need basal melt rates for the shelf dampening*/
		bmb_core(femmodel);
	}

	/*Compute SIA velocities*/
	if(isSIA){

		/*Take the last velocity into account so that the velocity on the SSA domain is not zero*/
		if(isSSA || isL1L2 || isHO ){
			femmodel->SetCurrentConfiguration(StressbalanceSIAAnalysisEnum);
			ResetBoundaryConditions(femmodel,StressbalanceSIAAnalysisEnum);
		}

		analysis = new StressbalanceSIAAnalysis();
		analysis->Core(femmodel);
		delete analysis;

		/*Reset velocities for other ice flow models*/
		if(isSSA || isL1L2 || isHO){
			femmodel->SetCurrentConfiguration(StressbalanceAnalysisEnum);
			ResetBoundaryConditions(femmodel,StressbalanceAnalysisEnum);
		}
	}

	/*Compute stressbalance for SSA L1L2 HO and FS*/
	if(isSSA || isL1L2 || isMOLHO || isHO || isFS){
		analysis = new StressbalanceAnalysis();
		analysis->Core(femmodel);
		delete analysis;
	}

	/*Compute vertical velocities*/
	if (domaintype==Domain3DEnum && (isSIA || isSSA || isL1L2 || isMOLHO || isHO)){

		/*We need basal melt rates for vertical velocity*/
		bmb_core(femmodel);

		analysis = new StressbalanceVerticalAnalysis();
		analysis->Core(femmodel);
		delete analysis;
	}

	if(save_results){
		femmodel->RequestedOutputsx(&femmodel->results,requested_outputs,numoutputs);
	}

	if(solution_type==StressbalanceSolutionEnum && !control_analysis)femmodel->RequestedDependentsx();

	/*Free resources:*/
	if(numoutputs){for(int i=0;i<numoutputs;i++){xDelete<char>(requested_outputs[i]);} xDelete<char*>(requested_outputs);}

	/*End profiler*/
	femmodel->profiler->Stop(STRESSBALANCECORE);
}

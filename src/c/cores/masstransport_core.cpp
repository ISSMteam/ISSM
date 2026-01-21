/*!\file: masstransport_core.cpp
 * \brief: core of the masstransport solution 
 */ 

#include "./cores.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"
#include "../solutionsequences/solutionsequences.h"
#include "../classes/Inputs/TransientInput.h"
void SolidEarthIceUpdates(FemModel* femmodel);
void masstransport_core(FemModel* femmodel){ /*{{{*/

	/*Start profiler*/
	femmodel->profiler->Start(MASSTRANSPORTCORE);

	/*parameters: */
	int    numoutputs,domaintype;
	bool   save_results;
	bool   isFS,isfreesurface,dakota_analysis;
	int    solution_type,stabilization;
	char** requested_outputs = NULL;

	/*activate configuration*/
	femmodel->SetCurrentConfiguration(MasstransportAnalysisEnum);

	/*recover parameters: */
	femmodel->parameters->FindParam(&save_results,SaveResultsEnum);
	femmodel->parameters->FindParam(&isFS,FlowequationIsFSEnum);
	femmodel->parameters->FindParam(&isfreesurface,MasstransportIsfreesurfaceEnum);
	femmodel->parameters->FindParam(&solution_type,SolutionTypeEnum);
	femmodel->parameters->FindParam(&domaintype,DomainTypeEnum);
	femmodel->parameters->FindParam(&numoutputs,MasstransportNumRequestedOutputsEnum);
	femmodel->parameters->FindParam(&stabilization,MasstransportStabilizationEnum);
	if(numoutputs) femmodel->parameters->FindParam(&requested_outputs,&numoutputs,MasstransportRequestedOutputsEnum);

	if(VerboseSolution()) _printf0_("   computing mass transport\n");

	/*Transport mass or free surface*/
	if(isFS && isfreesurface){
		if(VerboseSolution()) _printf0_("   call free surface computational core\n");
		femmodel->SetCurrentConfiguration(FreeSurfaceBaseAnalysisEnum);
		solutionsequence_linear(femmodel);
		femmodel->parameters->SetParam(BaseEnum,InputToExtrudeEnum);
		extrudefrombase_core(femmodel);

		femmodel->SetCurrentConfiguration(FreeSurfaceTopAnalysisEnum);
		solutionsequence_linear(femmodel);
		femmodel->parameters->SetParam(SurfaceEnum,InputToExtrudeEnum);
		extrudefromtop_core(femmodel);
		femmodel->parameters->SetParam(ThicknessEnum,InputToExtrudeEnum);
		extrudefromtop_core(femmodel);
		femmodel->parameters->SetParam(BaseEnum,InputToExtrudeEnum);
		extrudefromtop_core(femmodel);
	}
	else{
		femmodel->parameters->SetParam(VxEnum,InputToDepthaverageInEnum);
		femmodel->parameters->SetParam(VxAverageEnum,InputToDepthaverageOutEnum);
		depthaverage_core(femmodel);
		if(domaintype!=Domain2DverticalEnum){
			femmodel->parameters->SetParam(VyEnum,InputToDepthaverageInEnum);
			femmodel->parameters->SetParam(VyAverageEnum,InputToDepthaverageOutEnum);
			depthaverage_core(femmodel);
		}
		femmodel->SetCurrentConfiguration(MasstransportAnalysisEnum);
		InputDuplicatex(femmodel,ThicknessEnum,ThicknessOldEnum);
		InputDuplicatex(femmodel,BaseEnum,BaseOldEnum);
		InputDuplicatex(femmodel,SurfaceEnum,SurfaceOldEnum);
		if(stabilization==4){
			solutionsequence_fct(femmodel);
		}
		else{
			solutionsequence_linear(femmodel);
			// ThicknessAverage: method not totally tested
			//if(stabilization==3){
			//	if(VerboseSolution()) _printf0_("   call thickness average core\n");
			//	femmodel->ThicknessAverage();
			//}
		}
		SolidEarthIceUpdates(femmodel);

		femmodel->parameters->SetParam(ThicknessEnum,InputToExtrudeEnum);
		extrudefrombase_core(femmodel);
		femmodel->parameters->SetParam(BaseEnum,InputToExtrudeEnum);
		extrudefrombase_core(femmodel);
		femmodel->parameters->SetParam(SurfaceEnum,InputToExtrudeEnum);
		extrudefrombase_core(femmodel);

	}

	if(save_results){
		femmodel->RequestedOutputsx(&femmodel->results,requested_outputs,numoutputs);
	}

	if(solution_type==MasstransportSolutionEnum)femmodel->RequestedDependentsx();

	/*Free resources:*/
	if(numoutputs){for(int i=0;i<numoutputs;i++){xDelete<char>(requested_outputs[i]);} xDelete<char*>(requested_outputs);}

	/*profiler*/
	femmodel->profiler->Stop(MASSTRANSPORTCORE);
} /*}}}*/
void SolidEarthIceUpdates(FemModel* femmodel){ /*{{{*/

	int isgrd;
	int grdmodel;
	IssmDouble time;
	int frequency,count;

	/*retrieve parameters:*/
	femmodel->parameters->FindParam(&isgrd,SolidearthSettingsGRDEnum);
	femmodel->parameters->FindParam(&grdmodel,GrdModelEnum);
	femmodel->parameters->FindParam(&time,TimeEnum);
	femmodel->parameters->FindParam(&frequency,SolidearthSettingsRunFrequencyEnum);
	femmodel->parameters->FindParam(&count,SealevelchangeRunCountEnum);

	/*early return?:*/
	if(!isgrd)return;

	/*From old and new thickness, create delta thickness  and accumulate:*/
	femmodel->inputs->ZAXPY(-1, ThicknessOldEnum,ThicknessEnum,DeltaIceThicknessEnum);
	femmodel->inputs->AXPY(+1, DeltaIceThicknessEnum,AccumulatedDeltaIceThicknessEnum);

	/*for Ivins deformation model, keep history of ice thickness changes inside TransientAccumulatedDeltaIceThicknessEnum:*/
	if(grdmodel==IvinsEnum){

		TransientInput* transientinput = femmodel->inputs->GetTransientInput(TransientAccumulatedDeltaIceThicknessEnum);

		for(Object* & object : femmodel->elements->objects){
			int *vertexlids = NULL;
			int *vertexsids= NULL;
			Element*   element=xDynamicCast<Element*>(object);
			const int numvertices = element->GetNumberOfVertices();
			IssmDouble* cumdeltathickness=NULL;

			/*Get values and lid list and recover vertices ids needed to initialize inputs*/
			vertexlids      = xNew<int>(numvertices);
			vertexsids      = xNew<int>(numvertices);
			cumdeltathickness=xNew<IssmDouble>(numvertices);
			element->GetVerticesLidList(&vertexlids[0]);
			element->GetVerticesSidList(&vertexsids[0]);
			element->GetInputListOnVertices(&cumdeltathickness[0],AccumulatedDeltaIceThicknessEnum);

			/*Add the current time cumdeltathickness to the existing time series: */
			switch(element->ObjectEnum()){
				case TriaEnum:  transientinput->AddTriaTimeInput( time,numvertices,vertexlids,cumdeltathickness,P1Enum); break;
				default: _error_("Not implemented yet");
			}
			xDelete<int>(vertexlids);
			xDelete<int>(vertexsids);
			xDelete<IssmDouble>(cumdeltathickness);
		}
	}	

	/*compute total ice thickness change between two sea-level solver time steps, ie. every frequency*dt:*/
	if(count==frequency){
		femmodel->inputs->ZAXPY(-1, OldAccumulatedDeltaIceThicknessEnum,AccumulatedDeltaIceThicknessEnum,DeltaIceThicknessEnum);
		femmodel->inputs->DuplicateInput(AccumulatedDeltaIceThicknessEnum,OldAccumulatedDeltaIceThicknessEnum);
	}
	return;
}/*}}}*/

/*!\file: hydrology_core.cpp
 * \brief: core of the hydrology solution
 */

#include "./cores.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"
#include "../solutionsequences/solutionsequences.h"
void SolidEarthWaterUpdates(FemModel* femmodel);

void hydrology_core(FemModel* femmodel){ /*{{{*/

	/*Start profiler*/
	femmodel->profiler->Start(HYDROLOGYCORE);

	/*intermediary*/
	int          hydrology_model;
	int          solution_type;
	int          numoutputs        = 0;
	bool         save_results;
	bool         modify_loads      = true;
	char       **requested_outputs = NULL;
	IssmDouble   ThawedNodes;

	/*first recover parameters common to all solutions*/
	femmodel->parameters->FindParam(&save_results,SaveResultsEnum);
	femmodel->parameters->FindParam(&hydrology_model,HydrologyModelEnum);
	femmodel->parameters->FindParam(&solution_type,SolutionTypeEnum);
	femmodel->parameters->FindParam(&numoutputs,HydrologyNumRequestedOutputsEnum);
	if(numoutputs) femmodel->parameters->FindParam(&requested_outputs,&numoutputs,HydrologyRequestedOutputsEnum);

	/*Using the Shreve based Model*/
	if (hydrology_model==HydrologyshreveEnum){
		if(VerboseSolution()) _printf0_("   computing water heads\n");
		/*first compute slopes: */
		surfaceslope_core(femmodel);
		bedslope_core(femmodel);
		/*and then go to water column*/
		if(VerboseSolution()) _printf0_("   computing water column\n");
		femmodel->SetCurrentConfiguration(HydrologyShreveAnalysisEnum);
		solutionsequence_nonlinear(femmodel,modify_loads);
		/*transfer water column thickness to old water column thickness: */
		InputDuplicatex(femmodel,WatercolumnEnum,WaterColumnOldEnum);
		/*solid earth considerations:*/
		SolidEarthWaterUpdates(femmodel);

	}
	/*Using the Tws based Model*/
	if (hydrology_model==HydrologyTwsEnum){
		if(VerboseSolution()) _printf0_("   computing water column\n");

		femmodel->SetCurrentConfiguration(HydrologyTwsAnalysisEnum);

		/*save current tws  before updating:*/
		InputDuplicatex(femmodel,WatercolumnEnum,WaterColumnOldEnum);

		/*grab tws from the hydrology.spcwatercolumn field input and update
		 * the solution with it:*/
		Vector<IssmDouble>*  ug  = NULL;
		GetVectorFromInputsx(&ug,femmodel,HydrologyTwsSpcEnum,VertexPIdEnum);
		InputUpdateFromSolutionx(femmodel,ug);

		/*solid earth considerations:*/
		SolidEarthWaterUpdates(femmodel);
		delete ug;
	}

	/*Using the double continuum model*/
	else if (hydrology_model==HydrologydcEnum){
		/*intermediary: */
		bool       isefficientlayer;
		/*recover parameters: */
		femmodel->parameters->FindParam(&isefficientlayer,HydrologydcIsefficientlayerEnum);

		/*first we exclude frozen nodes of the solved nodes*/
		femmodel->SetCurrentConfiguration(HydrologyDCInefficientAnalysisEnum);
		femmodel->HydrologyIDSupdateDomainx(&ThawedNodes);

		if(ThawedNodes>0){
			/*check if we need sub steps*/
			int  dtslices;
			bool sliceadapt;
         bool conv_fail=false;
         femmodel->parameters->FindParam(&dtslices,HydrologyStepsPerStepEnum);
         femmodel->parameters->FindParam(&sliceadapt,HydrologyStepAdaptEnum);

			if(dtslices>1 || sliceadapt){
				int        step, substep, numaveragedinput, hydro_averaging, remainingslices;
				IssmDouble global_time, subtime, yts, remainingtime, averagetime;
				IssmDouble dt, subdt;

            femmodel->parameters->FindParam(&global_time,TimeEnum);
				femmodel->parameters->FindParam(&step,StepEnum);
            femmodel->parameters->FindParam(&dt,TimesteppingTimeStepEnum);
            femmodel->parameters->FindParam(&yts,ConstantsYtsEnum);
				femmodel->parameters->FindParam(&hydro_averaging,HydrologyAveragingEnum);

				subtime=global_time-dt; //getting the time back to the start of the timestep
				subdt=dt/dtslices; //computing hydro dt from dt and a divider
				substep=0;
				femmodel->parameters->SetParam(subdt,TimesteppingTimeStepEnum);
				femmodel->parameters->SetParam(subtime,TimeEnum);

				/*intermiedaries to deal with averaging*/
				static const int substeplist[4] = {EffectivePressureSubstepEnum,SedimentHeadSubstepEnum,EplHeadSubstepEnum,HydrologydcEplThicknessSubstepEnum};
				static const int transientlist[4] = {EffectivePressureTransientEnum,SedimentHeadTransientEnum,EplHeadTransientEnum,HydrologydcEplThicknessTransientEnum};
				static const int averagelist[4] = {EffectivePressureEnum,SedimentHeadEnum,EplHeadEnum,HydrologydcEplThicknessEnum};
				std::vector<int> substepinput;
				std::vector<int> transientinput;
				std::vector<int> averagedinput;

				if (isefficientlayer){
					/*define which variable needs to be averaged on the sub-timestep and initialize as needed*/
					numaveragedinput = 4;
					substepinput.assign(substeplist,substeplist+4);
					transientinput.assign(transientlist,transientlist+4);
					averagedinput.assign(averagelist,averagelist+4);
				}
				else{
					numaveragedinput = 2;
					substepinput.assign(substeplist,substeplist+2);
					transientinput.assign(transientlist,transientlist+2);
					averagedinput.assign(averagelist,averagelist+2);
				}
				femmodel->InitTransientInputx(&transientinput[0],numaveragedinput);
				averagetime=0;
				while(substep<dtslices){ //loop on hydro dts
					substep+=1;
					subtime+=subdt;
					averagetime+=subtime*subdt;
					/*Setting substep time as global time*/
					femmodel->parameters->SetParam(subtime,TimeEnum);
					if(VerboseSolution()) _printf0_("sub iteration " << substep << "/" << dtslices << "  time [yr]: " << setprecision(4) << subtime/yts << " (time step: " << subdt/yts << ")\n");
					if(VerboseSolution()) _printf0_("   computing water heads\n");
					/*save preceding timestep*/
					InputDuplicatex(femmodel,SedimentHeadSubstepEnum,SedimentHeadOldEnum);
					if (isefficientlayer){
						InputDuplicatex(femmodel,EplHeadSubstepEnum,EplHeadOldEnum);
						InputDuplicatex(femmodel,HydrologydcEplThicknessSubstepEnum,HydrologydcEplThicknessOldEnum);
					}
					/*Proceed now to heads computations*/
					solutionsequence_hydro_nonlinear(femmodel, &conv_fail);
					if(conv_fail){
                  /*convergence failed, we want to go back to the begining of the main step and increase the number of subslices*/
                  /*First we get teh time and step counter back to the begining of the step that did not converge*/
                  averagetime-=subtime*subdt;
                  subtime-=subdt;
                  substep-=1;
                  /*compute the number of slice that are remaining and the time left in the timestep*/
                  remainingslices=dtslices-substep;
                  remainingtime=global_time-subtime;
                  /*We double the number of remaining slices and compute their duration*/
                  dtslices=dtslices-remainingslices+(2*remainingslices);
                  subdt=remainingtime/(2*remainingslices);
                  if(VerboseSolution())_printf0_("convergence failed for sub-step "<< substep <<" total number of slice is now "<< dtslices <<" for step "<<step<<"\n");
                  if(VerboseSolution())_printf0_("next slice duration is "<< subdt/yts <<" years\n");
                  conv_fail = false;  //re-initialize the control keyword
                  if (dtslices>500){
                     _error_("   We reached (" << dtslices << ") which exceeds the hard limit of 500");
                  }

                  femmodel->parameters->SetParam(subdt,TimesteppingTimeStepEnum);
                  femmodel->parameters->SetParam(subtime,TimeEnum);

               }
               else{
						/*If we have a sub-timestep we store the substep inputs in a transient input here*/
						femmodel->StackTransientInputonBasex(&substepinput[0],&transientinput[0],subtime,numaveragedinput);
					}
				}
				/*averaging the stack*/
				femmodel->AverageTransientInputonBasex(&transientinput[0],&averagedinput[0],global_time-dt,subtime,numaveragedinput,hydro_averaging);
				/*And reseting to global time*/
				femmodel->parameters->SetParam(dt,TimesteppingTimeStepEnum);
				femmodel->parameters->SetParam(global_time,TimeEnum);
				if(save_results){
               femmodel->results->AddResult(new GenericExternalResult<int>(femmodel->results->Size()+1,HydrologySubstepsEnum,dtslices,step,global_time));
               femmodel->results->AddResult(new GenericExternalResult<IssmDouble>(femmodel->results->Size()+1,HydrologySubTimeEnum,(averagetime/dt)/yts,step,global_time));
            }
			}
			else{
				InputDuplicatex(femmodel,SedimentHeadSubstepEnum,SedimentHeadOldEnum);
				if (isefficientlayer){
					InputDuplicatex(femmodel,EplHeadSubstepEnum,EplHeadOldEnum);
					InputDuplicatex(femmodel,HydrologydcEplThicknessSubstepEnum,HydrologydcEplThicknessOldEnum);
				}
				/*Proceed now to heads computations*/
				if(VerboseSolution()) _printf0_("   computing water heads\n");
				solutionsequence_hydro_nonlinear(femmodel, &conv_fail);
				/*If no substeps are present we want to duplicate the results for coupling purposes*/
				InputDuplicatex(femmodel,SedimentHeadSubstepEnum,SedimentHeadEnum);
				InputDuplicatex(femmodel,EffectivePressureSubstepEnum,EffectivePressureEnum);
				if (isefficientlayer){
					InputDuplicatex(femmodel,EplHeadSubstepEnum,EplHeadEnum);
					InputDuplicatex(femmodel,HydrologydcEplThicknessSubstepEnum,HydrologydcEplThicknessEnum);
				}
			}
		}
		if(VerboseSolution())_printf0_("   hydroDC done\n");
	}

	/*Using the SHAKTI model*/
	else if (hydrology_model==HydrologyshaktiEnum){

		/*Get second derivatives of gap height*/
		femmodel->SetCurrentConfiguration(L2ProjectionBaseAnalysisEnum);
		femmodel->parameters->SetParam(HydrologyGapHeightXEnum,InputToL2ProjectEnum);
		solutionsequence_linear(femmodel);
		femmodel->parameters->SetParam(HydrologyGapHeightYEnum,InputToL2ProjectEnum);
		solutionsequence_linear(femmodel);
		femmodel->parameters->SetParam(HydrologyGapHeightXXEnum,InputToL2ProjectEnum);
		solutionsequence_linear(femmodel);
		femmodel->parameters->SetParam(HydrologyGapHeightYYEnum,InputToL2ProjectEnum);
		solutionsequence_linear(femmodel);

		/*Update Effective pressure*/
		if(VerboseSolution()) _printf0_("   computing effective pressure\n");
		HydrologyShaktiAnalysis* analysis = new HydrologyShaktiAnalysis();
		analysis->UpdateEffectivePressure(femmodel);

		/*Get new head*/
		femmodel->SetCurrentConfiguration(HydrologyShaktiAnalysisEnum);
		InputDuplicatex(femmodel,HydrologyHeadEnum,HydrologyHeadOldEnum);
		solutionsequence_shakti_nonlinear(femmodel);

		/*Update Gap Height*/
		if(VerboseSolution()) _printf0_("   updating gap height\n");
		analysis->UpdateGapHeight(femmodel);
		delete analysis;
	}

	/*Using the GlaDS model*/
	else if (hydrology_model==HydrologyGlaDSEnum){
		HydrologyGlaDSAnalysis* analysis = new HydrologyGlaDSAnalysis();
		femmodel->SetCurrentConfiguration(HydrologyGlaDSAnalysisEnum);

		/*Set fields as old*/
		InputDuplicatex(femmodel,HydraulicPotentialEnum,HydraulicPotentialOldEnum);
		InputDuplicatex(femmodel,HydrologySheetThicknessEnum,HydrologySheetThicknessOldEnum);
		analysis->SetChannelCrossSectionOld(femmodel);

		/*Solve for new potential*/
		solutionsequence_glads_nonlinear(femmodel);

		if(VerboseSolution()) _printf0_("   updating effective pressure\n");
		analysis->UpdateEffectivePressure(femmodel);
		delete analysis;
	}

	/*Using the PISM hydrology model*/
	else if (hydrology_model==HydrologypismEnum){
		femmodel->SetCurrentConfiguration(HydrologyPismAnalysisEnum);
		if(VerboseSolution()) _printf0_("   updating water column\n");
		HydrologyPismAnalysis* analysis = new HydrologyPismAnalysis();
		InputDuplicatex(femmodel,WatercolumnEnum,WaterColumnOldEnum);
		analysis->UpdateWaterColumn(femmodel);
		delete analysis;
	}

	/*Using the armaPw hydrology model*/
   else if (hydrology_model==HydrologyarmapwEnum){
      femmodel->SetCurrentConfiguration(HydrologyArmapwAnalysisEnum);
      if(VerboseSolution()) _printf0_("   updating subglacial water pressure\n");
      HydrologyArmapwAnalysis* analysis = new HydrologyArmapwAnalysis();
      analysis->UpdateSubglacialWaterPressure(femmodel);
      delete analysis;
   }
	else{
		_error_("Hydrology model "<< EnumToStringx(hydrology_model) <<" not supported yet");
	}
	if(save_results){
		if(hydrology_model==HydrologydcEnum && ThawedNodes==0){
			if(VerboseSolution()) _printf0_("   No thawed node hydro is skiped \n");}
		else{
			femmodel->RequestedOutputsx(&femmodel->results,requested_outputs,numoutputs);
		}
	}
	/*Free resources:*/
	if(numoutputs){
		for(int i=0;i<numoutputs;i++){
			xDelete<char>(requested_outputs[i]);
		}
		xDelete<char*>(requested_outputs);
	}
	/*End profiler*/
	femmodel->profiler->Stop(HYDROLOGYCORE);
} /*}}}*/
void SolidEarthWaterUpdates(FemModel* femmodel){ /*{{{*/

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
	femmodel->inputs->ZAXPY(-1, WaterColumnOldEnum,WatercolumnEnum,DeltaTwsEnum);
	femmodel->inputs->AXPY(+1, DeltaTwsEnum,AccumulatedDeltaTwsEnum);

	/*compute total water column change between two sea-level solver time steps, ie. every frequency*dt:*/
	if(count==frequency){
		femmodel->inputs->ZAXPY(-1, OldAccumulatedDeltaTwsEnum,AccumulatedDeltaTwsEnum,DeltaTwsEnum);
		femmodel->inputs->DuplicateInput(AccumulatedDeltaTwsEnum,OldAccumulatedDeltaTwsEnum);
	}
	return;
}/*}}}*/

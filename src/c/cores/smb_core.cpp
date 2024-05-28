/*!\file: smb_core.cpp
 * \brief: core of the smb solution
 */

#include "./cores.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"
#include "../solutionsequences/solutionsequences.h"

void smb_core(FemModel* femmodel){

	/*Start profiler*/
	femmodel->profiler->Start(SMBCORE);

	/*parameters: */
	Analysis* analysis=NULL;
	int    smb_model;
	int    numoutputs;
	bool   save_results;
	int    solution_type;
	char** requested_outputs = NULL;

	/*activate configuration*/
	femmodel->SetCurrentConfiguration(SmbAnalysisEnum);

	/*recover parameters: */
	femmodel->parameters->FindParam(&smb_model,SmbEnum);
	femmodel->parameters->FindParam(&save_results,SaveResultsEnum);
	femmodel->parameters->FindParam(&solution_type,SolutionTypeEnum);
	femmodel->parameters->FindParam(&numoutputs,SmbNumRequestedOutputsEnum);
	if(numoutputs) femmodel->parameters->FindParam(&requested_outputs,&numoutputs,SmbRequestedOutputsEnum);

	/*sub steping specifics*/
	int dtslices;
	int numaveragedinput;
	femmodel->parameters->FindParam(&dtslices,SmbStepsPerStepEnum);
	/*intermediaries to deal with averaging*/
	static const int substeplist[2] = {SmbMassBalanceSubstepEnum,SmbRunoffSubstepEnum};
	static const int transientlist[2] = {SmbMassBalanceTransientEnum,SmbRunoffTransientEnum};
	static const int averagelist[2] = {SmbMassBalanceEnum,SmbRunoffEnum};
	std::vector<int> substepinput;
	std::vector<int> transientinput;
	std::vector<int> averagedinput;

	/*define which variable needs to be averaged on the sub-timestep and initialize as needed*/
	if(smb_model==SMBgradientscomponentsEnum){
		numaveragedinput = 2;
		substepinput.assign(substeplist,substeplist+2);
		transientinput.assign(transientlist,transientlist+2);
		averagedinput.assign(averagelist,averagelist+2);
	}

	/*if yes compute necessary intermediaries and start looping*/
	if (dtslices>1){
		int        substep,smb_averaging;
		IssmDouble global_time,subtime,yts;
		IssmDouble dt,subdt;

		femmodel->parameters->FindParam(&global_time,TimeEnum);
		femmodel->parameters->FindParam(&dt,TimesteppingTimeStepEnum);
		femmodel->parameters->FindParam(&yts,ConstantsYtsEnum);
		femmodel->parameters->FindParam(&smb_averaging,SmbAveragingEnum);

		subtime=global_time-dt; //getting the time back to the start of the timestep
		subdt=dt/dtslices; //computing substep from dt and a divider
		substep=0;
		femmodel->parameters->SetParam(subdt,TimesteppingTimeStepEnum);

		femmodel->InitTransientInputx(&transientinput[0],numaveragedinput);
		analysis = new SmbAnalysis();
		while(substep<dtslices){ //loop on sub dts
			substep+=1;
			subtime+=subdt;
			femmodel->parameters->SetParam(subtime,TimeEnum);
         if(VerboseSolution()) _printf0_("sub iteration " << substep << "/" << dtslices << "  time [yr]: " << setprecision(4) << subtime/yts << " (time step: " << subdt/yts << ")\n");
         if(VerboseSolution()) _printf0_("   computing smb\n");
			if(VerboseSolution()) _printf0_("   Calling core\n");
			analysis->Core(femmodel);
         /*If we have a sub-timestep we store the substep inputs in a transient input here*/
         femmodel->StackTransientInputx(&substepinput[0],&transientinput[0],subtime,numaveragedinput);
		}
		delete analysis;
      /*averaging the transient input*/
		femmodel->AverageTransientInputx(&transientinput[0],&averagedinput[0],global_time-dt,subtime,numaveragedinput,smb_averaging);
		/*and reset timesteping variables to original*/
		femmodel->parameters->SetParam(global_time,TimeEnum);
		femmodel->parameters->SetParam(dt,TimesteppingTimeStepEnum);
	}
	else{
      if(VerboseSolution()) _printf0_("   computing smb \n");
      analysis = new SmbAnalysis();
		analysis->Core(femmodel);
		/*If no substeps are present we want to duplicate the computed substep enum for coupling purposes*/
		if(smb_model==SMBgradientscomponentsEnum){
			for(int i=0;i<numaveragedinput;i++){
				InputDuplicatex(femmodel,substepinput[i],averagedinput[i]);
			}
		}
		delete analysis;
	}

	if(save_results){
		femmodel->RequestedOutputsx(&femmodel->results,requested_outputs,numoutputs);
	}

	if(solution_type==SmbSolutionEnum)femmodel->RequestedDependentsx();

	/*Free resources:*/
	if(numoutputs){for(int i=0;i<numoutputs;i++){xDelete<char>(requested_outputs[i]);} xDelete<char*>(requested_outputs);}

	/*End profiler*/
	femmodel->profiler->Stop(SMBCORE);
}

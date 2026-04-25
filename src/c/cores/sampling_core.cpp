/*!\file: sampling_core.cpp
 * \brief: core of the sampling solution
 */

#include "../classes/classes.h"
#include "../solutionsequences/solutionsequences.h"
#include "../analyses/analyses.h" // new
#include "../modules/modules.h"

void sampling_core(FemModel* femmodel){

	/*Start profiler*/
	femmodel->profiler->Start(SAMPLINGCORE);

	/*parameters: */
	int    numoutputs;
	bool   save_results;
	int    solution_type;
	char** requested_outputs = NULL;

	/*activate configuration*/
	femmodel->SetCurrentConfiguration(SamplingAnalysisEnum);

	/*recover parameters: */
	femmodel->parameters->FindParam(&save_results,SaveResultsEnum);
	femmodel->parameters->FindParam(&solution_type,SolutionTypeEnum);
	femmodel->parameters->FindParam(&numoutputs,SamplingNumRequestedOutputsEnum);
	if(numoutputs) femmodel->parameters->FindParam(&requested_outputs,&numoutputs,SamplingRequestedOutputsEnum);

	if(VerboseSolution()) _printf0_("   Generating random samples\n");

	/*Generate random sample*/
	SamplingAnalysis* analysis = new SamplingAnalysis();
	femmodel->SetCurrentConfiguration(SamplingAnalysisEnum);

	if(solution_type==TransientSolutionEnum){
		InputDuplicatex(femmodel,SampleEnum,SampleOldEnum);

    int seed;
		femmodel->parameters->FindParam(&seed,SamplingSeedEnum);
		if(seed>=0){
			int step;
			femmodel->parameters->FindParam(&step,StepEnum);
			seed = seed + 13923272*step; // change default seed for transient simulations (by considering an arbitrary shift based on the step number)
			femmodel->parameters->SetParam(seed,SamplingSeedEnum);
		}

	}

	solutionsequence_sampling(femmodel);

	if(solution_type==TransientSolutionEnum){

		InputDuplicatex(femmodel,SampleEnum,SampleNoiseEnum);

		analysis->UpdateTransientSample(femmodel);

	}

	delete analysis;

	if(save_results){
		if(VerboseSolution()) _printf0_("   saving results\n");
		int outputs = SampleEnum;
		femmodel->RequestedOutputsx(&femmodel->results,&outputs,1);
	}

	if(solution_type==SamplingSolutionEnum)femmodel->RequestedDependentsx();

	/*Free resources:*/
	if(numoutputs){for(int i=0;i<numoutputs;i++){xDelete<char>(requested_outputs[i]);} xDelete<char*>(requested_outputs);}

	/*profiler*/
	femmodel->profiler->Stop(SAMPLINGCORE);

}

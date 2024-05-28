/*!\file:  see IssmParallelDirectApplicInterface.h for documentation.  */ 

/*Issm Configuration: {{{*/
#ifdef HAVE_CONFIG_H
#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif
/*}}}*/

#if !defined(_WRAPPERS_) && defined(_HAVE_DAKOTA_) && _DAKOTA_MAJOR_ >= 6

#include "../classes.h"
#include "../../cores/cores.h"
#include "../../modules/modules.h"

namespace SIM {
	IssmParallelDirectApplicInterface::IssmParallelDirectApplicInterface(const Dakota::ProblemDescDB& problem_db, const MPI_Comm& evaluation_comm, int argc, char** argv) :Dakota::DirectApplicInterface(problem_db){ /*{{{*/

		int world_rank;
		ISSM_MPI_Comm_rank(ISSM_MPI_COMM_WORLD,&world_rank);

		/*Build an femmodel if you are a slave, using the corresponding communicator:*/
		if(world_rank!=0){
			femmodel_init= new FemModel(argc,argv,evaluation_comm);
			femmodel_init->profiler->Start(CORE);
		}

	}
	/*}}}*/
	IssmParallelDirectApplicInterface::~IssmParallelDirectApplicInterface(){ /*{{{*/

		int world_rank;
		ISSM_MPI_Comm_rank(ISSM_MPI_COMM_WORLD,&world_rank);

		if(world_rank!=0){

			/*Wrap up: */
			femmodel_init->profiler->Stop(CORE);
			femmodel_init->CleanUp(); //only close file pointers on rank 0 of slave 1!

			/*Delete Model: */
			delete femmodel_init;
		}
	}
	/*}}}*/
	int IssmParallelDirectApplicInterface::derived_map_ac(const Dakota::String& ac_name){/*{{{*/

		FemModel* femmodel;

		char     **responses_descriptors    = NULL;      //these are ours! there are only numresponsedescriptors of them, not d_numresponses!!!
		char      *response_descriptor      = NULL;
		int        numresponsedescriptors;
		int        solution_type;
		bool       control_analysis         = false;
		void     (*solutioncore)(FemModel*) = NULL;
		void     (*solutionprecore)(FemModel*) = NULL;
		bool       nodakotacore             = true;

		int world_rank;
		ISSM_MPI_Comm_rank(ISSM_MPI_COMM_WORLD,&world_rank);

		/*Only have slaves work!:*/
		if(world_rank==0)return 0;

		#ifdef MPI_DEBUG
		_printf0_("eval server id" << evalServerId << " invoking " << ac_name << " within SIM::IssmParallelDirectApplicInterface." << std::endl);
		_printf0_("evalServerId " << evalServerId << "evaluation_id " << currEvalId <<  "\n");
		#endif // MPI_DEBUG

		int i;
		IssmDouble  *variables            = NULL;
		char       **variable_descriptors = NULL;
		char        *variable_descriptor  = NULL;
		IssmDouble  *responses            = NULL;

		/*Before launching evaluation, we need to transfer the dakota inputs into Issm readable variables: */

		/*First, the variables: */
		variables=xNew<IssmDouble>(numACV);
		for(i=0;i<numACV;i++){
			variables[i]=xC[i];
		}
		/*The descriptors: */
		variable_descriptors=xNew<char*>(numACV);
		for(i=0;i<numACV;i++){
			std::string label=xCLabels[i];
			variable_descriptor=xNew<char>(strlen(label.c_str())+1);
			memcpy(variable_descriptor,label.c_str(),(strlen(label.c_str())+1)*sizeof(char));

			variable_descriptors[i]=variable_descriptor;
		}

		/*Initialize responses: */
		responses=xNewZeroInit<IssmDouble>(numFns);

		/*Launch cores that are not used during the uncertainty quantification: */
		WrapperPreCorePointerFromSolutionEnum(&solutionprecore,femmodel_init->parameters,solution_type);
		if(solutionprecore)solutionprecore(femmodel_init);

		/*Make a copy of femmodel, so we start this new evaluation run for this specific sample with a brand 
		 * new copy of the model, which has not been tempered with by previous evaluation runs: */
		femmodel=femmodel_init->copy();

		/*retrieve parameters: */
		femmodel->parameters->FindParam(&responses_descriptors,&numresponsedescriptors,QmuResponsedescriptorsEnum);
		femmodel->parameters->FindParam(&solution_type,SolutionTypeEnum);
		femmodel->parameters->FindParam(&control_analysis,InversionIscontrolEnum);

		/*include currEvalId in parameters:*/
		femmodel->parameters->SetParam(currEvalId,QmuCurrEvalIdEnum);

		/*Modify core inputs in objects contained in femmodel, to reflect the dakota variables inputs: */
		InputUpdateFromDakotax(femmodel,variables,variable_descriptors,numACV);

		/*Determine solution sequence: */
		if(VerboseQmu()) _printf0_("Starting " << EnumToStringx(solution_type) << " core:\n");
		WrapperCorePointerFromSolutionEnum(&solutioncore,femmodel->parameters,solution_type,nodakotacore);

		/*Run the core solution sequence: */
		solutioncore(femmodel);

		/*compute responses: */
		if(VerboseQmu()) _printf0_("compute dakota responses:\n");
		femmodel->DakotaResponsesx(responses,responses_descriptors,numresponsedescriptors,numFns);

		/*Output results for this iteration: */
		if(VerboseQmu()) _printf0_("output results for this iteration: \n");
		OutputResultsx(femmodel);

		/*populate responses: */
		for(i=0;i<numFns;i++){
			fnVals[i]=responses[i];
		}

		/*Free resources:*/
		xDelete<IssmDouble>(variables);
		for(i=0;i<numACV;i++){
			variable_descriptor=variable_descriptors[i];
			xDelete<char>(variable_descriptor);
		}
		xDelete<char*>(variable_descriptors);
		for(i=0;i<numresponsedescriptors;i++){
			response_descriptor=responses_descriptors[i];
			xDelete<char>(response_descriptor);
		}
		if(responses_descriptors) xDelete<char*>(responses_descriptors);
		xDelete<IssmDouble>(responses);
		delete femmodel;

		return 0;
	}/*}}}*/
}
#endif

/*!\file:  dakota_core.cpp
 * \brief: wrapper to the Dakota capabilities. qmu fires up Dakota, and registers a Dakota Plugin
 * which will be in charge of running the solution sequences repeatedly, to garner statistics.
 *
 * This routine deals with running ISSM and Dakota in library mode. In library mode, Dakota does not
 * run as an executable. Its capabilities are linked into the ISSM software. ISSM calls Dakota routines
 * directly from the Dakota library. qmu.cpp is the code that is in charge of calling those routines.
 *
 * Dakota has its own way of running in parallel (for embarrassingly parallel jobs). We do not want that,
 * as ISSM knows exactly how to run "really parallel" jobs that use all CPUS. To bypass Dakota's parallelism,
 * we overloaded the constructor for the parallel library (see the Dakota patch in the externalpackages/dakota
 * directory). This overloaded constructor fires up Dakota serially on CPU 0 only! We take care of broadcasting
 * to the other CPUS, hence ISSM is running in parallel, and Dakota serially on CPU0.
 *
 * Now, how does CPU 0 drive all other CPUS to carry out sensitivity analyses? By synchronizing its call to
 * our ISSM cores (stressbalance_core, thermal_core, transient_core, etc ...) on CPU 0 with all other CPUS.
 * This explains the structure of qmu.cpp, where cpu 0 runs Dakota, the Dakota pluggin fires up DakotaSpawnCore.cpp,
 * while the other CPUS are waiting for a broadcast from CPU0, once they get it, they also fire up
 * DakotaSpawnCore. In the end, DakotaSpawnCore is fired up on all CPUS, with CPU0 having Dakota inputs, that it will
 * broadcast to other CPUS.
 *
 * Now, how does Dakota call the DakotaSpawnCore routine? The DakotaSpawnCore is embedded into the DakotaPlugin object
 * which is derived from the Direct Interface Dakota object. This is the only way to run Dakota in library
 * mode (see their developer guide for more info). Dakota registers the DakotaPlugin object into its own
 * database, and calls on the embedded DakotaSpawnCore from CPU0.
 *
 */

 /* \brief: run core ISSM solution using Dakota inputs coming from CPU 0.
 * \sa qmu.cpp DakotaPlugin.cpp
 *
 * This routine needs to be understood simultaneously with qmu.cpp and DakotaPlugin.
 * DakotaSpawnCoreParallel is called by all CPUS, with CPU 0 holding Dakota variable values, along
 * with variable descriptors.
 *
 * DakotaSpawnCoreParallel takes care of broadcasting the variables and their descriptors across the MPI
 * ring. Once this is done, we use the variables to modify the inputs for the solution core.
 * For ex, if "rho_ice" is provided, for ex 920, we include "rho_ice" in the inputs, then
 * call the core with the modified inputs. This is the way we get Dakota to explore the parameter
 * space of the core.
 *
 * Once the core is called, we process the results of the core, and using the processed results,
 * we compute response functions. The responses are computed on all CPUS, but they are targeted
 * for CPU 0, which will get these values back to the Dakota engine.
 *
 */

/*include config: {{{*/
#ifdef HAVE_CONFIG_H
#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif
/*}}}*/

/*include ISSM files: */
#include "./cores.h"
#include "../toolkits/toolkits.h"
#include "../shared/shared.h"
#include "../classes/classes.h"
#include "../modules/modules.h"

#if defined(_HAVE_DAKOTA_) && (_DAKOTA_MAJOR_ <= 5) //this only works for Dakota <=5, which had no effective parallel capabilities yet.
/*Dakota include files:{{{*/
#if (_DAKOTA_MAJOR_ < 5 || (_DAKOTA_MAJOR_ == 5 && _DAKOTA_MINOR_ < 3))
#include <ParallelLibrary.H>
#include <ProblemDescDB.H>
#include <DakotaStrategy.H>
#include <DakotaModel.H>
#include <DakotaInterface.H>
#else
#include <ParallelLibrary.hpp>
#include <ProblemDescDB.hpp>
#include <DakotaStrategy.hpp>
#include <DakotaModel.hpp>
#include <DakotaInterface.hpp>
#endif
/*}}}*/

void DakotaFree(double** pvariables,char*** pvariables_descriptors,char*** presponses_descriptors,int numvariables,int numresponses){ /*{{{*/

	/*\brief DakotaFree: free allocations on other CPUs, not done by Dakota.*/

	int i;
	int my_rank;

	double  *variables             = NULL;
	char   **variables_descriptors = NULL;
	char   **responses_descriptors = NULL;
	char    *string                = NULL;

	/*recover pointers: */
	variables=*pvariables;
	variables_descriptors=*pvariables_descriptors;
	responses_descriptors=*presponses_descriptors;

	/*recover my_rank:*/
	my_rank=IssmComm::GetRank();

	/*Free variables and variables_descriptors only on CPU !=0*/
	if(my_rank!=0){
		xDelete<double>(variables);
		for(i=0;i<numvariables;i++){
			string=variables_descriptors[i];
			xDelete<char>(string);
		}
		xDelete<char*>(variables_descriptors);
	}

	//responses descriptors on every CPU
	for(i=0;i<numresponses;i++){
		string=responses_descriptors[i];
		xDelete<char>(string);
	}
	//rest of dynamic allocations.
	xDelete<char*>(responses_descriptors);

	/*Assign output pointers:*/
	*pvariables=variables;
	*pvariables_descriptors=variables_descriptors;
	*presponses_descriptors=responses_descriptors;
} /*}}}*/
void DakotaMPI_Bcast(double** pvariables, char*** pvariables_descriptors,int* pnumvariables, int* pnumresponses){ /*{{{*/

	/* * \brief: broadcast variables_descriptors, variables, numvariables and numresponses
	 * from CPU 0 to all other CPUs.
	 */

	int i;
	int my_rank;

	/*inputs and outputs: */
	double* variables=NULL;
	char**  variables_descriptors=NULL;
	int     numvariables;
	int     numresponses;

	/*intermediary: */
	char* string=NULL;
	int   string_length;

	/*recover my_rank:*/
	my_rank=IssmComm::GetRank();

	/*recover inputs from pointers: */
	variables=*pvariables;
	variables_descriptors=*pvariables_descriptors;
	numvariables=*pnumvariables;
	numresponses=*pnumresponses;

	/*numvariables: */
	ISSM_MPI_Bcast(&numvariables,1,ISSM_MPI_INT,0,IssmComm::GetComm());

	/*variables:*/
	if(my_rank!=0)variables=xNew<double>(numvariables);
	ISSM_MPI_Bcast(variables,numvariables,MPI_DOUBLE,0,IssmComm::GetComm());

	/*variables_descriptors: */
	if(my_rank!=0){
		variables_descriptors=xNew<char*>(numvariables);
	}
	for(i=0;i<numvariables;i++){
		if(my_rank==0){
			string=variables_descriptors[i];
			string_length=(strlen(string)+1)*sizeof(char);
		}
		ISSM_MPI_Bcast(&string_length,1,ISSM_MPI_INT,0,IssmComm::GetComm());
		if(my_rank!=0)string=xNew<char>(string_length);
		ISSM_MPI_Bcast(string,string_length,ISSM_MPI_CHAR,0,IssmComm::GetComm());
		if(my_rank!=0)variables_descriptors[i]=string;
	}

	/*numresponses: */
	ISSM_MPI_Bcast(&numresponses,1,ISSM_MPI_INT,0,IssmComm::GetComm());

	/*Assign output pointers:*/
	*pnumvariables=numvariables;
	*pvariables=variables;
	*pvariables_descriptors=variables_descriptors;
	*pnumresponses=numresponses;
} /*}}}*/
int  DakotaSpawnCore(double* d_responses, int d_numresponses, double* d_variables, char** d_variables_descriptors,int d_numvariables, void* void_femmodel,int counter){ /*{{{*/

	/*Notice the d_, which prefixes anything that is being provided to us by the Dakota plugin. Careful: some things are ours; some are DDkota's!: */

	char     **responses_descriptors    = NULL;      //these are our! there are only numresponsedescriptors of them, not d_numresponses!!!
	int        numresponsedescriptors;
	int        solution_type;
	bool       control_analysis         = false;
	void     (*solutioncore)(FemModel*) = NULL;
	FemModel  *femmodel                 = NULL;
	bool       nodakotacore             = true;

	/*If counter==-1 on CPU 0, it means that the Dakota runs are done. In which case, bail out and return 0: */
	ISSM_MPI_Bcast(&counter,1,ISSM_MPI_INT,0,IssmComm::GetComm());
	if(counter==-1)return 0;

	/*cast void_femmodel to FemModel, and at the same time, make a copy, so we start this new core run for this specific sample
	 *with a brand new copy of the model, which has not been tempered with by previous Dakota runs: */
	femmodel=(reinterpret_cast<FemModel*>(void_femmodel))->copy();

	/*retrieve parameters: */
	femmodel->parameters->FindParam(&responses_descriptors,&numresponsedescriptors,QmuResponsedescriptorsEnum);
	femmodel->parameters->FindParam(&solution_type,SolutionTypeEnum);
	femmodel->parameters->FindParam(&control_analysis,InversionIscontrolEnum);

	if(VerboseQmu()) _printf0_("qmu iteration: " << counter << "\n");

	/* only CPU 0, running Dakota is providing us with variables and variables_descriptors and numresponses: broadcast onto other CPUs: */
	DakotaMPI_Bcast(&d_variables,&d_variables_descriptors,&d_numvariables,&d_numresponses);

	/*Modify core inputs in objects contained in femmodel, to reflect the dakota variables inputs: */
	InputUpdateFromDakotax(femmodel,d_variables,d_variables_descriptors,d_numvariables);

	/*Determine solution sequence: */
	if(VerboseQmu()) _printf0_("Starting " << EnumToStringx(solution_type) << " core:\n");
	WrapperCorePointerFromSolutionEnum(&solutioncore,femmodel->parameters,solution_type,nodakotacore);

	/*Run the core solution sequence: */
	solutioncore(femmodel);

	/*compute responses: */
	if(VerboseQmu()) _printf0_("compute dakota responses:\n");
	femmodel->DakotaResponsesx(d_responses,responses_descriptors,numresponsedescriptors,d_numresponses);

	/*output for this core:*/
	if(VerboseQmu()) _printf0_("outputing results for this core:\n");
	OutputResultsx(femmodel);

	/*Free resources:*/
	DakotaFree(&d_variables,&d_variables_descriptors,&responses_descriptors, d_numvariables, numresponsedescriptors);

	/*Avoid leaks here: */
	delete femmodel;

	return 1; //this is critical! do not return 0, otherwise, dakota_core will stop running!
}
/*}}}*/
void dakota_core(FemModel* femmodel){  /*{{{*/

	int                my_rank;
	char              *dakota_input_file  = NULL;
	char              *dakota_output_file = NULL;
	char              *dakota_error_file  = NULL;
	Dakota::ModelLIter ml_iter;

	/*Recover dakota_input_file, dakota_output_file and dakota_error_file, in the parameters dataset in parallel */
	femmodel->parameters->FindParam(&dakota_input_file,QmuInNameEnum);
	femmodel->parameters->FindParam(&dakota_output_file,QmuOutNameEnum);
	femmodel->parameters->FindParam(&dakota_error_file,QmuErrNameEnum);

	/*recover my_rank:*/
	my_rank=IssmComm::GetRank();

	if(my_rank==0){

		// Instantiate/initialize the parallel library and problem description
		// database objects.
		char* dakotamode=xNew<char>(strlen("serial")+1);
		xMemCpy<char>(dakotamode,"serial",strlen("serial")+1);
		Dakota::ParallelLibrary parallel_lib(dakotamode); //use our own ISSM Dakota library mode constructor, which only fires up Dakota on CPU 0.
		Dakota::ProblemDescDB problem_db(parallel_lib);
		xDelete<char>(dakotamode);

		// Manage input file parsing, output redirection, and restart processing
		// without a CommandLineHandler.  This version relies on parsing of an
		// input file.
		problem_db.manage_inputs(dakota_input_file);
		// specify_outputs_restart() is only necessary if specifying non-defaults
		parallel_lib.specify_outputs_restart(dakota_output_file,dakota_error_file,NULL,NULL);

		// Instantiate the Strategy object (which instantiates all Model and
		// Iterator objects) using the parsed information in problem_db.
		Dakota::Strategy selected_strategy(problem_db);

		// convenience function for iterating over models and performing any
		// interface plug-ins
		Dakota::ModelList& models = problem_db.model_list();

		for (ml_iter = models.begin(); ml_iter != models.end(); ml_iter++) {

			Dakota::Interface& interface = ml_iter->interface();

			//set DB nodes to the existing Model specification
			problem_db.set_db_model_nodes(ml_iter->model_id());

			// Serial case: plug in derived Interface object without an analysisComm
			interface.assign_rep(new SIM::IssmDirectApplicInterface(problem_db,(void*)femmodel), false);
		}

		// Execute the strategy
		problem_db.lock(); // prevent run-time DB queries
		selected_strategy.run_strategy();

		//Warn other CPUs that we are done running the Dakota iterator, by setting the counter to -1:
		DakotaSpawnCore(NULL,0, NULL,NULL,0,femmodel,-1);

	}
	else{

		for(;;){
			if(!DakotaSpawnCore(NULL,0, NULL,NULL,0,femmodel,0)) break; //counter came in at -1 on CPU 0, bail out.
		}
	}

	/*Free resources:*/
	xDelete<char>(dakota_input_file);
	xDelete<char>(dakota_error_file);
	xDelete<char>(dakota_output_file);

} /*}}}*/
#else
void dakota_core(FemModel* femmodel){
	_error_("dakota_core for versions of Dakota >=6 should not be used anymore! Use instead the issm_dakota executable!");
}
#endif

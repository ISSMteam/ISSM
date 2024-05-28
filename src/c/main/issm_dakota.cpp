/*!\file:  issm_dakota.cpp
 * \brief: ISSM DAKOTA main program
 */ 

#include "./issm.h"

/*Dakota includes: */
#if defined(_HAVE_DAKOTA_) && _DAKOTA_MAJOR_ >= 6
#include "ParallelLibrary.hpp"
#include "ProblemDescDB.hpp"
#include "LibraryEnvironment.hpp"
#include "DakotaModel.hpp"
#include "DakotaInterface.hpp"
#endif

int main(int argc,char **argv){ /*{{{*/

	#if defined(_HAVE_DAKOTA_) && _DAKOTA_MAJOR_ >= 6

	bool parallel=true;
	char* dakota_input_file=NULL;
	char* dakota_output_file = NULL;
	char* dakota_error_file = NULL;
	bool statistics=false;

	/*Define MPI_DEBUG in dakota_global_defs.cpp to cause a hold here*/
	Dakota::mpi_debug_hold();

	/*Initialize MPI: */
	ISSM_MPI_Init(&argc, &argv); // initialize MPI

	/*Recover file name for dakota input file:*/
	dakota_input_file=xNew<char>((strlen(argv[2])+strlen(argv[3])+strlen(".qmu.in")+2));
	sprintf(dakota_input_file,"%s/%s%s",argv[2],argv[3],".qmu.in");

	dakota_output_file=xNew<char>((strlen(argv[2])+strlen(argv[3])+strlen(".qmu.out")+2));
	sprintf(dakota_output_file,"%s/%s%s",argv[2],argv[3],".qmu.out");

	dakota_error_file=xNew<char>((strlen(argv[2])+strlen(argv[3])+strlen(".qmu.err")+2));
	sprintf(dakota_error_file,"%s/%s%s",argv[2],argv[3],".qmu.err");

	/*Create directory structure for model outputs:*/
	#if _SYSTEM_HAS_FMEMOPEN_ == 1
	statistics=DakotaDirStructure(argc,argv);
	#endif

	/* Parse input and construct Dakota LibraryEnvironment, performing input data checks*/
	Dakota::ProgramOptions opts;
	opts.input_file(dakota_input_file);
	opts.output_file(dakota_output_file);
	opts.error_file(dakota_error_file);

	/* Defaults constructs the MPIManager, which assumes COMM_WORLD*/
	Dakota::LibraryEnvironment env(opts);

	/* get the list of all models matching the specified model, interface, driver:*/
	Dakota::ModelList filt_models = env.filtered_model_list("single", "direct", "matlab");
	if (filt_models.empty()) {
		Cerr << "Error: no parallel interface plugin performed.  Check compatibility "
			<< "between parallel\n       configuration and selected analysis_driver."
			<< std::endl;
		Dakota::abort_handler(-1);
	}

	Dakota::ProblemDescDB& problem_db = env.problem_description_db();
	Dakota::ModelLIter ml_iter;
	size_t model_index = problem_db.get_db_model_node(); // for restoration
	for (ml_iter = filt_models.begin(); ml_iter != filt_models.end(); ++ml_iter) {
		// set DB nodes to input specification for this Model
		problem_db.set_db_model_nodes(ml_iter->model_id());

		Dakota::Interface& model_interface = ml_iter->derived_interface();

		// Parallel case: plug in derived Interface object with an analysisComm.
		// Note: retrieval and passing of analysisComm is necessary only if
		// parallel operations will be performed in the derived constructor.

		// retrieve the currently active analysisComm from the Model.  In the most
		// general case, need an array of Comms to cover all Model configurations.
		const MPI_Comm& analysis_comm = ml_iter->analysis_comm();

		// don't increment ref count since no other envelope shares this letter
		model_interface.assign_rep(new
				SIM::IssmParallelDirectApplicInterface(problem_db, analysis_comm, argc, argv), false);
	}
	problem_db.set_db_model_nodes(model_index);            // restore

	/* Execute the environment:*/
	env.execute();

	/* Run statistics if requested:*/
	#if _SYSTEM_HAS_FMEMOPEN_ == 1
	if(statistics)DakotaStatistics(argc,argv);
	#endif

	/*free allocations:*/
	xDelete<char>(dakota_input_file);
	xDelete<char>(dakota_output_file);
	xDelete<char>(dakota_error_file);

	/*Return unix success: */
	return 0; 
	#else 
	Cout <<  "ISSM Dakota executable was compiled without support of Dakota! Will just return now!" << "\n";
	return 1;
	#endif

} /*}}}*/

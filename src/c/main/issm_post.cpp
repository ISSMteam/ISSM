/*!\file:  issm_post.cpp
 * \brief: ISSM DAKOTA post-processing of statistics
 */ 

#include "./issm.h"
#include <sys/stat.h>

int main(int argc,char **argv){ /*{{{*/

	char* dakota_input_file=NULL;
	char* dakota_output_file = NULL;
	char* dakota_error_file = NULL;
	bool statistics=false;

	/*Initialize MPI: */
	ISSM_MPI_Init(&argc, &argv); // initialize MPI

	/*Run statistics:*/
	#if _SYSTEM_HAS_FMEMOPEN_ == 1
	DakotaStatistics(argc,argv);
	#endif

	/*Return unix success: */
	return 0; 

} /*}}}*/

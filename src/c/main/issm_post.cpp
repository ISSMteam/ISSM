/*!\file:  issm_post.cpp
 * \brief: ISSM DAKOTA post-processing of statistics
 */ 

#include "./issm.h"
#include <sys/stat.h>

int main(int argc,char **argv){ /*{{{*/

	/*Initialize MPI: */
	ISSM_MPI_Init(&argc, &argv); // initialize MPI

	/*Run statistics:*/
	#if _SYSTEM_HAS_FMEMOPEN_ == 1
	DakotaStatistics(argc,argv);
	#endif

	/*Return unix success: */
	return 0; 

} /*}}}*/

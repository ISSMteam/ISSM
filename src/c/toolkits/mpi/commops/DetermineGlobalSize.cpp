/* \file DetermineGlobalSize.cpp
 * \brief: routine to determine global size from local size 
 */

#include <stdio.h>
#include <math.h>
#include "../../../shared/shared.h"
#include "../../../shared/Numerics/types.h"

int DetermineGlobalSize(int local_size,ISSM_MPI_Comm comm){

	/*output: */
	int  global_size;

	ISSM_MPI_Reduce(&local_size, &global_size, 1, ISSM_MPI_INT, ISSM_MPI_SUM, 0, comm);
	ISSM_MPI_Bcast(&global_size,1,ISSM_MPI_INT,0,comm);

	return global_size;

}

/* \file DetermineLocalSize.cpp
 * \brief: routine to determine local size of a global petsc matrix or vector. 
 */

#include <stdio.h>
#include <math.h>
#include "../../../shared/shared.h"

int DetermineLocalSize(int global_size,ISSM_MPI_Comm comm){

	/*output: */
	int  local_size;

	/*intermediary: */
	int  i;
	int  row_rest;
	int* num_local_rows=NULL;

	/*from MPI: */
	int num_procs;
	int my_rank;

	/*recover my_rank*/
	ISSM_MPI_Comm_rank(comm,&my_rank);
	ISSM_MPI_Comm_size(comm,&num_procs);

	/* TODO replace the following with ->

	local_size=global_size/num_procs; // integer division
	if (global_size%num_procs>my_rank) local_size++; // distribute the remainder
	return local_size;

	<- to here  */

	/*We are  not bound by any library, just use what seems most logical*/
	num_local_rows=xNew<int>(num_procs);    

	for (i=0;i<num_procs;i++){
		/*Here, we use floor. We under distribute rows. The rows 
		  left  are then redistributed, therefore resulting in a 
		  more even distribution.*/
		num_local_rows[i]=(int)floor((double)global_size/(double)num_procs); 
	}

	/*There may be some rows left. Distribute evenly.*/ 
	row_rest=global_size - num_procs*(int)floor((double)global_size/(double)num_procs);
	for (i=0;i<row_rest;i++){
		num_local_rows[i]++;
	}
	local_size=num_local_rows[my_rank];

	/*Free resources: */
	xDelete<int>(num_local_rows);

	/*return size: */
	return local_size;

}

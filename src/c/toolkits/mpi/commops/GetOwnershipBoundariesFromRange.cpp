/*! \file:  GetOwnershipBoundariesFromRange.cpp
 *  \brief from a local range on each cpu, we determine what 
 *  lower row and upper row from a matrix a cpu owns.
 */

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include <stdio.h>
#include "../../../shared/MemOps/MemOps.h"
#include "../../../shared/io/Comm/IssmComm.h"

void GetOwnershipBoundariesFromRange(int* plower_row,int* pupper_row,int range,ISSM_MPI_Comm comm){

	/*externals :*/
	int my_rank;
	int num_procs;

	/*recover my_rank and num_procs:*/
	ISSM_MPI_Comm_size(comm,&num_procs);
	ISSM_MPI_Comm_rank(comm,&my_rank);

	/*output: */
	int lower_row,upper_row;

	/*Gather all range values into allranges, for all nodes*/
	int* allranges=xNew<int>(num_procs);
	ISSM_MPI_Allgather(&range,1,ISSM_MPI_INT,allranges,1,ISSM_MPI_INT,comm);

	/*From all ranges, get lower row and upper row*/
	lower_row=0;
	upper_row=lower_row+allranges[0];
	for(int i=1;i<=my_rank;i++){
		lower_row=lower_row+allranges[i-1];
		upper_row=upper_row+allranges[i];
	}

	/*Assign output*/
	*plower_row=lower_row;
	*pupper_row=upper_row;
	xDelete<int>(allranges);
}

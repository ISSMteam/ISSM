/* \file DetermineRowRankFromLocalSize.cpp
 * \brief: routine to determine, from local size of a matrix or vector (in terms of number 
 * of rows), a vector of global size  which for each row of the matrix or vector, determines
 * what cpu this row belong to.
 */

#include <stdio.h>
#include <math.h>
#include "../../../shared/shared.h"
#include "../../../shared/Numerics/types.h"

int* DetermineRowRankFromLocalSize(int global_size,int localsize,ISSM_MPI_Comm comm){

	/*intermediary: */
	int i,j;
	int my_rank=0;
	int num_procs=0;
	int lower_row,upper_row;

	/*output: */
	int* RowRank=NULL;

	ISSM_MPI_Comm_rank(comm,&my_rank);
	ISSM_MPI_Comm_size(comm,&num_procs);

	/*allocate: */
	RowRank=xNew<int>(global_size);

	/*Gather all local_size values into alllocalsizes, for all cpus*/
	int* alllocalsizes=xNew<int>(num_procs);
	ISSM_MPI_Allgather(&localsize,1,ISSM_MPI_INT,alllocalsizes,1,ISSM_MPI_INT,comm);

	/*From all localsizes, get lower row and upper row*/
	lower_row=0;
	upper_row=lower_row+alllocalsizes[0];
	for(j=lower_row;j<upper_row;j++)RowRank[j]=0;

	for(i=1;i<num_procs;i++){
		lower_row=lower_row+alllocalsizes[i-1];
		upper_row=upper_row+alllocalsizes[i];
		for(j=lower_row;j<upper_row;j++)RowRank[j]=i;
	}
	xDelete<int>(alllocalsizes);

	return RowRank;
}

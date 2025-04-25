/*!\file MatToMPISerial.cpp
 * \brief gather a Petsc Mat matrix onto all cpus
 */

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "petscpatches.h"
#include "../petscincludes.h"
#include "../../../shared/shared.h"

void MatToMPISerial(IssmDouble** poutmatrix,PMat matrix,ISSM_MPI_Comm comm,bool broadcast){

	int i;
	int my_rank;
	int num_procs;

	/*Petsc variables*/
	PetscInt lower_row,upper_row; 
	int range;
	int M,N; //size of matrix
	ISSM_MPI_Status status;
	int* idxm=NULL;
	int* idxn=NULL; 
	IssmDouble* local_matrix=NULL; /*matrix local to each node used for temporary holding matrix values*/
	int buffer[3];

	/*recover my_rank and num_procs:*/
	ISSM_MPI_Comm_size(comm,&num_procs);
	ISSM_MPI_Comm_rank(comm,&my_rank);

	/*Output*/
	IssmDouble* outmatrix=NULL;

	/*get matrix size: */
	MatGetSize(matrix,&M,&N);

	/*partition: */
	MatGetOwnershipRange(matrix,&lower_row,&upper_row);    
	upper_row--; 
	range=upper_row-lower_row+1;

	/*Local and global allocation*/
	if(broadcast || my_rank==0){ 
		outmatrix=xNew<IssmDouble>(M*N);
	}

	if (range){
		local_matrix=xNew<IssmDouble>(N*range);
		idxm=xNew<int>(range);  
		idxn=xNew<int>(N);  

		for (i=0;i<N;i++){
			*(idxn+i)=i;
		}
		for (i=0;i<range;i++){
			*(idxm+i)=lower_row+i;
		}

		MatGetValues(matrix,range,idxm,N,idxn,local_matrix);     
	}

	/*Now each node holds its local_matrix containing range rows. 
	 * We send these rows to the matrix on node 0*/

	for (i=1;i<num_procs;i++){
		if (my_rank==i){ 
			buffer[0]=my_rank;
			buffer[1]=lower_row;
			buffer[2]=range;
			ISSM_MPI_Send(buffer,3,ISSM_MPI_INT,0,1,comm);   
			if (range)ISSM_MPI_Send(local_matrix,N*range,ISSM_MPI_PDOUBLE,0,1,comm); 
		}
		if (my_rank==0){
			ISSM_MPI_Recv(buffer,3,ISSM_MPI_INT,i,1,comm,&status); 
			if (buffer[2])ISSM_MPI_Recv(outmatrix+(buffer[1]*N),N*buffer[2],ISSM_MPI_PDOUBLE,i,1,comm,&status);
		}
	} 
	if (my_rank==0){ 
		//Still have the local_matrix on node 0 to take care of.
		if (range) memcpy(outmatrix,local_matrix,N*range*sizeof(IssmDouble));
	} 

	if(broadcast){
		/*Broadcast:*/
		ISSM_MPI_Bcast(outmatrix,M*N,ISSM_MPI_PDOUBLE,0,comm);
	}

	/*Assign output pointer: */
	xDelete<int>(idxm);
	xDelete<int>(idxn);
	xDelete<IssmDouble>(local_matrix);
	*poutmatrix=outmatrix;
}

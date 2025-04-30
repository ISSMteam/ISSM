/*!\file VecToMPISerial.cpp
 * \brief gather a Petsc Vector spread across the cluster, onto node 0, and then broadcast to all nodes. 
 */

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "../petscincludes.h"
#include "../../../shared/shared.h"

//template<typename doubletype, typename vectype>
//int VecToMPISerialNew(doubletype** pserial_vector, vectype vector,ISSM_MPI_Comm comm,bool broadcast){
//
//	/*Output*/
//	const doubletype *vec_array     = NULL;
//	doubletype       *serial_vector = NULL;
//
//	/*Sequential Vector*/
//	int        n;
//	vectype    vector_seq = NULL;
//	VecScatter ctx        = NULL;
//
//	if(broadcast){
//		VecScatterCreateToAll(vector, &ctx, &vector_seq);
//	}
//	else{
//		VecScatterCreateToZero(vector, &ctx, &vector_seq);
//	}
//
//  /*scatter as many times as you need*/
//  VecScatterBegin(ctx, vector, vector_seq, INSERT_VALUES, SCATTER_FORWARD);
//  VecScatterEnd(  ctx, vector, vector_seq, INSERT_VALUES, SCATTER_FORWARD);
//
//  /*Get pointer to array and copy*/
//  VecGetArrayRead(vector_seq, &vec_array);
//
//  /* Use memcpy to copy data*/
//  VecGetSize(vector_seq, &n);
//  memcpy(serial_vector, vec_array, n*sizeof(doubletype));
//
//  /* Restore and destroy the PETSc Vec array*/
//  VecRestoreArrayRead(vector_seq, &vec_array);
//
//  /* destroy scatter context and local vector when no longer needed*/
//  VecScatterDestroy(&ctx);
//  VecDestroy(&vector_seq);
//
//  /*Assign output pointer*/
//  *pserial_vector = serial_vector;
//}

template<typename doubletype, typename vectype>
int VecToMPISerial(doubletype** pgathered_vector, vectype vector,ISSM_MPI_Comm comm,bool broadcast){

	int i;
	int num_procs; 
	int my_rank;

	/*Petsc*/
	ISSM_MPI_Status status;
	PetscInt lower_row,upper_row; 
	int range;
	int * idxn=NULL; 
	int buffer[3];

	/*intermediary results*/
	doubletype* local_vector=NULL;

	/*input*/
	int vector_size;

	/*Output*/
	doubletype* gathered_vector=NULL; //Global vector holding the final assembled vector on all nodes.

	/*recover my_rank and num_procs*/
	ISSM_MPI_Comm_size(comm,&num_procs);
	ISSM_MPI_Comm_rank(comm,&my_rank);

	VecGetSize(vector,&vector_size);
	if(vector_size==0){
		*pgathered_vector=NULL;
		return 1;
	}

	/*Allocate gathered vector on all nodes .*/
	if(broadcast || my_rank==0){ 
		gathered_vector=xNew<doubletype>(vector_size);
	}

	/*Allocate local vectors*/
	VecGetOwnershipRange(vector,&lower_row,&upper_row);
	upper_row--;
	range=upper_row-lower_row+1;    

	if (range){
		idxn=xNew<int>(range); 
		for (i=0;i<range;i++){
			*(idxn+i)=lower_row+i;
		} 
		local_vector=xNew<doubletype>(range);
		/*Extract values from MPI vector to serial local_vector on each node*/
		VecGetValues(vector,range,idxn,local_vector); 
	}

	/*Now each node holds its local_vector containing range rows. 
	 * We send this local_vector  to the gathered_vector on node 0*/
	for (i=1;i<num_procs;i++){
		if (my_rank==i){ 
			buffer[0]=my_rank;
			buffer[1]=lower_row;
			buffer[2]=range;
			ISSM_MPI_Send(buffer,3,ISSM_MPI_INT,0,1,comm);
			if (range)ISSM_MPI_Send(local_vector,range,TypeToMPIType<doubletype>(),0,1,comm);
		}
		if (my_rank==0){
			ISSM_MPI_Recv(buffer,3,ISSM_MPI_INT,i,1,comm,&status); 
			if (buffer[2])ISSM_MPI_Recv(gathered_vector+buffer[1],buffer[2],TypeToMPIType<doubletype>(),i,1,comm,&status);
		}
	}

	if (my_rank==0){ 
		//Still have the local_vector on node 0 to take care of.
		if (range) {
			xMemCpy<doubletype>(&gathered_vector[lower_row], local_vector, range);
		}
	}

	if(broadcast){
		/*Now, broadcast gathered_vector from node 0 to other nodes: */
		ISSM_MPI_Bcast(gathered_vector,vector_size,TypeToMPIType<doubletype>(),0,comm);
	}

	/*Assign output pointers: */
	*pgathered_vector=gathered_vector;

	/*Free resources: */
	xDelete<int>(idxn);
	xDelete<doubletype>(local_vector);

	return 1;
}

//template int VecToMPISerialNew(IssmDouble** pserial_vector, PVec vector,ISSM_MPI_Comm comm,bool broadcast);
template int VecToMPISerial(IssmDouble** pgathered_vector, PVec vector,ISSM_MPI_Comm comm,bool broadcast);
#if _HAVE_CODIPACK_
//template int VecToMPISerialNew(IssmPDouble** pserial_vector, Vec vector,ISSM_MPI_Comm comm,bool broadcast);
template int VecToMPISerial(IssmPDouble** pgathered_vector, Vec vector,ISSM_MPI_Comm comm,bool broadcast);
#endif

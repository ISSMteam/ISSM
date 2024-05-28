/*!\file:  MatMultPatch
 * \brief: relocalize vector when MatMult yields non conforming object sizes errors.
 */ 

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

/*Petsc includes: */
#include <petscksp.h>

#include "../../mpi/issmmpi.h"
#include "../../../shared/shared.h"

/*Function prototypes: */
int MatMultCompatible(Mat A,Vec x,ISSM_MPI_Comm comm);
void VecRelocalize(Vec* outvector,Vec vector,int m,ISSM_MPI_Comm comm);

void MatMultPatch(Mat A,Vec X, Vec AX,ISSM_MPI_Comm comm){ //same prototype as MatMult in Petsc

	int m,n;
	Vec X_rel=NULL;

	_assert_(A); _assert_(X);

	if (MatMultCompatible(A,X,comm)){
		MatMult(A,X,AX); 
	}
	else{
		MatGetLocalSize(A,&m,&n);;
		VecRelocalize(&X_rel,X,n,comm);
		MatMult(A,X_rel,AX); ;
		#if PETSC_VERSION_LT(3,2,0)
		VecDestroy(X_rel);
		#else
		VecDestroy(&X_rel);
		#endif
	}
}

int MatMultCompatible(Mat A,Vec x,ISSM_MPI_Comm comm){

	/*error management*/

	int local_m,local_n;
	int range;
	int result=1;
	int sumresult;
	int num_procs;

	/*recover num_procs:*/
	ISSM_MPI_Comm_size(comm,&num_procs);

	MatGetLocalSize(A,&local_m,&local_n);;
	VecGetLocalSize(x,&range);;

	if (local_n!=range)result=0;

	/*synchronize result: */
	ISSM_MPI_Reduce (&result,&sumresult,1,ISSM_MPI_INT,ISSM_MPI_SUM,0,comm );
	ISSM_MPI_Bcast(&sumresult,1,ISSM_MPI_INT,0,comm);                
	if (sumresult!=num_procs){
		result=0;
	}
	else{
		result=1;
	}
	return result;
}

void VecRelocalize(Vec* poutvector,Vec vector,int m,ISSM_MPI_Comm comm){

	/*vector index and vector values*/
	int* index=NULL;
	double* values=NULL;
	int lower_row,upper_row,range;

	/*output: */
	Vec outvector=NULL;

	/*Create outvector with local size m*/
	VecCreate(comm,&outvector); ; 
	VecSetSizes(outvector,m,PETSC_DECIDE); ; 
	VecSetFromOptions(outvector); ; 

	/*Go through vector, get values, and plug them into outvector*/
	VecGetOwnershipRange(vector,&lower_row,&upper_row); ; 
	upper_row--;
	range=upper_row-lower_row+1;
	if (range){
		index=xNew<int>(range);
		values=xNew<double>(range);
		for (int i=0;i<range;i++){
			*(index+i)=lower_row+i;
		}
		VecGetValues(vector,range,index,values);
		VecSetValues(outvector,range,index,values,INSERT_VALUES);
	}

	/*Assemble outvector*/
	VecAssemblyBegin(outvector);; 
	VecAssemblyEnd(outvector);; 

	/*Free resources:*/
	xDelete<int>(index);
	xDelete<double>(values);	

	/*Assign output pointers:*/
	*poutvector=outvector;

}

/*!\file PetscMat.cpp
 * \brief: implementation of the PetscMat object
 */

/*Headers:*/
/*{{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include <stdio.h>
#include <string.h>
#include "../petscincludes.h"
#include "../../../shared/shared.h"

/*}}}*/

/*PetscMat constructors and destructor*/
PetscMat::PetscMat(){/*{{{*/
	this->matrix=NULL;
	#ifdef _HAVE_AD_
	this->amatrix=NULL;
	#endif

}
/*}}}*/
PetscMat::PetscMat(int M,int N){/*{{{*/

	this->matrix=NewMat(M,N,IssmComm::GetComm());
}
/*}}}*/
PetscMat::PetscMat(int M,int N, IssmDouble sparsity){/*{{{*/

	this->matrix=NewMat(M,N,sparsity,IssmComm::GetComm());
}
/*}}}*/
PetscMat::PetscMat(int m,int n,int M,int N,int* d_nnz,int* o_nnz){/*{{{*/

	MatCreate(IssmComm::GetComm(),&this->matrix);
	MatSetSizes(this->matrix,m,n,M,N);
	MatSetFromOptions(this->matrix);

	/* 
	 * Versions of Petsc beyond 3.3 have changed the use of preallocation 
	 * routines to distinguish between parallel builds and sequential. Since
	 * our Windows builds are currently only sequential, we need to change
	 * the way we use these functions.
	 *
	 * The following code computes the total number of non-zeroes per row of the
	 * matrix in question. In parallel builds it is nescessary to kep track of 
	 * diagonal non zeros and off-diagonal (d_nnz and o_nnz). Sequential does
	 * not make that distinction.
	*/
	#ifdef _HAVE_PETSC_MPI_
		int* nnz = new int[M];
		for(int i = 0; i < M; i++)
			nnz[i] = o_nnz[i] + d_nnz[i];

		PetscErrorCode ierr = MatSeqAIJSetPreallocation(this->matrix,0,nnz);
		delete[] nnz;
	#else
		PetscErrorCode ierr = MatMPIAIJSetPreallocation(this->matrix,0,d_nnz,0,o_nnz);
	#endif
	if(ierr) _error_("PETSc could not allocate matrix (probably not enough memory)");
//	MatSetOption(this->matrix,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE);

}
/*}}}*/
PetscMat::PetscMat(IssmDouble* serial_mat,int M,int N,IssmDouble sparsity){/*{{{*/

	int     i;
	int* idxm=NULL;
	int* idxn=NULL;

	if(M)idxm=xNew<int>(M);
	if(N)idxn=xNew<int>(N);

	for(i=0;i<M;i++)idxm[i]=i;
	for(i=0;i<N;i++)idxn[i]=i;

	this->matrix=NewMat(M,N,sparsity,IssmComm::GetComm());
	MatSetValues(this->matrix,M,idxm,N,idxn,serial_mat,INSERT_VALUES);
	MatAssemblyBegin(this->matrix,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(this->matrix,MAT_FINAL_ASSEMBLY);

	xDelete<int>(idxm);
	xDelete<int>(idxn);

}
/*}}}*/
PetscMat::PetscMat(int M,int N, int connectivity,int numberofdofspernode){/*{{{*/

	this->matrix=NewMat(M,N,connectivity,numberofdofspernode,IssmComm::GetComm());

}
/*}}}*/
PetscMat::~PetscMat(){/*{{{*/
	MatFree(&this->matrix);
}
/*}}}*/

/*PetscMat specific routines: */
void PetscMat::AllocationInfo(void){/*{{{*/

	MatInfo info;
	MatGetInfo(this->matrix,MAT_GLOBAL_SUM,&info);
	_printf0_("=========================== Stiffness matrix allocation info ===========================\n");
	_printf0_("\n");
   _printf0_(" Block size  : "<<info.block_size << "\n");
	_printf0_(" nz_allocated: "<<info.nz_allocated << "\n");
	_printf0_(" nz_used     : "<<info.nz_used << "\n");
	_printf0_(" nz_unneeded : "<<info.nz_unneeded<<" ("<<double(info.nz_unneeded)/double(info.nz_allocated)*100.<<"%)\n");
	_printf0_("\n");
	_printf0_("========================================================================================\n");
}
/*}}}*/
void PetscMat::Echo(void){/*{{{*/

	MatView(this->matrix,PETSC_VIEWER_STDOUT_WORLD);
}
/*}}}*/
void PetscMat::Assemble(void){/*{{{*/

	_assert_(this->matrix);
	MatAssemblyBegin(this->matrix,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(this->matrix,MAT_FINAL_ASSEMBLY);

}
/*}}}*/
IssmDouble PetscMat::Norm(NormMode mode){/*{{{*/

	IssmDouble norm=0;
	_assert_(this->matrix);
	MatNorm(this->matrix,ISSMToPetscNormMode(mode),&norm);

	return norm;

}
/*}}}*/
void PetscMat::GetSize(int* pM,int* pN){/*{{{*/

	_assert_(this->matrix);
	MatGetSize(this->matrix,pM,pN);
}
/*}}}*/
void PetscMat::GetLocalSize(int* pM,int* pN){/*{{{*/

	_assert_(this->matrix);
	MatGetLocalSize(this->matrix,pM,pN);

}
/*}}}*/
void PetscMat::MatMult(PetscVec* X,PetscVec* AX){/*{{{*/

	_assert_(this->matrix);
	_assert_(X->vector);
	MatMultPatch(this->matrix,X->vector,AX->vector,IssmComm::GetComm());

}
/*}}}*/
PetscMat* PetscMat::Duplicate(void){/*{{{*/

	PetscMat* output=new PetscMat();
	_assert_(this->matrix);
	MatDuplicate(this->matrix,MAT_COPY_VALUES,&output->matrix);

	return output;

}
/*}}}*/
IssmDouble* PetscMat::ToMPISerial(void){/*{{{*/

	 IssmDouble* output=NULL;
	 MatToMPISerial(&output,this->matrix,IssmComm::GetComm(),true);
	 return output;

}
/*}}}*/
IssmDouble* PetscMat::ToMPISerial0(void){/*{{{*/

	 IssmDouble* output=NULL;
	 MatToMPISerial(&output,this->matrix,IssmComm::GetComm(),false);
	 return output;

}
/*}}}*/
void PetscMat::SetValues(int m,int* idxm,int n,int* idxn,IssmDouble* values,InsMode mode){/*{{{*/

	PetscErrorCode ierr = MatSetValues(this->matrix,m,idxm,n,idxn,values,ISSMToPetscInsertMode(mode));
	if(ierr) _error_("PETSc's MatSetValues reported an error");

}
/*}}}*/
void PetscMat::Convert(MatrixType type){/*{{{*/

	MatConvert(this->matrix,ISSMToPetscMatrixType(type),MAT_REUSE_MATRIX,&this->matrix);

}
/*}}}*/
void PetscMat::SetZero(void){/*{{{*/
	MatZeroEntries(this->matrix);
}
/*}}}*/

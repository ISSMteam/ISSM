/*!\file PetscMat.cpp
 * \brief: implementation of the PetscMat object
 */

/*Headers:*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif
#include <stdio.h>
#include <string.h>
#include "../petscincludes.h"
#include "../../../shared/shared.h"

#ifdef _HAVE_CODIPACK_
#include "../../codipack/CoDiPackDebug.h"
#endif

/*PetscMat constructors and destructor*/
template<typename doubletype>
PetscMat<doubletype>::PetscMat(){/*{{{*/
	this->matrix=NULL;
}
/*}}}*/
template<typename doubletype>
PetscMat<doubletype>::PetscMat(int M,int N){/*{{{*/

	this->matrix=NewMat(M,N,IssmComm::GetComm());
}
/*}}}*/
template<typename doubletype>
PetscMat<doubletype>::PetscMat(int M,int N, IssmPDouble sparsity){/*{{{*/

	this->matrix=NewMat(M,N,sparsity,IssmComm::GetComm());
}
/*}}}*/
template<typename doubletype>
PetscMat<doubletype>::PetscMat(int m,int n,int M,int N,int* d_nnz,int* o_nnz){/*{{{*/

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
template<typename doubletype>
PetscMat<doubletype>::PetscMat(doubletype* serial_mat,int M,int N,IssmPDouble sparsity){/*{{{*/

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
template<typename doubletype>
PetscMat<doubletype>::PetscMat(int M,int N, int connectivity,int numberofdofspernode){/*{{{*/

	this->matrix=NewMat(M,N,connectivity,numberofdofspernode,IssmComm::GetComm());

}
/*}}}*/
template<typename doubletype>
PetscMat<doubletype>::~PetscMat(){/*{{{*/
	MatFree(&this->matrix);
}
/*}}}*/

/*PetscMat specific routines: */
template<typename doubletype>
void PetscMat<doubletype>::AllocationInfo(void){/*{{{*/

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
template<typename doubletype>
void PetscMat<doubletype>::Echo(void){/*{{{*/

	MatView(this->matrix,PETSC_VIEWER_STDOUT_WORLD);
}
/*}}}*/
template<typename doubletype>
void PetscMat<doubletype>::EchoDebug(std::string message){/*{{{*/
#if defined(_HAVE_CODIPACK_) & defined(_HAVE_ADJOINTPETSC_)
	if (std::is_same<doubletype, IssmDouble>::value && CoDiIsDebugOutput()) {
		adjoint_petsc::ADMatDebugOutput(this->matrix, message, CoDiGetUniqueID());
	}
#endif
}
/*}}}*/
template<typename doubletype>
void PetscMat<doubletype>::Assemble(void){/*{{{*/

	_assert_(this->matrix);
	MatAssemblyBegin(this->matrix,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(this->matrix,MAT_FINAL_ASSEMBLY);

}
/*}}}*/
template<typename doubletype>
doubletype PetscMat<doubletype>::Norm(NormMode mode){/*{{{*/

	doubletype norm=0;
	_assert_(this->matrix);
	MatNorm(this->matrix,ISSMToPetscNormMode(mode),&norm);

	return norm;

}
/*}}}*/
template<typename doubletype>
void PetscMat<doubletype>::GetSize(int* pM,int* pN){/*{{{*/

	_assert_(this->matrix);
	MatGetSize(this->matrix,pM,pN);
}
/*}}}*/
template<typename doubletype>
void PetscMat<doubletype>::GetLocalSize(int* pM,int* pN){/*{{{*/

	_assert_(this->matrix);
	MatGetLocalSize(this->matrix,pM,pN);

}
/*}}}*/
template<typename doubletype>
void PetscMat<doubletype>::MatMult(PetscVec<doubletype>* X,PetscVec<doubletype>* AX){/*{{{*/

	_assert_(this->matrix);
	_assert_(X->vector);

  using ::MatMult;
#if _HAVE_CODIPACK_
  using ::adjoint_petsc::MatMult;
#endif

	MatMult(this->matrix, X->vector, AX->vector);
}/*}}}*/
template<typename doubletype>
PetscMat<doubletype>* PetscMat<doubletype>::Duplicate(void){/*{{{*/

	_assert_(this->matrix);

	/*Instantiate output Matrix*/
	PetscMat* output=new PetscMat();

	/*Duplicate matrix*/
	MatDuplicate(this->matrix,MAT_COPY_VALUES,&output->matrix);

	/*Return new matrix*/
	return output;
}
/*}}}*/
template<typename doubletype>
doubletype* PetscMat<doubletype>::ToMPISerial(void){/*{{{*/

	 doubletype* output=NULL;
	 MatToMPISerial(&output,this->matrix,IssmComm::GetComm(),true);
	 return output;

}
/*}}}*/
template<typename doubletype>
doubletype* PetscMat<doubletype>::ToMPISerial0(void){/*{{{*/

	 doubletype* output=NULL;
	 MatToMPISerial(&output,this->matrix,IssmComm::GetComm(),false);
	 return output;

}
/*}}}*/
template<typename doubletype>
void PetscMat<doubletype>::SetValues(int m,int* idxm,int n,int* idxn,doubletype* values,InsMode mode){/*{{{*/

	PetscErrorCode ierr = MatSetValues(this->matrix,m,idxm,n,idxn,values,ISSMToPetscInsertMode(mode));
	if(ierr) _error_("PETSc's MatSetValues reported an error");

}
/*}}}*/
template<typename doubletype>
void PetscMat<doubletype>::Convert(MatrixType type){/*{{{*/

	MatConvert(this->matrix,ISSMToPetscMatrixType(type),MAT_REUSE_MATRIX,&this->matrix);

}
/*}}}*/
template<typename doubletype>
void PetscMat<doubletype>::SetZero(void){/*{{{*/
	MatZeroEntries(this->matrix);
}
/*}}}*/

// Explicit instantiations.
template class PetscMat<IssmDouble>;

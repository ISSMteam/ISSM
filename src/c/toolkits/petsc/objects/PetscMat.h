/*!\file:  PetscMat.h
 * \brief wrapper to our own PetscMat object, which is needed to add AD capabilities (using ADOLC) 
 * to a C-coded Petsc API. We are just wrapping the Petsc objects into C++ equivalent, so that 
 * later, we can map all of the Petsc routines into Adolc equivalents.
 */ 

#ifndef _PETSCMAT_H_
#define _PETSCMAT_H_

/*Headers:*/
/*{{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "../petscincludes.h"
#include "../../../shared/Numerics/types.h"

/*}}}*/
class PetscVec;

class PetscMat{

	public:
		Mat matrix;

		#ifdef _HAVE_AD_
		IssmDouble* amatrix;
		#endif

		/*PetscMat constructors, destructors*/
		PetscMat();
		PetscMat(int M,int N);
		PetscMat(int M,int N,IssmDouble sparsity);
		PetscMat(int m,int n,int M,int N,int* d_nnz,int* o_nnz);
		PetscMat(IssmDouble* serial_mat,int M,int N,IssmDouble sparsity);
		PetscMat(int M,int N,int connectivity,int numberofdofspernode);
		~PetscMat();

		/*PetscMat specific routines*/
		void AllocationInfo(void);
		void Echo(void);
		void Assemble(void);
		IssmDouble Norm(NormMode norm_type);
		void GetSize(int* pM,int* pN);
		void GetLocalSize(int* pM,int* pN);
		void MatMult(PetscVec* X,PetscVec* AX);
		PetscMat* Duplicate(void);
		IssmDouble* ToMPISerial(void);
		IssmDouble* ToMPISerial0(void);
		void SetValues(int m,int* idxm,int n,int* idxn,IssmDouble* values,InsMode mode);
		void Convert(MatrixType type);
		void SetZero(void);
};

#endif //#ifndef _PETSCMAT_H_

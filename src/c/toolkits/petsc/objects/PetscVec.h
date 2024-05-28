/*!\file:  PetscVec.h
 * \brief wrapper to our own PetscVec object, which is needed to add AD capabilities (using ADOLC)
 * to a C-coded Petsc API. We are just wrapping the Petsc objects into C++ equivalent, so that
 * later, we can map all of the Petsc routines into Adolc equivalents.
 */

#ifndef _PETSCVEC_H_
#define _PETSCVEC_H_

/*Headers:*/
/*{{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "../petscincludes.h"

/*}}}*/

class PetscVec{

	public:
		Vec vector;

		#ifdef _HAVE_AD_
		IssmDouble* avector;
		#endif

		/*PetscVec constructors, destructors*/
		PetscVec();
		PetscVec(int M,bool fromlocalsize=false);
		PetscVec(int m,int M);
		PetscVec(IssmDouble* buffer, int M);
		PetscVec(Vec petsc_vec);
		~PetscVec();

		/*PetscVec specific routines*/
		void        Echo(void);
		void        Assemble(void);
		void        SetValues(int ssize, int* list, IssmDouble* values, InsMode mode);
		void        SetValue(int dof, IssmDouble value, InsMode  mode);
		void        GetValue(IssmDouble* pvalue, int dof);
		void        GetSize(int* pM);
		void        GetLocalSize(int* pM);
		void        GetLocalVector(IssmDouble** pvector,int** pindices);
		PetscVec*   Duplicate(void);
		void        Set(IssmDouble value);
		void        AXPY(PetscVec* X, IssmDouble a);
		void        AYPX(PetscVec* X, IssmDouble a);
		IssmDouble* ToMPISerial(void);
		IssmDouble* ToMPISerial0(void);
		void        Shift(IssmDouble shift);
		void        Copy(PetscVec* to);
		IssmDouble  Norm(NormMode norm_type);
		IssmDouble  Max(void);
		void        Scale(IssmDouble scale_factor);
		void        Pow(IssmDouble scale_factor);
		void        Sum(IssmDouble* pvalue);
		void        PointwiseDivide(PetscVec* x,PetscVec* y);
		void        PointwiseMult(PetscVec* x,PetscVec* y);
		IssmDouble  Dot(PetscVec* vector);
};

#endif //#ifndef _PETSCVEC_H_

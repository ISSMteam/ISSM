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

#if _HAVE_CODIPACK_
#include <adjoint_petsc/vec.h>
#endif

template<typename doubletype>
class PetscVec {

	public:
#if _HAVE_CODIPACK_
		using PVec = typename std::conditional<std::is_same<doubletype, IssmDouble>::value, adjoint_petsc::ADVec, Vec>::type;
#else
		using PVec = Vec;
#endif

		PVec vector;

		/*PetscVec constructors, destructors*/
		PetscVec();
		PetscVec(int M,bool fromlocalsize=false);
		PetscVec(int m,int M);
		PetscVec(doubletype* buffer, int M);
		PetscVec(PVec petsc_vec);
		~PetscVec();

		/*PetscVec specific routines*/
		void        Echo(void);
		void        EchoDebug(std::string message);
		void        Assemble(void);
		void        SetValues(int ssize, int* list, doubletype* values, InsMode mode);
		void        SetValue(int dof, doubletype value, InsMode  mode);
		void        GetValue(doubletype* pvalue, int dof);
		void        GetSize(int* pM);
		void        GetLocalSize(int* pM);
		void        GetLocalVector(doubletype** pvector,int** pindices);
		PetscVec*   Duplicate(void);
		void        Set(doubletype value);
		void        AXPY(PetscVec* X, doubletype a);
		void        AYPX(PetscVec* X, doubletype a);
		doubletype* ToMPISerial(void);
		doubletype* ToMPISerial0(void);
		void        Shift(doubletype shift);
		void        Copy(PetscVec* to);
		doubletype  Norm(NormMode norm_type);
		doubletype  Max(void);
		void        Scale(doubletype scale_factor);
		void        Pow(doubletype scale_factor);
		void        Sum(doubletype* pvalue);
		void        PointwiseDivide(PetscVec* x,PetscVec* y);
		void        PointwiseMult(PetscVec* x,PetscVec* y);
		doubletype  Dot(PetscVec* vector);

};

#endif //#ifndef _PETSCVEC_H_

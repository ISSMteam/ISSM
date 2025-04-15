/*!\file:  IssmAbsMat.h
 * \brief Main abstract class for the ISSM matrices.  This abstract class defines the pure virtual
 * functions that each of its descendants need to implement, such as contructors, destructors, as well
 * as matrix specific routines, such as SetValue, Assemple, MatMult, etc ...
 * Descendants include among others:
 * IssmDenseMat and IssmMpiDenseMat
 *
 */ 

#ifndef _ISSM_ABS_MAT_H_
#define _ISSM_ABS_MAT_H_

/*Headers:*/
#include "../toolkitsenums.h"
#include "../../shared/Numerics/types.h"

/*We need to template this class, in case we want to create Matrices that hold
  IssmDouble* matrix or IssmPDouble* matrix. 
  Such vectors are useful for use without or with the matlab or python
  interface (which do not care for IssmDouble types, but only rely on
  IssmPDouble types) */

template <class doubletype> class IssmAbsVec;
class Parameters;

template <class doubletype> 
class IssmAbsMat{

	public:

		/*IssmAbsMat constructors, destructors*/
		virtual ~IssmAbsMat(){};

		/*Functionality: */
		virtual void EchoDebug(std::string message) = 0;
		virtual void Echo(void)=0;
		virtual void Assemble(void)=0;
		virtual doubletype Norm(NormMode mode)=0;
		virtual void GetSize(int* pM,int* pN)=0;
		virtual void GetLocalSize(int* pM,int* pN)=0;
		virtual void MatMult(IssmAbsVec<doubletype>* X,IssmAbsVec<doubletype>* AX)=0;
		virtual IssmAbsMat<doubletype>* Duplicate(void)=0;
		virtual doubletype* ToSerial(void)=0;
		virtual void SetValues(int m,int* idxm,int n,int* idxn,doubletype* values,InsMode mode)=0;
		virtual void Convert(MatrixType type)=0;
		virtual void SetZero(void)=0;
		#ifndef _HAVE_WRAPPERS_
		virtual IssmAbsVec<IssmDouble>* Solve(IssmAbsVec<IssmDouble>* pf, Parameters* parameters)=0;
		#endif
};

#endif //#ifndef _ISSM_ABS_MAT_H_

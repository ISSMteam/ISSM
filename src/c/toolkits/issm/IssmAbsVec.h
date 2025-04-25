/*!\file:  IssmAbsVec.h
 * \brief Main abstract class for the ISSM vectors.  This abstract class defines the pure virtual
 * functions that each of its descendants need to implement, such as contructors, destructors, as well
 * as matrix specific routines, such as SetValue, Assemple, VecMult, etc ...
 * Descendants include among others:
 *	  IssmSeqVec and IssmMpiVec
 *
 */

#ifndef _ISSM_ABS_VEC_H_
#define _ISSM_ABS_VEC_H_

/*Headers:*/
#include "../toolkitsenums.h"
#include "../../shared/Numerics/types.h"

/*We need to template this class, in case we want to create Vectors that hold
  IssmDouble* vector or IssmPDouble* vector.
  Such vectors are useful for use without or with the matlab or python
  interface (which do not care for IssmDouble types, but only rely on
  IssmPDouble types)
*/
template <class doubletype>
class IssmAbsVec{

	public:

		/*IssmAbsVec constructors, destructors*/
		virtual ~IssmAbsVec(){};

		/*IssmAbsVec specific routines*/
		virtual void Echo(void)=0;
		virtual void EchoDebug(std::string message)=0;
		virtual void Assemble(void)=0;
		virtual void SetValues(int ssize, int* list, doubletype* values, InsMode mode)=0;
		virtual void SetValue(int dof, doubletype value, InsMode mode)=0;
		virtual void GetValue(doubletype* pvalue,int dof)=0;
		virtual void GetSize(int* pM)=0;
		virtual void GetLocalSize(int* pM)=0;
		virtual void GetLocalVector(doubletype** pvector,int** pindices)=0;
		virtual IssmAbsVec<doubletype>* Duplicate(void)=0;
		virtual void Set(doubletype value)=0;
		virtual void AXPY(IssmAbsVec* X, doubletype a)=0;
		virtual void AYPX(IssmAbsVec* X, doubletype a)=0;
		virtual doubletype* ToMPISerial(void)=0;
		virtual doubletype* ToMPISerial0(void)=0;
		virtual void Shift(doubletype shift)=0;
		virtual void Copy(IssmAbsVec* to)=0;
		virtual doubletype Norm(NormMode mode)=0;
		virtual void Scale(doubletype scale_factor)=0;
		virtual doubletype Dot(IssmAbsVec* input)=0;
		virtual void PointwiseDivide(IssmAbsVec* x,IssmAbsVec* y)=0;
		virtual void PointwiseMult(IssmAbsVec* x,IssmAbsVec* y)=0;
		virtual void Pow(doubletype scale_factor)=0;
		virtual void Sum(doubletype* pvalue)=0;
};

#endif //#ifndef _ISSM_ABS_VEC_H_

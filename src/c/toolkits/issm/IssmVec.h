/*!\file:  IssmVec.h
 * \brief Main Vector class for the Issm toolkit.
 */

#ifndef _ISSMVEC_H_
#define _ISSMVEC_H_

/*Headers:*/
/*{{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "../../shared/Enum/Enum.h"
#include "../../shared/Exceptions/exceptions.h"
#include "../../shared/MemOps/MemOps.h"
#include "./IssmToolkitUtils.h"
#include <math.h>
/*}}}*/

/*We need to template this class, in case we want to create Vectors that hold
  IssmDouble* matrix or IssmPDouble* matrix.
  Such vectors are useful for use without or with the matlab or python
  interface (which do not care for IssmDouble types, but only rely on
  IssmPDouble types)
*/
int IssmVecTypeFromToolkitOptions(void);
template <class doubletype> class IssmSeqVec;
template <class doubletype> class IssmMpiVec;

template <class doubletype>
class IssmVec{

	public:

		IssmAbsVec<doubletype>* vector; /*abstract vector, which implements object orientation*/

		/*IssmVec constructors, destructors*/
		IssmVec(){ /*{{{*/
			this->vector=NULL;
		}
		/*}}}*/
		IssmVec(int M){/*{{{*/

			switch(IssmVecTypeFromToolkitOptions()){

				case SeqEnum:
					this->vector=new IssmSeqVec<doubletype>(M);
					break;
				case MpiEnum:
					#ifdef _HAVE_MPI_
					this->vector=new IssmMpiVec<doubletype>(M);
					#else
					_error_("Mpi vector requires compilation of MPI!");
					#endif
					break;
				default:
					_error_("vector type not supported yet!");
			}
		}
		/*}}}*/
		IssmVec(int m,int M){/*{{{*/

			switch(IssmVecTypeFromToolkitOptions()){

				case SeqEnum:
					this->vector=new IssmSeqVec<doubletype>(m,M);
					break;
				case MpiEnum:
					#ifdef _HAVE_MPI_
					this->vector=new IssmMpiVec<doubletype>(m,M);
					#else
					_error_("Mpi vector requires compilation of MPI!");
					#endif
					break;
				default:
					_error_("vector type not supported yet!");
			}
		}
		/*}}}*/
		IssmVec(int M,bool fromlocalsize){/*{{{*/

			switch(IssmVecTypeFromToolkitOptions()){

				case SeqEnum:
					this->vector=new IssmSeqVec<doubletype>(M,fromlocalsize);
					break;
				case MpiEnum:
					#ifdef _HAVE_MPI_
					this->vector=new IssmMpiVec<doubletype>(M,fromlocalsize);
					#else
					_error_("Mpi vector requires compilation of MPI!");
					#endif
					break;
				default:
					_error_("vector type not supported yet!");
			}
		}
		/*}}}*/
		IssmVec(doubletype* buffer,int M){/*{{{*/

			switch(IssmVecTypeFromToolkitOptions()){

				case SeqEnum:
					this->vector=new IssmSeqVec<doubletype>(buffer,M);
					break;
				case MpiEnum:
					#ifdef _HAVE_MPI_
					this->vector=new IssmMpiVec<doubletype>(buffer,M);
					#else
					_error_("Mpi vector requires compilation of MPI!");
					#endif
					break;
				default:
					_error_("vector type not supported yet!");
			}
		}
		/*}}}*/
		~IssmVec(){/*{{{*/
			delete this->vector;
		}
		/*}}}*/

		/*IssmVec specific routines*/
		void Echo(void){/*{{{*/
			vector->Echo();
		}
		void EchoDebug(std::string message){/*{{{*/
			vector->EchoDebug(message);
		}
		/*}}}*/
		void Assemble(void){/*{{{*/
			vector->Assemble();
		}
		/*}}}*/
		void SetValues(int ssize, int* list, doubletype* values, InsMode mode){/*{{{*/
			vector->SetValues(ssize,list,values,mode);
		}
		/*}}}*/
		void SetValue(int dof, doubletype value, InsMode mode){/*{{{*/
			vector->SetValue(dof,value,mode);
		}
		/*}}}*/
		void GetValue(doubletype* pvalue,int dof){/*{{{*/
			vector->GetValue(pvalue,dof);
		}
		/*}}}*/
		void GetSize(int* pM){/*{{{*/
			vector->GetSize(pM);
		}
		/*}}}*/
		void GetLocalSize(int* pM){/*{{{*/
			vector->GetLocalSize(pM);
		}
		/*}}}*/
		void GetLocalVector(doubletype** pvector,int** pindices){/*{{{*/
			vector->GetLocalVector(pvector,pindices);
		} /*}}}*/
		IssmVec<doubletype>* Duplicate(void){/*{{{*/

			_assert_(this);
			IssmVec<doubletype>* issmvector=NULL;

			issmvector=new IssmVec<doubletype>();
			issmvector->vector=this->vector->Duplicate();
			this->vector->Copy(issmvector->vector);

			return issmvector;
		}
		/*}}}*/
		void Set(doubletype value){/*{{{*/
			vector->Set(value);
		}
		/*}}}*/
		void AXPY(IssmVec* X, doubletype a){/*{{{*/
			vector->AXPY(X->vector,a);
		}
		/*}}}*/
		void AYPX(IssmVec* X, doubletype a){/*{{{*/
			vector->AYPX(X->vector,a);
		}
		/*}}}*/
		doubletype* ToMPISerial(void){/*{{{*/
			return vector->ToMPISerial();
		}
		/*}}}*/
		doubletype* ToMPISerial0(void){/*{{{*/
			return vector->ToMPISerial0();
		}
		/*}}}*/
		void Shift(doubletype shift){/*{{{*/
			vector->Shift(shift);
		}
		/*}}}*/
		void Copy(IssmVec* to){/*{{{*/
			vector->Copy(to->vector);
		}
		/*}}}*/
		doubletype Norm(NormMode mode){/*{{{*/
			return vector->Norm(mode);
		}
		/*}}}*/
		void Scale(doubletype scale_factor){/*{{{*/
			vector->Scale(scale_factor);
		}
		/*}}}*/
		doubletype Dot(IssmVec* input){/*{{{*/
			return vector->Dot(input->vector);
		}
		/*}}}*/
		void PointwiseDivide(IssmVec* x,IssmVec* y){/*{{{*/
			vector->PointwiseDivide(x->vector,y->vector);
		}
		/*}}}*/
		void PointwiseMult(IssmVec* x,IssmVec* y){/*{{{*/
			vector->PointwiseMult(x->vector,y->vector);
		}
		/*}}}*/
		void Pow(doubletype scale_factor){/*{{{*/
			vector->Pow(scale_factor);
		}
		/*}}}*/
		void Sum(doubletype*  pvalue){/*{{{*/
			vector->Sum(pvalue);
		}
		/*}}}*/
};

#endif //#ifndef _ISSMVEC_H_

/*!\file:  Vector.h
 * \brief wrapper to vector objects. The goal is to control which API (PETSc,Scalpack, Plapack?)
 * implements our underlying vector format.
 */

#ifndef _VECTOR_H_
#define _VECTOR_H_

/*Headers:*/
/*{{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif
#include <cstring>
#include "../../shared/Enum/Enum.h"
#include "../petsc/petscincludes.h"
#include "../issm/issmtoolkit.h"
/*}}}*/

enum vectortype { PetscVecType, IssmVecType };

template <class doubletype>
class Vector{

	public:

		int  type;
		#ifdef _HAVE_PETSC_
		PetscVec<doubletype>* pvector;
		#endif
		IssmVec<doubletype>* ivector;

		/*Vector constructors, destructors */
		Vector(){ /*{{{*/

			InitCheckAndSetType();
		}
		/*}}}*/
		Vector(int M,bool fromlocalsize=false){ /*{{{*/

			InitCheckAndSetType();

			if(type==PetscVecType){
				#ifdef _HAVE_PETSC_
				this->pvector=new PetscVec<doubletype>(M,fromlocalsize);
				#endif
			}
			else this->ivector=new IssmVec<doubletype>(M,fromlocalsize);

		}
		/*}}}*/
		Vector(int m,int M){ /*{{{*/

			InitCheckAndSetType();

			if(type==PetscVecType){
				#ifdef _HAVE_PETSC_
					this->pvector=new PetscVec<doubletype>(m,M);
				 #endif
			}
			else this->ivector=new IssmVec<doubletype>(m,M);
		}
		/*}}}*/
		Vector(doubletype* serial_vec,int M){ /*{{{*/

			InitCheckAndSetType();

			if(type==PetscVecType){
				#ifdef _HAVE_PETSC_
				this->pvector=new PetscVec<doubletype>(serial_vec,M);
				#endif
			}
			else this->ivector=new IssmVec<doubletype>(serial_vec,M);
		}
		/*}}}*/
		~Vector(){ /*{{{*/

			if(type==PetscVecType){
				#ifdef _HAVE_PETSC_
				delete this->pvector;
				#endif
			}
			else delete this->ivector;
		}
		/*}}}*/
		#ifdef _HAVE_PETSC_
		Vector(PVec petsc_vector){ /*{{{*/

			this->type=PetscVecType;
			this->ivector=NULL;
			this->pvector=new PetscVec<doubletype>(petsc_vector);

		}
		/*}}}*/
		#endif
		void InitCheckAndSetType(void){ /*{{{*/

			#ifdef _HAVE_PETSC_
			pvector=NULL;
			#endif
			ivector=NULL;

			/*retrieve toolkittype: */
			char* toolkittype=ToolkitOptions::GetToolkitType();
			_assert_(toolkittype);

			/*set vector type: */
			if(strcmp(toolkittype,"petsc")==0){
				#ifdef _HAVE_PETSC_
				type=PetscVecType;
				#else
				_error_("cannot create petsc vector without PETSC compiled!");
				#endif
			}
			else if(strcmp(toolkittype,"issm")==0){
				/*let this choice stand:*/
				type=IssmVecType;
			}
			else{
				_error_("unknow toolkit type ");
			}

			/*Free resources: */
			xDelete<char>(toolkittype);
		}
		/*}}}*/

		/*Vector specific routines*/
		void Echo(void){_assert_(this);/*{{{*/

			if(type==PetscVecType){
				#ifdef _HAVE_PETSC_
				this->pvector->Echo();
				#endif
			}
			else this->ivector->Echo();

		}
		/*}}}*/
		void EchoDebug(std::string message){_assert_(this);/*{{{*/

			if(type==PetscVecType){
#ifdef _HAVE_PETSC_
				this->pvector->EchoDebug(message);
#endif
			}
			else this->ivector->EchoDebug(message);
		}
		/*}}}*/
		void Assemble(void){_assert_(this);/*{{{*/

			if(type==PetscVecType){
				#ifdef _HAVE_PETSC_
				this->pvector->Assemble();
				#endif
			}
			else this->ivector->Assemble();

		}
		/*}}}*/
		void SetValues(int ssize, int* list, doubletype* values, InsMode mode){ _assert_(this);/*{{{*/
			if(type==PetscVecType){
				#ifdef _HAVE_PETSC_
				this->pvector->SetValues(ssize,list,values,mode);
				#endif
			}
			else this->ivector->SetValues(ssize,list,values,mode);

		}
		/*}}}*/
		void SetValue(int dof, doubletype value, InsMode mode){_assert_(this);/*{{{*/

			if(type==PetscVecType){
				#ifdef _HAVE_PETSC_
				this->pvector->SetValue(dof,value,mode);
				#endif
			}
			else this->ivector->SetValue(dof,value,mode);

		}
		/*}}}*/
		void GetValue(doubletype* pvalue,int dof){_assert_(this);/*{{{*/

			if(type==PetscVecType){
				#ifdef _HAVE_PETSC_
				this->pvector->GetValue(pvalue,dof);
				#endif
			}
			else this->ivector->GetValue(pvalue,dof);

		}
		/*}}}*/
		void GetSize(int* pM){_assert_(this);/*{{{*/

			if(type==PetscVecType){
				#ifdef _HAVE_PETSC_
				this->pvector->GetSize(pM);
				#endif
			}
			else this->ivector->GetSize(pM);

		}
		/*}}}*/
		bool IsEmpty(void){/*{{{*/
			int M;

			_assert_(this);
			this->GetSize(&M);

			if(M==0)
				return true;
			else
				return false;
		}
		/*}}}*/
		void GetLocalSize(int* pM){_assert_(this);/*{{{*/

			if(type==PetscVecType){
				#ifdef _HAVE_PETSC_
				this->pvector->GetLocalSize(pM);
				#endif
			}
			else this->ivector->GetLocalSize(pM);

		}
		/*}}}*/
		void GetLocalVector(doubletype** pvector,int** pindices){_assert_(this);/*{{{*/

			if(type==PetscVecType){
				#ifdef _HAVE_PETSC_
				this->pvector->GetLocalVector(pvector,pindices);
				#endif
			}
			else this->ivector->GetLocalVector(pvector,pindices);

		}
		/*}}}*/
		Vector<doubletype>* Duplicate(void){_assert_(this);/*{{{*/

			Vector<doubletype>* output=NULL;

			output=new Vector<doubletype>();

			if(type==PetscVecType){
				#ifdef _HAVE_PETSC_
				output->pvector=this->pvector->Duplicate();
				#endif
			}
			else output->ivector=this->ivector->Duplicate();

			return output;
		} /*}}}*/
		void Set(doubletype value){_assert_(this);/*{{{*/

			if(type==PetscVecType){
				#ifdef _HAVE_PETSC_
				this->pvector->Set(value);
				#endif
			}
			else this->ivector->Set(value);

		}
		/*}}}*/
		void AXPY(Vector* X, doubletype a){_assert_(this);/*{{{*/

			if(type==PetscVecType){
				#ifdef _HAVE_PETSC_
				this->pvector->AXPY(X->pvector,a);
				#endif
			}
			else this->ivector->AXPY(X->ivector,a);

		}
		/*}}}*/
		void AYPX(Vector* X, doubletype a){_assert_(this);/*{{{*/

			if(type==PetscVecType){
				#ifdef _HAVE_PETSC_
				this->pvector->AYPX(X->pvector,a);
				#endif
			}
			else this->ivector->AYPX(X->ivector,a);
		}
		/*}}}*/
		doubletype* ToMPISerial(void){/*{{{*/

			doubletype* vec_serial=NULL;

			_assert_(this);
			if(type==PetscVecType){
				#ifdef _HAVE_PETSC_
				vec_serial=this->pvector->ToMPISerial();
				#endif
			}
			else vec_serial=this->ivector->ToMPISerial();

			return vec_serial;

		}
		/*}}}*/
		doubletype* ToMPISerial0(void){/*{{{*/

			doubletype* vec_serial=NULL;

			_assert_(this);
			if(type==PetscVecType){
				#ifdef _HAVE_PETSC_
				vec_serial=this->pvector->ToMPISerial0();
				#else
				_error_("Cannot serialize PETSc Vec without PETSc");
				#endif
			}
			else vec_serial=this->ivector->ToMPISerial0();

			return vec_serial;

		}
		/*}}}*/
		void Shift(doubletype shift){_assert_(this);/*{{{*/

			if(type==PetscVecType){
				#ifdef _HAVE_PETSC_
				this->pvector->Shift(shift);
				#endif
			}
			else this->ivector->Shift(shift);
		}
		/*}}}*/
		void Copy(Vector* to){_assert_(this);/*{{{*/

			if(type==PetscVecType){
				#ifdef _HAVE_PETSC_
				this->pvector->Copy(to->pvector);
				#endif
			}
			else this->ivector->Copy(to->ivector);
		}
		/*}}}*/
		doubletype Max(void){_assert_(this);/*{{{*/

			doubletype max=0;

			if(type==PetscVecType){
				#ifdef _HAVE_PETSC_
				max=this->pvector->Max();
				#endif
			}
			else _error_("operation not supported yet");
			return max;
		}
		/*}}}*/
		doubletype Norm(NormMode norm_type){_assert_(this);/*{{{*/

			doubletype norm=0;

			if(type==PetscVecType){
				#ifdef _HAVE_PETSC_
				norm=this->pvector->Norm(norm_type);
				#endif
			}
			else norm=this->ivector->Norm(norm_type);
			return norm;
		}
		/*}}}*/
		void Scale(doubletype scale_factor){_assert_(this);/*{{{*/

			if(type==PetscVecType){
				#ifdef _HAVE_PETSC_
				this->pvector->Scale(scale_factor);
				#endif
			}
			else this->ivector->Scale(scale_factor);
		}
		/*}}}*/
		doubletype Dot(Vector* vector){_assert_(this);/*{{{*/

			doubletype dot;

			if(type==PetscVecType){
				#ifdef _HAVE_PETSC_
				dot=this->pvector->Dot(vector->pvector);
				#endif
			}
			else dot=this->ivector->Dot(vector->ivector);
			return dot;
		}
		/*}}}*/
		void PointwiseDivide(Vector* x,Vector* y){_assert_(this);/*{{{*/

			if(type==PetscVecType){
				#ifdef _HAVE_PETSC_
				this->pvector->PointwiseDivide(x->pvector,y->pvector);
				#endif
			}
			else this->ivector->PointwiseDivide(x->ivector,y->ivector);
		}
		/*}}}*/
		void PointwiseMult(Vector* x,Vector* y){_assert_(this);/*{{{*/

			if(type==PetscVecType){
				#ifdef _HAVE_PETSC_
				this->pvector->PointwiseMult(x->pvector,y->pvector);
				#endif
			}
			else this->ivector->PointwiseMult(x->ivector,y->ivector);
		}
		/*}}}*/
		void Pow(doubletype scale_factor){_assert_(this);/*{{{*/

			if(type==PetscVecType){
				#ifdef _HAVE_PETSC_
				this->pvector->Pow(scale_factor);
				#endif
			}
			else this->ivector->Pow(scale_factor);
		}
		/*}}}*/
void Sum(doubletype* pvalue){ /*{{{*/
	_assert_(this);/*{{{*/

	if(type==PetscVecType){
		#ifdef _HAVE_PETSC_
		this->pvector->Sum(pvalue);
		#endif
	}
	else this->ivector->Sum(pvalue);
}
/*}}}*/
}; /*}}}*/
#endif //#ifndef _VECTOR_H_

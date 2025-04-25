/*!\file:  IssmSeqVec.h
 * \brief implementation of an ISSM vector which run serially (1 cpu only), which is made of a fully dense
 * vector. Internally, this dense vector is just a linear buffer of type doubletype.
 * This object needs to answer the API defined by the virtual functions in IssmAbsVec,
 * and the contructors required by IssmVec (see IssmVec.h)
 */

#ifndef _ISSM_SEQ_VEC_H_
#define _ISSM_SEQ_VEC_H_

/*Headers:*/
/*{{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "../../shared/shared.h"
#include <math.h>

/*}}}*/

/*We need to template this class, in case we want to create vectors that hold IssmDouble* vector or IssmPDouble* vector.
  Such vectors would be useful for use without or with the matlab or python interface (which do not care for IssmDouble types,
  but only rely on IssmPDouble types)*/

template <class doubletype> class IssmAbsVec;

template <class doubletype>
class IssmSeqVec: public IssmAbsVec<doubletype>{

	public:

		doubletype* vector;
		int M;

		/*IssmSeqVec constructors, destructors*/
		IssmSeqVec(){/*{{{*/

			this->M=0;
			this->vector=NULL;
		}
		/*}}}*/
		IssmSeqVec(int pM){/*{{{*/

			this->M=pM;
			this->vector=NULL;
			if(this->M) this->vector=xNewZeroInit<doubletype>(pM);
		}
		/*}}}*/
		IssmSeqVec(int pm,int pM){/*{{{*/

			this->M=pM;
			this->vector=NULL;
			if(this->M) this->vector=xNewZeroInit<doubletype>(pM);
		}
		/*}}}*/
		IssmSeqVec(int pM,bool fromlocalsize){/*{{{*/

			this->M=pM;
			this->vector=NULL;
			if(this->M) this->vector=xNewZeroInit<doubletype>(pM);
		}
		/*}}}*/
		IssmSeqVec(doubletype* buffer,int pM){/*{{{*/

			this->M=pM;
			this->vector=NULL;
			if(this->M){
				this->vector=xNew<doubletype>(pM);
				xMemCpy<doubletype>(this->vector,buffer,pM);
			}
		}
		/*}}}*/
		~IssmSeqVec(){/*{{{*/
			if(this->M)xDelete<doubletype>(this->vector);
			M=0;
		}
		/*}}}*/

		/*IssmSeqVec specific routines*/
		void Echo(void){/*{{{*/

			int i;
			_printf_("IssmSeqVec size " << this->M << "\n");
			for(i=0;i<M;i++){
				_printf_(vector[i] << "\n ");
			}
		}
		/*}}}*/
		void EchoDebug(std::string message){/*{{{*/
			_printf_("Error: EchoeReverse SeqVec\n");
		}
		void Assemble(void){/*{{{*/

			/*do nothing*/

		}
		/*}}}*/
		void SetValues(int ssize, int* list, doubletype* values, InsMode mode){/*{{{*/

			int i;
			switch(mode){
				case ADD_VAL:
					for(i=0;i<ssize;i++) if(list[i]>=0) this->vector[list[i]]+=values[i];
					break;
				case INS_VAL:
					for(i=0;i<ssize;i++) if(list[i]>=0) this->vector[list[i]]=values[i];
					break;
				default:
					_error_("unknown insert mode!");
					break;
			}

		}
		/*}}}*/
		void SetValue(int dof, doubletype value,InsMode mode){/*{{{*/

			_assert_(dof>=0);

			switch(mode){
				case ADD_VAL:
					this->vector[dof]+=value;
					break;
				case INS_VAL:
					this->vector[dof]=value;
					break;
				default:
					_error_("unknown insert mode!");
					break;
			}
		}
		/*}}}*/
		void GetValue(doubletype* pvalue,int dof){/*{{{*/

			*pvalue=this->vector[dof];

		}
		/*}}}*/
		void GetSize(int* pM){/*{{{*/

			*pM=this->M;

		}
		/*}}}*/
		void GetLocalSize(int* pM){/*{{{*/
			*pM=this->M;
		}/*}}}*/
		void GetLocalVector(doubletype** pvector,int** pindices){/*{{{*/

			/*First, check that vector size is not 0*/
			int vector_size;
			this->GetSize(&vector_size);
			if(vector_size==0){
				*pvector=NULL;
				*pindices=NULL;
				return;
			}

			/*Build indices*/
			int* indices=xNew<int>(vector_size);
			for(int i=0;i<vector_size;i++) indices[i]=i;

			/*Get vector*/
			doubletype* values =xNew<doubletype>(vector_size);
			xMemCpy<doubletype>(values,this->vector,vector_size);

			*pvector  = values;
			*pindices = indices;
		} /*}}}*/
		IssmSeqVec<doubletype>* Duplicate(void){/*{{{*/

			return new IssmSeqVec<doubletype>(this->vector,this->M);

		}
		/*}}}*/
		void Set(doubletype value){/*{{{*/

			for(int i=0;i<this->M;i++)this->vector[i]=value;

		}
		/*}}}*/
		void AXPY(IssmAbsVec<doubletype>* Xin, doubletype a){/*{{{*/

			int i;

			/*Assume X is of the correct type, and downcast: */
			IssmSeqVec* X=NULL;

			X=(IssmSeqVec<doubletype>*)Xin;

			/*y=a*x+y where this->vector is y*/
			for(i=0;i<this->M;i++)this->vector[i]=a*X->vector[i]+this->vector[i];

		}
		/*}}}*/
		void AYPX(IssmAbsVec<doubletype>* Xin, doubletype a){/*{{{*/

			int i;

			/*Assume X is of the correct type, and downcast: */
			IssmSeqVec* X=NULL;

			X=(IssmSeqVec<doubletype>*)Xin;

			/*y=x+a*y where this->vector is y*/
			for(i=0;i<this->M;i++)this->vector[i]=X->vector[i]+a*this->vector[i];

		}
		/*}}}*/
		doubletype* ToMPISerial(void){/*{{{*/

			doubletype* buffer=NULL;

			if(this->M){
				buffer=xNew<doubletype>(this->M);
				xMemCpy<doubletype>(buffer,this->vector,this->M);
			}
			return buffer;

		}
		/*}}}*/
		doubletype* ToMPISerial0(void){/*{{{*/

			return this->ToMPISerial();

		}
		/*}}}*/
		void Shift(doubletype shift){/*{{{*/

			for(int i=0;i<this->M;i++)this->vector[i]+=shift;

		}
		/*}}}*/
		void Copy(IssmAbsVec<doubletype>* toin){/*{{{*/

			int i;

			/*Assume toin is of the correct type, and downcast: */
			IssmSeqVec* to=NULL;

			to=(IssmSeqVec<doubletype>*)toin;

			to->M=this->M;
			for(i=0;i<this->M;i++)to->vector[i]=this->vector[i];

		}
		/*}}}*/
		doubletype Norm(NormMode mode){/*{{{*/

			doubletype norm;
			int i;

			switch(mode){
				case NORM_INF:
					norm=0.; for(i=0;i<this->M;i++)norm=max(norm,fabs(this->vector[i]));
					//norm=0.; for(i=0;i<this->M;i++)norm=max(norm,this->vector[i]);
					return norm;
					break;
				case NORM_TWO:
					norm=0.;
					for(i=0;i<this->M;i++)norm+=this->vector[i]*this->vector[i];
					return sqrt(norm);
					break;
				default:
					_error_("unknown norm !");
					break;
			}

			return 0.;
		}
		/*}}}*/
		void Scale(doubletype scale_factor){/*{{{*/

			int i;
			for(i=0;i<this->M;i++)this->vector[i]=scale_factor*this->vector[i];

		}
		/*}}}*/
		doubletype Dot(IssmAbsVec<doubletype>* inputin){/*{{{*/

			int i;

			/*Assume inputin is of the correct type, and downcast: */
			IssmSeqVec* input=NULL;

			input=(IssmSeqVec<doubletype>*)inputin;

			doubletype dot=0;
			for(i=0;i<this->M;i++)dot+=this->vector[i]*input->vector[i];
			return dot;

		}
		/*}}}*/
		void PointwiseDivide(IssmAbsVec<doubletype>* xin,IssmAbsVec<doubletype>* yin){/*{{{*/

			int i;

			/*Assume xin and yin are of the correct type, and downcast: */
			IssmSeqVec* x=NULL;
			IssmSeqVec* y=NULL;

			x=(IssmSeqVec<doubletype>*)xin;
			y=(IssmSeqVec<doubletype>*)yin;

			/*pointwise w=x/y where this->vector is w: */
			for(i=0;i<this->M;i++)this->vector[i]=x->vector[i]/y->vector[i];
		}
		/*}}}*/
		void PointwiseMult(IssmAbsVec<doubletype>* xin,IssmAbsVec<doubletype>* yin){/*{{{*/

			int i;

			/*Assume xin and yin are of the correct type, and downcast: */
			IssmSeqVec* x=NULL;
			IssmSeqVec* y=NULL;

			x=(IssmSeqVec<doubletype>*)xin;
			y=(IssmSeqVec<doubletype>*)yin;

			/*pointwise w=x*y where this->vector is w: */
			for(i=0;i<this->M;i++)this->vector[i]=x->vector[i]*y->vector[i];
		}
		/*}}}*/
		void Pow(doubletype scale_factor){/*{{{*/

			int i;
			for(i=0;i<this->M;i++)this->vector[i]=pow(this->vector[i],scale_factor);

		}
		/*}}}*/
		void Sum(doubletype* pvalue){/*{{{*/

			doubletype value=0;
			int i;
			for(i=0;i<this->M;i++)value+=this->vector[i];
			*pvalue=value;
		}
		/*}}}*/

};
#endif //#ifndef _ISSM_SEQ_VEC_H_

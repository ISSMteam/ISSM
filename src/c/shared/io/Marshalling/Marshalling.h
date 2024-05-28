/*\file Marshalling.h
 *\brief: macros to help automate the marshalling, demarshalling, and marshalling size routines. 
 */

#ifndef _MARSHALLING_H_
#define _MARSHALLING_H_

#include <string.h>
#include "../../Exceptions/exceptions.h"
#include "../../MemOps/MemOps.h"
#include "../../Numerics/recast.h"

/*Define Marshall operation Enums first*/
enum MarshallOpEnum{
	MARSHALLING_WRITE,
	MARSHALLING_LOAD,
	MARSHALLING_SIZE,
#if defined(_HAVE_CODIPACK_) && !defined(_WRAPPERS_)
	AD_COUNTDOUBLES,
	AD_REGISTERINPUT,
	AD_REGISTEROUTPUT,
	AD_SETADJOINT,
#endif
};

/*Define virtual Marshall Handle*/
class MarshallHandle{ /*{{{*/
	public:
		MarshallOpEnum operation_enum;
		MarshallHandle(MarshallOpEnum op_in) : operation_enum(op_in){}
		~MarshallHandle(){}
		virtual void Echo(void)=0;
		int OperationNumber(){return operation_enum;}
		template<typename T> void call(T  & value);
		template<typename T> void call(T* & value,int size);
}; /*}}}*/
/* !! Make sure to pass all fields by reference !! */
class WriteCheckpointFunctor: public MarshallHandle{ /*{{{*/

	private:
		char** pmarshalled_data;

	public:
		WriteCheckpointFunctor(char** pmarshalled_data_in);
		void Echo(void);
		template<typename T> void call(T & value){
			memcpy(*pmarshalled_data,&value,sizeof(T));
			*pmarshalled_data+=sizeof(T);
		}
		void call(char* & value);
		template<typename T> void call(T* & value,int size){
			bool pointer_null = true;
			if(value) pointer_null = false;
			this->call<bool>(pointer_null);
			if(value){
				memcpy(*pmarshalled_data,value,size*sizeof(T));
				*pmarshalled_data+=size*sizeof(T);
			}
		}
};/*}}}*/
class LoadCheckpointFunctor:  public MarshallHandle{ /*{{{*/

	private:
		char** pmarshalled_data;

	public:
		LoadCheckpointFunctor(char** pmarshalled_data_in);
		void Echo(void);
		void call(char* & value);
		template<typename T> void call(T & value){
			memcpy(&value,*pmarshalled_data,sizeof(T));
			*pmarshalled_data+=sizeof(T);
		}
		template<typename T> void call(T* & value,int size){
			bool pointer_null;
			call(pointer_null);
			if(!pointer_null){
				value=xNew<T>(size);
				memcpy(value,*pmarshalled_data,size*sizeof(T));
				*pmarshalled_data+=size*sizeof(T);
			}
			else{
				value = NULL;
			}
		}
};/*}}}*/
class SizeCheckpointFunctor:  public MarshallHandle{ /*{{{*/

	private:
		int marshalled_data_size;

	public:
		SizeCheckpointFunctor(void);
		int MarshalledSize(void);
		void Echo(void);
		template<typename T> void call(T & value){
			marshalled_data_size+=sizeof(T);
		}
		void call(char* & value);
		template<typename T> void call(T* & value,int size){
			bool pointer_null = true;
			if(value) pointer_null = false;
			this->call(pointer_null);
			if(!pointer_null){
				marshalled_data_size+=size*sizeof(T);
			}
		}
};/*}}}*/
#if defined(_HAVE_CODIPACK_) && !defined(_WRAPPERS_)
#if _CODIPACK_MAJOR_==2
using Tape = typename IssmDouble::Tape;
#elif _CODIPACK_MAJOR_==1
using Tape = typename IssmDouble::TapeType;
#else
#error "_CODIPACK_MAJOR_ not supported"
#endif
class CountDoublesFunctor:    public MarshallHandle{ /*{{{*/

	private:
		int double_count;

	public:
		CountDoublesFunctor(void);
		int  DoubleCount(void);
		void Echo(void);
		template<typename T> void call(T & value){/*General case: do nothing*/}
		template<typename T> void call(T* & value,int size){/*General case: do nothing*/}
		void call(IssmDouble & value);
		void call(IssmDouble* & value,int size);
}; /*}}}*/
class RegisterInputFunctor:   public MarshallHandle{ /*{{{*/

	private:
		int  double_count;
		int *identifiers;
		int  size_max;
		Tape *tape_codi;

	public:
		RegisterInputFunctor(int* identifiers_in,int size_max_in);
		void Echo(void);
		template<typename T> void call(T & value){/*General case: do nothing*/}
		template<typename T> void call(T* & value,int size){/*General case: do nothing*/}
		void call(IssmDouble & value);
		void call(IssmDouble* & value,int size);
}; /*}}}*/
class RegisterOutputFunctor:  public MarshallHandle{ /*{{{*/

	private:
		int   double_count;
		Tape *tape_codi;

	public:
		RegisterOutputFunctor(void);
		void Echo(void);
		template<typename T> void call(T & value){/*General case: do nothing*/}
		template<typename T> void call(T* & value,int size){/*General case: do nothing*/}
		void call(IssmDouble & value);
		void call(IssmDouble* & value,int size);
}; /*}}}*/
class SetAdjointFunctor:      public MarshallHandle{ /*{{{*/

	private:
		int     double_count;
		int     size_max;
		Tape   *tape_codi;
		double *adjoint;

	public:
		SetAdjointFunctor(double* adjoint_in,int size_max_in);
		void Echo(void);
		template<typename T> void call(T & value){/*General case: do nothing*/}
		template<typename T> void call(T* & value,int size){/*General case: do nothing*/}
		void call(IssmDouble & value);
		void call(IssmDouble* & value,int size);
}; /*}}}*/
#endif

template<typename T> void MarshallHandle::call(T & value){
	switch(OperationNumber()){
		case MARSHALLING_WRITE:{WriteCheckpointFunctor* temp = xDynamicCast<WriteCheckpointFunctor*>(this); temp->call(value); break;}
		case MARSHALLING_LOAD: {LoadCheckpointFunctor*  temp = xDynamicCast<LoadCheckpointFunctor*>(this);  temp->call(value); break;}
		case MARSHALLING_SIZE: {SizeCheckpointFunctor*  temp = xDynamicCast<SizeCheckpointFunctor*>(this);  temp->call(value); break;}
#if defined(_HAVE_CODIPACK_) && !defined(_WRAPPERS_)
		case AD_COUNTDOUBLES:  {CountDoublesFunctor*   temp = xDynamicCast<CountDoublesFunctor*>(this);    temp->call(value); break;}
		case AD_REGISTERINPUT: {RegisterInputFunctor*  temp = xDynamicCast<RegisterInputFunctor*>(this);   temp->call(value); break;}
		case AD_REGISTEROUTPUT:{RegisterOutputFunctor* temp = xDynamicCast<RegisterOutputFunctor*>(this);  temp->call(value); break;}
		case AD_SETADJOINT:    {SetAdjointFunctor*     temp = xDynamicCast<SetAdjointFunctor*>(this);      temp->call(value); break;}
#endif
		default: _error_("Operation "<<OperationNumber()<<" not supported yet");
	}
}
template<typename T> void MarshallHandle::call(T* & value,int size){
	switch(OperationNumber()){
		case MARSHALLING_WRITE:{WriteCheckpointFunctor* temp = xDynamicCast<WriteCheckpointFunctor*>(this); temp->call(value,size); break;}
		case MARSHALLING_LOAD: {LoadCheckpointFunctor*  temp = xDynamicCast<LoadCheckpointFunctor*>(this);  temp->call(value,size); break;}
		case MARSHALLING_SIZE: {SizeCheckpointFunctor*  temp = xDynamicCast<SizeCheckpointFunctor*>(this);  temp->call(value,size); break;}
#if defined(_HAVE_CODIPACK_) && !defined(_WRAPPERS_)
		case AD_COUNTDOUBLES:  {CountDoublesFunctor*   temp = xDynamicCast<CountDoublesFunctor*>(this);    temp->call(value,size); break;}
		case AD_REGISTERINPUT: {RegisterInputFunctor*  temp = xDynamicCast<RegisterInputFunctor*>(this);   temp->call(value,size); break;}
		case AD_REGISTEROUTPUT:{RegisterOutputFunctor* temp = xDynamicCast<RegisterOutputFunctor*>(this);  temp->call(value,size); break;}
		case AD_SETADJOINT:    {SetAdjointFunctor*     temp = xDynamicCast<SetAdjointFunctor*>(this);      temp->call(value,size); break;}
#endif
		default: _error_("Operation "<<OperationNumber() <<" not supported yet");
	}
}

#endif

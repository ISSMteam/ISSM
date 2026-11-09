/*!\file IntArrayInput.c
 * \brief: implementation of the IntArrayInput object
 */

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "../classes.h"
#include "../../shared/shared.h"
#include "./IntArrayInput.h"

/*IntArrayInput constructors and destructor*/
template<typename T> TIntArrayInput<T>::TIntArrayInput(void){/*{{{*/

	this->numberofelements_local = -1;
	this->N                      = NULL;
	this->values                 = NULL;

}/*}}}*/
template<typename T> TIntArrayInput<T>::TIntArrayInput(int nbe_in){/*{{{*/

	_assert_(nbe_in>0);
	_assert_(nbe_in<1e11);
	this->numberofelements_local = nbe_in;
	this->N                      = xNewZeroInit<int>(this->numberofelements_local);
	this->values                 = xNewZeroInit<T*>(this->numberofelements_local);

}/*}}}*/
template<typename T> TIntArrayInput<T>::~TIntArrayInput(){/*{{{*/
	if(this->values){
		for(int i=0;i<this->numberofelements_local;i++) if(this->values[i]) xDelete<T>(this->values[i]);
		xDelete<T*>(this->values);
	}
	if(this->N) xDelete<int>(this->N);
}
/*}}}*/

/*Object virtual functions definitions:*/
template<typename T> int TIntArrayInput<T>::ObjectEnum(){/*{{{*/
	_error_("TIntArrayInput does not support this integer type");
	return -1;
}/*}}}*/
template<> int TIntArrayInput<int>::ObjectEnum(){return IntArrayInputEnum;}
template<> int TIntArrayInput<int16_t>::ObjectEnum(){return Int16ArrayInputEnum;}

template<typename T> Input* TIntArrayInput<T>::copy() {/*{{{*/

	TIntArrayInput<T>* output = new TIntArrayInput<T>(this->numberofelements_local);

	xMemCpy<int>(output->N,this->N,this->numberofelements_local);

	for(int i=0;i<this->numberofelements_local;i++){
		if(this->values[i]){
			_assert_(this->N[i]>0);
			output->values[i] = xNew<T>(this->N[i]);
			xMemCpy<T>(output->values[i],this->values[i],this->N[i]);
		}
		else{
			output->values[i] = NULL;
		}
	}

	return output;
}
/*}}}*/
template<typename T> void TIntArrayInput<T>::DeepEcho(void){/*{{{*/
	_printf_("IntArrayInput Echo:\n");
	///_printf_("   Size:          "<<N<<"\n");
	//printarray(this->values,this->M,this->N);
	//_printf_(setw(15)<<"   IntArrayInput "<<setw(25)<<left<<EnumToStringx(this->enum_type)<<" "<<(value?"true":"false") << "\n");
}
/*}}}*/
template<typename T> void TIntArrayInput<T>::Echo(void){/*{{{*/
	this->DeepEcho();
}
/*}}}*/
template<typename T> void TIntArrayInput<T>::Marshall(MarshallHandle* marshallhandle){ /*{{{*/

	int object_enum = this->ObjectEnum();
	marshallhandle->call(object_enum);
	marshallhandle->call(this->numberofelements_local);

	/*Allocate memory if reading restart file*/
	if(marshallhandle->OperationNumber() == MARSHALLING_LOAD){
		_assert_(this->numberofelements_local>0);
		_assert_(this->numberofelements_local<1e11);
		if(this->numberofelements_local){
			this->N                      = xNewZeroInit<int>(this->numberofelements_local);
			this->values                 = xNewZeroInit<T*>(this->numberofelements_local);
		}
		else{
			this->N      = NULL;
			this->values = NULL;
		}
	}
	/*Marshall N*/
	if(this->numberofelements_local){
		marshallhandle->call(this->N,this->numberofelements_local);
	}

	/*Marshall individual arrays*/
	if(this->numberofelements_local){

		for(int i=0;i<this->numberofelements_local;i++){
			if(this->N[i]){

				/*Allocate if reading restart*/
				if(marshallhandle->OperationNumber() == MARSHALLING_LOAD){
					this->values[i] = xNew<T>(this->N[i]);
				}
				_assert_(this->values[i]);
				marshallhandle->call(this->values[i],this->N[i]);
			}
		}
	}

}
/*}}}*/

/*IntArrayInput management*/
template<typename T> void TIntArrayInput<T>::SetInput(int row,int numindices,T* values_in){/*{{{*/

	_assert_(this);
	_assert_(row>=0 && row<this->numberofelements_local);

	if(this->N[row] != numindices){
		if(this->values[row]) xDelete<T>(this->values[row]);
		this->values[row] = xNew<T>(numindices);
	}

	T *el_values = this->values[row];
	for(int i=0;i<numindices;i++) el_values[i] = values_in[i];

	this->N[row] = numindices;
}
/*}}}*/
template<typename T> void TIntArrayInput<T>::GetArray(int row,T** pvalues,int* pN){/*{{{*/

	_assert_(this);
	_assert_(row>=0 && row<this->numberofelements_local);
	if(pvalues){
		T* outvalues = xNew<T>(this->N[row]);
		xMemCpy<T>(outvalues,this->values[row],this->N[row]);
		*pvalues = outvalues;
	}
	if(pN){
		*pN = this->N[row];
	}
}
/*}}}*/
template<typename T> void TIntArrayInput<T>::GetArrayPtr(int row,T** pvalues,int* pN){/*{{{*/

	_assert_(this);
	_assert_(row>=0 && row<this->numberofelements_local);
	if(pvalues){
		*pvalues = this->values[row];
	}
	if(pN){
		*pN = this->N[row];
	}
}
/*}}}*/

/*Explicit instantiations (add new integer types here)*/
template class TIntArrayInput<int>;
template class TIntArrayInput<int16_t>;

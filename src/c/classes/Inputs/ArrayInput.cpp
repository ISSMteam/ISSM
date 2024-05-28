/*!\file ArrayInput.c
 * \brief: implementation of the ArrayInput object
 */

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "../classes.h"
#include "../../shared/shared.h"
#include "./ArrayInput.h"

/*ArrayInput constructors and destructor*/
ArrayInput::ArrayInput(void){/*{{{*/

	this->numberofelements_local = -1;
	this->N                      = NULL;
	this->values                 = NULL;

}/*}}}*/
ArrayInput::ArrayInput(int nbe_in){/*{{{*/

	_assert_(nbe_in>0);
	_assert_(nbe_in<1e11);
	this->numberofelements_local = nbe_in;
	this->N                      = xNewZeroInit<int>(this->numberofelements_local);
	this->values                 = xNewZeroInit<IssmDouble*>(this->numberofelements_local);

}/*}}}*/
ArrayInput::~ArrayInput(){/*{{{*/
	if(this->values){
		for(int i=0;i<this->numberofelements_local;i++) if(this->values[i]) xDelete<IssmDouble>(this->values[i]);
		xDelete<IssmDouble>(this->values);
	}
	if(this->N) xDelete<int>(this->N);
}
/*}}}*/

/*Object virtual functions definitions:*/
Input* ArrayInput::copy() {/*{{{*/

	ArrayInput* output = new ArrayInput(this->numberofelements_local);

	xMemCpy<int>(output->N,this->N,this->numberofelements_local);

	for(int i=0;i<this->numberofelements_local;i++){
		if(this->values[i]){
			_assert_(this->N[i]>0);
			output->values[i] = xNew<IssmDouble>(this->N[i]);
			xMemCpy<IssmDouble>(output->values[i],this->values[i],this->N[i]);
		}
		else{
			output->values[i] = NULL;
		}
	}

	return output;
}
/*}}}*/
void ArrayInput::DeepEcho(void){/*{{{*/
	_printf_("ArrayInput Echo:\n");
	///_printf_("   Size:          "<<N<<"\n");
	//printarray(this->values,this->M,this->N);
	//_printf_(setw(15)<<"   ArrayInput "<<setw(25)<<left<<EnumToStringx(this->enum_type)<<" "<<(value?"true":"false") << "\n");
}
/*}}}*/
void ArrayInput::Echo(void){/*{{{*/
	this->DeepEcho();
}
/*}}}*/
int  ArrayInput::Id(void){/*{{{*/
	return -1;
}/*}}}*/
void ArrayInput::Marshall(MarshallHandle* marshallhandle){ /*{{{*/

	int object_enum = ArrayInputEnum;
	marshallhandle->call(object_enum);
	marshallhandle->call(this->numberofelements_local);

	/*Allocate memory if reading restart file*/
	if(marshallhandle->OperationNumber() == MARSHALLING_LOAD){
		_assert_(this->numberofelements_local>0);
		_assert_(this->numberofelements_local<1e11);
		if(this->numberofelements_local){
			this->N                      = xNewZeroInit<int>(this->numberofelements_local);
			this->values                 = xNewZeroInit<IssmDouble*>(this->numberofelements_local);
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
					this->values[i] = xNew<IssmDouble>(this->N[i]);
				}
				_assert_(this->values[i]);
				marshallhandle->call(this->values[i],this->N[i]);
			}
		}
	}
}
/*}}}*/
int  ArrayInput::ObjectEnum(void){/*{{{*/
	return ArrayInputEnum;
}
/*}}}*/

/*ArrayInput management*/
void ArrayInput::SetInput(int row,int numindices,IssmDouble* values_in){/*{{{*/

	_assert_(this);
	_assert_(row>=0 && row<this->numberofelements_local);

	if(this->N[row] != numindices){
		if(this->values[row]) xDelete<IssmDouble>(this->values[row]);
		this->values[row] = xNew<IssmDouble>(numindices);
	}

	IssmDouble *el_values = this->values[row];
	for(int i=0;i<numindices;i++) el_values[i] = values_in[i];

	this->N[row] = numindices;
}
/*}}}*/
void ArrayInput::GetArray(int row,IssmDouble** pvalues,int* pN){/*{{{*/

	_assert_(this);
	_assert_(row>=0 && row<this->numberofelements_local);
	if(pvalues){
		IssmDouble* outvalues = xNew<IssmDouble>(this->N[row]);
		xMemCpy<IssmDouble>(outvalues,this->values[row],this->N[row]);
		*pvalues = outvalues;
	}
	if(pN){
		*pN = this->N[row];
	}
}
/*}}}*/
void ArrayInput::GetArrayPtr(int row,IssmDouble** pvalues,int* pN){/*{{{*/

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

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
IntArrayInput::IntArrayInput(void){/*{{{*/

	this->numberofelements_local = -1;
	this->N                      = NULL;
	this->values                 = NULL;

}/*}}}*/
IntArrayInput::IntArrayInput(int nbe_in){/*{{{*/

	_assert_(nbe_in>0);
	_assert_(nbe_in<1e11);
	this->numberofelements_local = nbe_in;
	this->N                      = xNewZeroInit<int>(this->numberofelements_local);
	this->values                 = xNewZeroInit<int*>(this->numberofelements_local);

}/*}}}*/
IntArrayInput::~IntArrayInput(){/*{{{*/
	if(this->values){
		for(int i=0;i<this->numberofelements_local;i++) if(this->values[i]) xDelete<int>(this->values[i]);
		xDelete<int>(this->values);
	}
	if(this->N) xDelete<int>(this->N);
}
/*}}}*/

/*Object virtual functions definitions:*/
Input* IntArrayInput::copy() {/*{{{*/

	IntArrayInput* output = new IntArrayInput(this->numberofelements_local);

	xMemCpy<int>(output->N,this->N,this->numberofelements_local);

	for(int i=0;i<this->numberofelements_local;i++){
		if(this->values[i]){
			_assert_(this->N[i]>0);
			output->values[i] = xNew<int>(this->N[i]);
			xMemCpy<int>(output->values[i],this->values[i],this->N[i]);
		}
		else{
			output->values[i] = NULL;
		}
	}

	return output;
}
/*}}}*/
void IntArrayInput::DeepEcho(void){/*{{{*/
	_printf_("IntArrayInput Echo:\n");
	///_printf_("   Size:          "<<N<<"\n");
	//printarray(this->values,this->M,this->N);
	//_printf_(setw(15)<<"   IntArrayInput "<<setw(25)<<left<<EnumToStringx(this->enum_type)<<" "<<(value?"true":"false") << "\n");
}
/*}}}*/
void IntArrayInput::Echo(void){/*{{{*/
	this->DeepEcho();
}
/*}}}*/
int  IntArrayInput::Id(void){/*{{{*/
	return -1;
}/*}}}*/
void IntArrayInput::Marshall(MarshallHandle* marshallhandle){ /*{{{*/

	int object_enum = IntArrayInputEnum;
	marshallhandle->call(object_enum);
	marshallhandle->call(this->numberofelements_local);
	if(this->numberofelements_local){
		marshallhandle->call(this->N,this->numberofelements_local);
		for(int i=0;i<this->numberofelements_local;i++){
			if(this->values[i]){
				marshallhandle->call(this->values[i],this->N[i]);
			}
		}
	}
	else{
		this->N      = NULL;
		this->values = NULL;
	}

}
/*}}}*/
int  IntArrayInput::ObjectEnum(void){/*{{{*/
	return IntArrayInputEnum;
}
/*}}}*/

/*IntArrayInput management*/
void IntArrayInput::SetInput(int row,int numindices,int* values_in){/*{{{*/

	_assert_(this);
	_assert_(row>=0 && row<this->numberofelements_local);

	if(this->N[row] != numindices){
		if(this->values[row]) xDelete<int>(this->values[row]);
		this->values[row] = xNew<int>(numindices);
	}

	int *el_values = this->values[row];
	for(int i=0;i<numindices;i++) el_values[i] = values_in[i];

	this->N[row] = numindices;
}
/*}}}*/
void IntArrayInput::GetArray(int row,int** pvalues,int* pN){/*{{{*/

	_assert_(this);
	_assert_(row>=0 && row<this->numberofelements_local);
	if(pvalues){
		int* outvalues = xNew<int>(this->N[row]);
		xMemCpy<int>(outvalues,this->values[row],this->N[row]);
		*pvalues = outvalues;
	}
	if(pN){
		*pN = this->N[row];
	}
}
/*}}}*/
void IntArrayInput::GetArrayPtr(int row,int** pvalues,int* pN){/*{{{*/

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

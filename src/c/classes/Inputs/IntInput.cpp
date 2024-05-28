/*!\file IntInput.c
 * \brief: implementation of the IntInput object
 */

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "../classes.h"
#include "../../shared/shared.h"
#include "./IntInput.h"

/*IntInput constructors and destructor*/
IntInput::IntInput(){/*{{{*/
	this->size   = -1;
	this->values = NULL;
}
/*}}}*/
IntInput::IntInput(int size_in){/*{{{*/
	_assert_(size_in>0);
	_assert_(size_in<1e11);
	this->size   = size_in;
	this->values = xNew<int>(size_in);
}
/*}}}*/
IntInput::~IntInput(){/*{{{*/
	xDelete<int>(this->values);
}
/*}}}*/

/*Object virtual functions definitions:*/
Input* IntInput::copy() {/*{{{*/

	IntInput* output = new IntInput(this->size);
	xMemCpy<int>(output->values,this->values,this->size);

	return output;
}
/*}}}*/
void IntInput::DeepEcho(void){/*{{{*/

	_printf_("IntInput Echo:\n");
	_printf_("   Size:          "<<size<<"\n");
	printarray(this->values,this->size);
	//_printf_(setw(15)<<"   IntInput "<<setw(25)<<left<<EnumToStringx(this->enum_type)<<" "<<(value?"true":"false") << "\n");
}
/*}}}*/
void IntInput::Echo(void){/*{{{*/
	this->DeepEcho();
}
/*}}}*/
int  IntInput::Id(void){ return -1; }/*{{{*/
/*}}}*/
void IntInput::Marshall(MarshallHandle* marshallhandle){ /*{{{*/

	int object_enum = IntInputEnum;
   marshallhandle->call(object_enum);

	marshallhandle->call(this->size);
	if(this->size > 0){
		marshallhandle->call(this->values,this->size);
	}
	else this->values = NULL;

}
/*}}}*/
int  IntInput::ObjectEnum(void){/*{{{*/

	return IntInputEnum;

}
/*}}}*/

/*IntInput management*/
void IntInput::GetInput(int* pvalue,int index){/*{{{*/

	_assert_(index>=0); 
	_assert_(index<this->size); 

	*pvalue = this->values[index];
}
/*}}}*/
void IntInput::SetInput(int index,int value){/*{{{*/

	_assert_(index>=0); 
	_assert_(index<this->size); 

	this->values[index] = value;
}
/*}}}*/

/*Object functions*/

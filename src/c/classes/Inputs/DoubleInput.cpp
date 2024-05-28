/*!\file DoubleInput.c
 * \brief: implementation of the DoubleInput object
 */

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "../classes.h"
#include "../../shared/shared.h"
#include "./DoubleInput.h"

/*DoubleInput constructors and destructor*/
DoubleInput::DoubleInput(){/*{{{*/
	this->size   = -1;
	this->values = NULL;
}
/*}}}*/
DoubleInput::DoubleInput(int size_in){/*{{{*/
	_assert_(size_in>0);
	_assert_(size_in<1e11);
	this->size   = size_in;
	this->values = xNew<IssmDouble>(size_in);
}
/*}}}*/
DoubleInput::~DoubleInput(){/*{{{*/
	xDelete<IssmDouble>(this->values);
}
/*}}}*/

/*Object virtual functions definitions:*/
Input* DoubleInput::copy() {/*{{{*/

	DoubleInput* output = new DoubleInput(this->size);
	xMemCpy<IssmDouble>(output->values,this->values,this->size);

	return output;
}
/*}}}*/
void DoubleInput::DeepEcho(void){/*{{{*/

	_printf_("DoubleInput Echo:\n");
	_printf_("   Size:          "<<size<<"\n");
	printarray(this->values,this->size);
	//_printf_(setw(15)<<"   DoubleInput "<<setw(25)<<left<<EnumToStringx(this->enum_type)<<" "<<(value?"true":"false") << "\n");
}
/*}}}*/
void DoubleInput::Echo(void){/*{{{*/
	this->DeepEcho();
}
/*}}}*/
int  DoubleInput::Id(void){ return -1; }/*{{{*/
/*}}}*/
void DoubleInput::Marshall(MarshallHandle* marshallhandle){ /*{{{*/

	int object_enum = DoubleInputEnum;
   marshallhandle->call(object_enum);

	marshallhandle->call(this->size);
	if(this->size > 0){
		marshallhandle->call(this->values,this->size);
	}
	else this->values = NULL;

}
/*}}}*/
int  DoubleInput::ObjectEnum(void){/*{{{*/

	return DoubleInputEnum;

}
/*}}}*/

/*DoubleInput management*/
void DoubleInput::GetInput(IssmDouble* pvalue,int index){/*{{{*/

	_assert_(index>=0); 
	_assert_(index<this->size); 

	*pvalue = this->values[index];
}
/*}}}*/
void DoubleInput::SetInput(int index,IssmDouble value){/*{{{*/

	_assert_(index>=0); 
	_assert_(index<this->size); 

	this->values[index] = value;
}
/*}}}*/

/*Object functions*/

/*!\file BoolInput.c
 * \brief: implementation of the BoolInput object
 */

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "../classes.h"
#include "./BoolInput.h"
#include "../../shared/shared.h"

/*BoolInput constructors and destructor*/
BoolInput::BoolInput(){/*{{{*/
	this->size   = -1;
	this->values = NULL;
}
/*}}}*/
BoolInput::BoolInput(int size_in){/*{{{*/
	_assert_(size_in>0);
	_assert_(size_in<1e11);
	this->size   = size_in;
	this->values = xNew<bool>(size_in);
}
/*}}}*/
BoolInput::~BoolInput(){/*{{{*/
	xDelete<bool>(this->values);
}
/*}}}*/

/*Object virtual functions definitions:*/
Input* BoolInput::copy() {/*{{{*/

	_assert_(this->size);
	BoolInput* output = new BoolInput(this->size);
	xMemCpy<bool>(output->values,this->values,this->size);

	return output;

}
/*}}}*/
void BoolInput::DeepEcho(void){/*{{{*/

	_printf_("BoolInput Echo:\n");
	_printf_("   Size:          "<<size<<"\n");
	printarray(this->values,this->size);
}
/*}}}*/
void BoolInput::Echo(void){/*{{{*/
	this->DeepEcho();
}
/*}}}*/
int  BoolInput::Id(void){ return -1; }/*{{{*/
/*}}}*/
void BoolInput::Marshall(MarshallHandle* marshallhandle){ /*{{{*/

	int object_enum = BoolInputEnum;
   marshallhandle->call(object_enum);

	marshallhandle->call(this->size);
	if(this->size > 0){
		marshallhandle->call(this->values,this->size);
	}
	else this->values = NULL;

}
/*}}}*/
int  BoolInput::ObjectEnum(void){/*{{{*/

	return BoolInputEnum;

}
/*}}}*/

/*BoolInput management*/
void BoolInput::GetInput(bool* pvalue,int index){/*{{{*/

	_assert_(index>=0); 
	_assert_(index<this->size); 

	*pvalue = this->values[index];
}
/*}}}*/
void BoolInput::SetInput(int index,bool value){/*{{{*/

	_assert_(index>=0); 
	_assert_(index<this->size); 

	this->values[index] = value;
}
/*}}}*/

/*Object functions*/

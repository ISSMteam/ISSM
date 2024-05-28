/*!\file IntParam.c
 * \brief: implementation of the IntParam object
 */

/*header files: */
/*{{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "../classes.h"
#include "../../shared/shared.h"
/*}}}*/

/*IntParam constructors and destructor*/
IntParam::IntParam(){/*{{{*/
	return;
}
/*}}}*/
IntParam::IntParam(int in_enum_type,IssmInt in_value){/*{{{*/

	enum_type=in_enum_type;
	value=in_value;
}
/*}}}*/
IntParam::~IntParam(){/*{{{*/
	return;
}
/*}}}*/

/*Object virtual functions definitions:*/
Param* IntParam::copy() {/*{{{*/

	return new IntParam(this->enum_type,this->value);

}
/*}}}*/
void IntParam::DeepEcho(void){/*{{{*/

	_printf_(setw(22)<<"   IntParam "<<setw(35)<<left<<EnumToStringx(this->enum_type)<<" "<<this->value<<"\n");
}
/*}}}*/
void IntParam::Echo(void){/*{{{*/
	this->DeepEcho();
}
/*}}}*/
int  IntParam::Id(void){ return -1; }/*{{{*/
/*}}}*/
void IntParam::Marshall(MarshallHandle* marshallhandle){ /*{{{*/

	int object_enum = IntParamEnum;
   marshallhandle->call(object_enum);

	marshallhandle->call(this->enum_type);
	marshallhandle->call(this->value);
}/*}}}*/
int  IntParam::ObjectEnum(void){/*{{{*/

	return IntParamEnum;

}
/*}}}*/

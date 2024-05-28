/*!\file BoolParam.c
 * \brief: implementation of the BoolParam object
 */

/*header files: */
/*{{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "../classes.h"
#include "shared/shared.h"
/*}}}*/

/*BoolParam constructors and destructor*/
BoolParam::BoolParam(){/*{{{*/
	return;
}
/*}}}*/
BoolParam::BoolParam(int in_enum_type,bool in_value){/*{{{*/

	enum_type=in_enum_type;
	value=in_value;
}
/*}}}*/
BoolParam::~BoolParam(){/*{{{*/
	return;
}
/*}}}*/

/*Object virtual functions definitions:*/
Param* BoolParam::copy() {/*{{{*/

	return new BoolParam(this->enum_type,this->value);

}
/*}}}*/
void BoolParam::DeepEcho(void){/*{{{*/
	_printf_(setw(22)<<"   BoolParam "<<setw(35)<<left<<EnumToStringx(this->enum_type)<<" "<<(value?"true":"false")<<"\n");
}
/*}}}*/
void BoolParam::Echo(void){/*{{{*/
	this->DeepEcho();
}
/*}}}*/
int    BoolParam::Id(void){ return -1; }/*{{{*/
/*}}}*/
void BoolParam::Marshall(MarshallHandle* marshallhandle){ /*{{{*/

	int object_enum = BoolParamEnum;
	marshallhandle->call(object_enum);
	marshallhandle->call(this->enum_type);
	marshallhandle->call(this->value);

}/*}}}*/
int BoolParam::ObjectEnum(void){/*{{{*/

	return BoolParamEnum;

}
/*}}}*/

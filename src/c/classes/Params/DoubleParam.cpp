/*!\file DoubleParam.c
 * \brief: implementation of the DoubleParam object
 */

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "../classes.h"
#include "shared/shared.h"

/*DoubleParam constructors and destructor*/
DoubleParam::DoubleParam(){/*{{{*/
	return;
}
/*}}}*/
DoubleParam::DoubleParam(int in_enum_type,IssmDouble in_value){/*{{{*/

	enum_type=in_enum_type;
	value=in_value;
}
/*}}}*/
DoubleParam::~DoubleParam(){/*{{{*/
	return;
}
/*}}}*/

/*Object virtual functions definitions:*/
Param* DoubleParam::copy() {/*{{{*/

	return new DoubleParam(this->enum_type,this->value);

}
/*}}}*/
void DoubleParam::DeepEcho(void){/*{{{*/

	_printf_(setw(22)<<"   DoubleParam "<<setw(35)<<left<<EnumToStringx(this->enum_type)<<" "<<this->value<<"\n");
}
/*}}}*/
void DoubleParam::Echo(void){/*{{{*/
	this->DeepEcho();
}
/*}}}*/
int  DoubleParam::Id(void){ return -1; }/*{{{*/
/*}}}*/
void DoubleParam::Marshall(MarshallHandle* marshallhandle){ /*{{{*/

	int object_enum = DoubleParamEnum;
   marshallhandle->call(object_enum);
	marshallhandle->call(this->enum_type);
	marshallhandle->call(this->value);

}
/*}}}*/
int  DoubleParam::ObjectEnum(void){/*{{{*/

	return DoubleParamEnum;

}
/*}}}*/

/*DoubleParam virtual functions definitions: */
void DoubleParam::GetParameterValue(int* pinteger){/*{{{*/
	_error_("Double param of enum " << enum_type << " (" << EnumToStringx(enum_type) << ") cannot return an integer");
}
/*}}}*/
void DoubleParam::GetParameterValue(bool* pbool){/*{{{*/
	_error_("Double param of enum " << enum_type << " (" << EnumToStringx(enum_type) << ") cannot return an bool");
}
/*}}}*/
void DoubleParam::GetParameterValue(int** pintarray,int* pM){/*{{{*/
	_error_("Double param of enum " << enum_type << " (" << EnumToStringx(enum_type) << ") cannot return an array of integers");
}
/*}}}*/
void DoubleParam::GetParameterValue(int** pintarray,int* pM,int* pN){/*{{{*/
	_error_("Double param of enum " << enum_type << " (" << EnumToStringx(enum_type) << ") cannot return an array of integers");
}
/*}}}*/
void DoubleParam::GetParameterValue(IssmDouble** pIssmDoublearray,int* pM){/*{{{*/
	_error_("Double param of enum " << enum_type << " (" << EnumToStringx(enum_type) << ") cannot return an array of IssmDouble");
}
/*}}}*/
void DoubleParam::GetParameterValue(IssmDouble** pIssmDoublearray,int* pM,int* pN){/*{{{*/
	_error_("Double param of enum " << enum_type << " (" << EnumToStringx(enum_type) << ") cannot return an array of IssmDouble");
}
/*}}}*/

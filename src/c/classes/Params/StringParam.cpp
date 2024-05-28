/*!\file StringParam.c
 * \brief: implementation of the StringParam object
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

/*StringParam constructors and destructor*/
StringParam::StringParam(){/*{{{*/
	return;
}
/*}}}*/
StringParam::StringParam(int in_enum_type,char* in_value){/*{{{*/

	enum_type=in_enum_type;
	value=xNew<char>(strlen(in_value)+1);
	xMemCpy<char>(value,in_value,(strlen(in_value)+1));

}
/*}}}*/
StringParam::~StringParam(){/*{{{*/
	xDelete<char>(value);
}
/*}}}*/

/*Object virtual functions definitions:*/
Param* StringParam::copy() {/*{{{*/

	return new StringParam(this->enum_type,this->value);

}
/*}}}*/
void StringParam::DeepEcho(void){/*{{{*/
	_printf_(setw(22)<<"   StringParam "<<setw(35)<<left<<EnumToStringx(this->enum_type)<<" "<<this->value<<"\n");
}
/*}}}*/
void StringParam::Echo(void){/*{{{*/
	this->DeepEcho();
}
/*}}}*/
int    StringParam::Id(void){ return -1; }/*{{{*/
/*}}}*/
void StringParam::Marshall(MarshallHandle* marshallhandle){ /*{{{*/

	int object_enum = StringParamEnum;
   marshallhandle->call(object_enum);
	marshallhandle->call(this->enum_type);
	marshallhandle->call(this->value);
}
/*}}}*/
int StringParam::ObjectEnum(void){/*{{{*/

	return StringParamEnum;

}
/*}}}*/

/*StringParam virtual functions definitions: */
void  StringParam::GetParameterValue(char** pstring){/*{{{*/

	char* outstring=NULL;
	int   stringsize;

	stringsize=strlen(this->value)+1;

	outstring=xNew<char>(stringsize);
	xMemCpy<char>(outstring,this->value,stringsize);

	*pstring=outstring;

}
/*}}}*/
void  StringParam::SetValue(char* string){/*{{{*/

	int   stringsize;

	/*avoid leak: */
	xDelete<char>(this->value);

	/*copy: */
	stringsize=strlen(string)+1;
	this->value=xNew<char>(stringsize);
	xMemCpy<char>(this->value,string,stringsize);

}
/*}}}*/

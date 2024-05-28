/*!\file FileParam.c
 * \brief: implementation of the FileParam object
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

/*FileParam constructors and destructor*/
FileParam::FileParam(){/*{{{*/
	return;
}
/*}}}*/
FileParam::FileParam(int in_enum_type,FILE* in_value){/*{{{*/

	enum_type=in_enum_type;
	value=in_value;
}
/*}}}*/
FileParam::~FileParam(){/*{{{*/
	return;
}
/*}}}*/

/*Object virtual functions definitions:*/
Param* FileParam::copy() {/*{{{*/

	return new FileParam(this->enum_type,this->value);

}
/*}}}*/
void FileParam::DeepEcho(void){/*{{{*/

	_printf_(setw(22)<<"   FileParam "<<setw(35)<<left<<EnumToStringx(this->enum_type)<<" "<<this->value<<"\n");
}
/*}}}*/
void FileParam::Echo(void){/*{{{*/
	this->DeepEcho();
}
/*}}}*/
int  FileParam::Id(void){ return -1; }/*{{{*/
/*}}}*/
void FileParam::Marshall(MarshallHandle* marshallhandle){ /*{{{*/

	int object_enum = FileParamEnum;
   marshallhandle->call(object_enum);

	marshallhandle->call(this->enum_type);
	marshallhandle->call(this->value);

	if(marshallhandle->OperationNumber()==MARSHALLING_LOAD){
		this->value=NULL; //meaningless file pointer!
	}

}
/*}}}*/
int  FileParam::ObjectEnum(void){/*{{{*/

	return FileParamEnum;

}
/*}}}*/

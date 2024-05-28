/*!\file StringArrayParam.c
 * \brief: implementation of the StringArrayParam object
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

/*StringArrayParam constructors and destructor*/
StringArrayParam::StringArrayParam(){/*{{{*/
	return;
}
/*}}}*/
StringArrayParam::StringArrayParam(int in_enum_type,char** in_values, int in_numstrings){/*{{{*/

	this->enum_type=in_enum_type;
	this->numstrings=in_numstrings;

	if(numstrings){
		this->value=xNew<char*>(numstrings);
		for(int i=0;i<numstrings;i++){
			int   size=strlen(in_values[i])+1;
			char* string=xNew<char>(size);
			xMemCpy<char>(string,in_values[i],size);
			this->value[i]=string;
		}
	}
	else{
		this->value=NULL;
	}

}
/*}}}*/
StringArrayParam::~StringArrayParam(){/*{{{*/

	int i;

	char* string=NULL;
	for(i=0;i<this->numstrings;i++){
		string=value[i];
		xDelete<char>(string);
	}
	xDelete<char*>(value);
}
/*}}}*/

/*Object virtual functions definitions:*/
Param* StringArrayParam::copy() {/*{{{*/

	return new StringArrayParam(this->enum_type,this->value,this->numstrings);

}
/*}}}*/
void StringArrayParam::DeepEcho(void){/*{{{*/

	_printf_(setw(22)<<"   StringArrayParam "<<setw(35)<<left<<EnumToStringx(this->enum_type)<<" {");
	for(int i=0;i<this->numstrings;i++) _printf_(" '"<<this->value[i]<<"'");
	_printf_("}\n");
}
/*}}}*/
void StringArrayParam::Echo(void){/*{{{*/
	this->DeepEcho();
}
/*}}}*/
int    StringArrayParam::Id(void){ return -1; }/*{{{*/
/*}}}*/
void StringArrayParam::Marshall(MarshallHandle* marshallhandle){ /*{{{*/

	int object_enum = StringArrayParamEnum;
   marshallhandle->call(object_enum);

	marshallhandle->call(this->enum_type);
	marshallhandle->call(this->numstrings);

	if(this->numstrings){
		if(marshallhandle->OperationNumber()==MARSHALLING_LOAD){
			this->value=xNew<char*>(this->numstrings);
		}
		for(int i=0;i<numstrings;i++){
			marshallhandle->call(this->value[i]);
		}
	}
	else{
		this->value=NULL;
	}
}
/*}}}*/
int StringArrayParam::ObjectEnum(void){/*{{{*/

	return StringArrayParamEnum;

}
/*}}}*/

/*StringArrayParam virtual functions definitions: */
void  StringArrayParam::GetParameterValue(char*** pstringarray,int* pM){/*{{{*/

	int   i;
	char** outstrings=NULL;
	int   M;
	char* string=NULL;
	char* string2=NULL;
	int   stringsize;

	M=this->numstrings;
	if(this->numstrings){
		outstrings=xNew<char*>(this->numstrings);

		for(i=0;i<this->numstrings;i++){
			string=this->value[i];
			stringsize=strlen(string)+1;

			string2=xNew<char>(stringsize);
			xMemCpy<char>(string2,string,stringsize);

			outstrings[i]=string2;
		}
	}
	else outstrings=NULL;

	/*Assign output pointers:*/
	if(pM) *pM=M;
	*pstringarray=outstrings;
}
/*}}}*/
void  StringArrayParam::SetValue(char** stringarray,int M){/*{{{*/

	int   i;
	char *string     = NULL;
	char *string2    = NULL;
	int   stringsize;

	/*first, avoid leak: */
	for(i=0;i<this->numstrings;i++){
		string=this->value[i];
		xDelete<char>(string);
	}
	xDelete<char*>(this->value);

	/*copy: */
	this->numstrings=M;
	this->value=xNew<char*>(this->numstrings);
	for(i=0;i<this->numstrings;i++){
		string=stringarray[i];
		stringsize=strlen(string)+1;

		string2=xNew<char>(stringsize);
		xMemCpy<char>(string2,string,stringsize);

		this->value[i]=string2;
	}
}
/*}}}*/

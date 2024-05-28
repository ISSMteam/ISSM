/*!\file IntVecParam.c
 * \brief: implementation of the IntVecParam object
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

/*IntVecParam constructors and destructor*/
IntVecParam::IntVecParam(){/*{{{*/
	return;
}
/*}}}*/
IntVecParam::IntVecParam(int in_enum_type,int* in_values, int in_M){/*{{{*/

	enum_type=in_enum_type;
	M=in_M;

	if(M){
		values=xNew<int>(M);
		xMemCpy<int>(values,in_values,M);
	}
	else values=NULL;
}
/*}}}*/
IntVecParam::IntVecParam(int in_enum_type,IssmDouble* in_values, int in_M){/*{{{*/

	enum_type=in_enum_type;
	M=in_M;

	if(M){
		values=xNew<int>(M);
		for(int i=0;i<in_M;i++) values[i]=reCast<int>(in_values[i]);
	}
	else values=NULL;
}
/*}}}*/
IntVecParam::~IntVecParam(){/*{{{*/
	xDelete<int>(values);
	return;
}
/*}}}*/

/*Object virtual functions definitions:*/
Param* IntVecParam::copy() {/*{{{*/

	return new IntVecParam(this->enum_type,this->values,this->M);

}
/*}}}*/
void IntVecParam::DeepEcho(void){/*{{{*/
	_printf_(setw(22)<<"   IntVecParam "<<setw(35)<<left<<EnumToStringx(this->enum_type)<<" [");
	for(int i=0;i<this->M;i++) _printf_(" "<<this->values[i]);
	_printf_("]\n");
}
/*}}}*/
void IntVecParam::Echo(void){/*{{{*/

	this->DeepEcho();
}
/*}}}*/
int  IntVecParam::Id(void){ return -1; }/*{{{*/
/*}}}*/
void IntVecParam::Marshall(MarshallHandle* marshallhandle){ /*{{{*/

	int object_enum = IntVecParamEnum;
   marshallhandle->call(object_enum);

	marshallhandle->call(this->enum_type);
	marshallhandle->call(this->M);
	if(M){ 
		marshallhandle->call(this->values,M);
	}
	else{
		this->values=NULL;
	}

}
/*}}}*/
int  IntVecParam::ObjectEnum(void){/*{{{*/

	return IntVecParamEnum;

}
/*}}}*/

/*IntVecParam virtual functions definitions: */
void  IntVecParam::GetParameterValue(int** pintarray,int* pM){/*{{{*/
	int* output=NULL;

	if(M){
		output=xNew<int>(M);
		xMemCpy<int>(output,values,M);
	}

	/*Assign output pointers:*/
	if(pM) *pM=M;
	*pintarray=output;
}
/*}}}*/
void  IntVecParam::SetValue(int* intarray,int in_M){/*{{{*/

	/*avoid leak: */
	xDelete<int>(this->values);

	if(in_M){
		this->values=xNew<int>(in_M);
		xMemCpy<int>(this->values,intarray,in_M);
	}
	else this->values=NULL;

	this->M=in_M;
}
/*}}}*/

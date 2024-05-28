/*!\file DoubleMatParam.c
 * \brief: implementation of the DoubleMatParam object
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

/*DoubleMatParam constructors and destructor*/
DoubleMatParam::DoubleMatParam(){/*{{{*/
	return;
}
/*}}}*/
DoubleMatParam::DoubleMatParam(int in_enum_type,IssmDouble* in_value, int in_M,int in_N){/*{{{*/

	enum_type=in_enum_type;
	M=in_M;
	N=in_N;

	value=xNew<IssmDouble>(M*N);
	xMemCpy<IssmDouble>(value,in_value,M*N);
}
/*}}}*/
DoubleMatParam::~DoubleMatParam(){/*{{{*/
	xDelete<IssmDouble>(value);
	return;
}
/*}}}*/

/*Object virtual functions definitions:*/
void DoubleMatParam::Echo(void){/*{{{*/

	_printf_("DoubleMatParam:\n");
	_printf_("   enum: " << this->enum_type << " (" << EnumToStringx(this->enum_type) << ")\n");
	_printf_("   matrix size: " << this->M << "x" << this->N << "\n");

}
/*}}}*/
void DoubleMatParam::DeepEcho(void){/*{{{*/

	int i,j;

	_printf_("DoubleMatParam:\n");
	_printf_("   enum: " << this->enum_type << " (" << EnumToStringx(this->enum_type) << ")\n");
	_printf_("   matrix size: " << this->M << "x" << this->N << "\n");
	for(i=0;i<this->M;i++){
		for(j=0;j<this->N;j++){
			_printf_(i << " " << j << " " << *(this->value+N*i+j) << "\n");
		}
	}
}
/*}}}*/
int    DoubleMatParam::Id(void){ return -1; }/*{{{*/
/*}}}*/
int DoubleMatParam::ObjectEnum(void){/*{{{*/

	return DoubleMatParamEnum;

}
/*}}}*/
Param* DoubleMatParam::copy() {/*{{{*/

	return new DoubleMatParam(this->enum_type,this->value,this->M,this->N);

}
/*}}}*/
void DoubleMatParam::Marshall(MarshallHandle* marshallhandle){ /*{{{*/

	int object_enum = DoubleMatParamEnum;
   marshallhandle->call(object_enum);

	marshallhandle->call(this->enum_type);
	marshallhandle->call(this->M);
	marshallhandle->call(this->N);
	marshallhandle->call(this->value,M*N);
}
/*}}}*/

/*DoubleMatParam virtual functions definitions: */
void  DoubleMatParam::GetParameterValue(IssmDouble** pIssmDoublearray,int* pM,int* pN){/*{{{*/
	IssmDouble* output=NULL;

	output=xNew<IssmDouble>(M*N);
	xMemCpy<IssmDouble>(output,value,M*N);

	/*Assign output pointers:*/
	if(pM) *pM=M;
	if(pN) *pN=N;
	*pIssmDoublearray=output;
}
/*}}}*/
void  DoubleMatParam::GetParameterValue(int** pintarray,int* pM,int* pN){/*{{{*/
	_error_("DoubleMat of enum " << enum_type << " (" << EnumToStringx(enum_type) << ") cannot return an array of int");
}
/*}}}*/
void  DoubleMatParam::SetValue(IssmDouble* IssmDoublearray,int in_M,int in_N){/*{{{*/

	/*avoid leak: */
	xDelete<IssmDouble>(this->value);

	this->value=xNew<IssmDouble>(in_M*in_N);
	xMemCpy<IssmDouble>(this->value,IssmDoublearray,in_M*in_N);

	this->M=in_M;
	this->N=in_N;
}
/*}}}*/

/*DoubleMatParam specific routines:*/
void  DoubleMatParam::GetParameterValueByPointer(IssmDouble** pIssmDoublearray,int* pM,int* pN){/*{{{*/

	/*Assign output pointers:*/
	if(pM) *pM=M;
	if(pN) *pN=N;
	*pIssmDoublearray=value;
}
/*}}}*/

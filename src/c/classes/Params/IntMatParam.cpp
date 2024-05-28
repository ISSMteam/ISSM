/*!\file IntMatParam.c
 * \brief: implementation of the IntMatParam object
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

/*IntMatParam constructors and destructor*/
IntMatParam::IntMatParam(){/*{{{*/
	return;
}
/*}}}*/
IntMatParam::IntMatParam(int in_enum_type,int* in_value, int in_M,int in_N){/*{{{*/

	enum_type=in_enum_type;
	M=in_M;
	N=in_N;

	value=xNew<int>(M*N);
	xMemCpy<int>(value,in_value,M*N);
}
/*}}}*/
IntMatParam::~IntMatParam(){/*{{{*/
	xDelete<int>(value);
	return;
}
/*}}}*/

/*Object virtual functions definitions:*/
Param* IntMatParam::copy() {/*{{{*/

	return new IntMatParam(this->enum_type,this->value,this->M,this->N);

}
/*}}}*/
void IntMatParam::DeepEcho(void){/*{{{*/

	int i,j;

	_printf_("IntMatParam:\n");
	_printf_("   enum: " << this->enum_type << " (" << EnumToStringx(this->enum_type) << ")\n");
	_printf_("   matrix size: " << this->M << "x" << this->N << "\n");
	for(i=0;i<this->M;i++){
		for(j=0;j<this->N;j++){
			_printf_("(" << i << "," << j << ") " << *(this->value+N*i+j) << "\n");
		}
	}
}
/*}}}*/
void IntMatParam::Echo(void){/*{{{*/

	_printf_("IntMatParam:\n");
	_printf_("   enum: " << this->enum_type << " (" << EnumToStringx(this->enum_type) << ")\n");
	_printf_("   matrix size: " << this->M << "x" << this->N << "\n");

}
/*}}}*/
int  IntMatParam::Id(void){ return -1; }/*{{{*/
/*}}}*/
void IntMatParam::Marshall(MarshallHandle* marshallhandle){ /*{{{*/

	int object_enum = IntMatParamEnum;
   marshallhandle->call(object_enum);

	marshallhandle->call(this->enum_type);
	marshallhandle->call(this->M);
	marshallhandle->call(this->N);
	marshallhandle->call(this->value,M*N);
}
/*}}}*/
int  IntMatParam::ObjectEnum(void){/*{{{*/

	return IntMatParamEnum;

}
/*}}}*/

/*IntMatParam virtual functions definitions: */
void  IntMatParam::GetParameterValue(int** pintarray,int* pM,int* pN){/*{{{*/
	int* output=NULL;

	output=xNew<int>(M*N);
	xMemCpy<int>(output,value,M*N);

	/*Assign output pointers:*/
	if(pM) *pM=M;
	if(pN) *pN=N;
	*pintarray=output;
}
/*}}}*/
void  IntMatParam::SetValue(int* intarray,int in_M,int in_N){/*{{{*/

	/*avoid leak: */
	xDelete<int>(this->value);

	this->value=xNew<int>(in_M*in_N);
	xMemCpy<int>(this->value,intarray,in_M*in_N);

	this->M=in_M;
	this->N=in_N;
}
/*}}}*/

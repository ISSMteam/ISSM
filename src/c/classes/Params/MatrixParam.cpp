/*!\file MatrixParam.c
 * \brief: implementation of the MatrixParam object
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

/*MatrixParam constructors and destructor*/
MatrixParam::MatrixParam(){/*{{{*/
	return;
}
/*}}}*/
MatrixParam::MatrixParam(int in_enum_type,Matrix<IssmDouble>* in_value){/*{{{*/

	enum_type=in_enum_type;
	value=NULL;

	if(in_value){
		value=in_value->Duplicate();
	}
}
/*}}}*/
MatrixParam::~MatrixParam(){/*{{{*/
	delete value;
}
/*}}}*/

/*Object virtual functions definitions:*/
Param* MatrixParam::copy() {/*{{{*/

	return new MatrixParam(this->enum_type,this->value);

}
/*}}}*/
void MatrixParam::DeepEcho(void){/*{{{*/

	_printf_("MatrixParam:\n");
	_printf_("   enum: " << this->enum_type << " (" << EnumToStringx(this->enum_type) << ")\n");
	this->value->Echo();
}
/*}}}*/
void MatrixParam::Echo(void){/*{{{*/

	_printf_("MatrixParam:\n");
	_printf_("   enum: " << this->enum_type << " (" << EnumToStringx(this->enum_type) << ")\n");

}
/*}}}*/
int  MatrixParam::Id(void){ return -1; }/*{{{*/
/*}}}*/
int MatrixParam::ObjectEnum(void){/*{{{*/

	return MatrixParamEnum;

}
/*}}}*/

/*MatrixParam virtual functions definitions: */
void  MatrixParam::GetParameterValue(Matrix<IssmDouble>** poutput){/*{{{*/
	Matrix<IssmDouble>* output=NULL;

	if(value){
		output=value->Duplicate();
	}
	*poutput=output;
}
/*}}}*/
void  MatrixParam::SetValue(Matrix<IssmDouble>* matrix){/*{{{*/

	/*avoid leak: */
	delete value;

	/*copy: */
	value=matrix->Duplicate();
}
/*}}}*/

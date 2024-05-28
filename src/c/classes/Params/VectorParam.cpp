/*!\file VectorParam.c
 * \brief: implementation of the VectorParam object
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

/*VectorParam constructors and destructor*/
VectorParam::VectorParam(){/*{{{*/
	return;
}
/*}}}*/
VectorParam::VectorParam(int in_enum_type,Vector<IssmDouble>* in_value){/*{{{*/

	enum_type=in_enum_type;

	value=NULL;

	if(in_value){
		value=in_value->Duplicate();
		in_value->Copy(value);
	}
}
/*}}}*/
VectorParam::~VectorParam(){/*{{{*/
	delete value;
}
/*}}}*/

/*Object virtual functions definitions:*/
Param* VectorParam::copy() {/*{{{*/

	return new VectorParam(this->enum_type,this->value);

}
/*}}}*/
void VectorParam::DeepEcho(void){/*{{{*/

	_printf_("VectorParam:\n");
	_printf_("   enum: " << this->enum_type << " (" << EnumToStringx(this->enum_type) << ")\n");
	value->Echo();
}
/*}}}*/
void VectorParam::Echo(void){/*{{{*/

	_printf_("VectorParam:\n");
	_printf_("   enum: " << this->enum_type << " (" << EnumToStringx(this->enum_type) << ")\n");

}
/*}}}*/
int  VectorParam::Id(void){ return -1; }/*{{{*/
/*}}}*/
int VectorParam::ObjectEnum(void){/*{{{*/

	return VectorParamEnum;

}
/*}}}*/

/*VectorParam virtual functions definitions: */
void  VectorParam::GetParameterValue(Vector<IssmDouble>** poutput){/*{{{*/
	Vector<IssmDouble>*  output=NULL;

	if(value){
		output=value->Duplicate();
		value->Copy(output);
	}
	*poutput=output;
}
/*}}}*/
void  VectorParam::SetValue(Vector<IssmDouble>* vector){/*{{{*/

	/*avoid leak: */
	delete value;

	/*copy: */
	value=vector->Duplicate();
	vector->Copy(value);
}
/*}}}*/

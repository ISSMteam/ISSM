/*!\file EmulatorParam.c
 * \brief: implementation of the EmulatorParam object
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
#include <pybind11/numpy.h>

/*EmulatorParam constructors and destructor*/
EmulatorParam::EmulatorParam(){/*{{{*/
	value=NULL;
	return;
}
/*}}}*/
EmulatorParam::EmulatorParam(int in_enum_type, char* pt_path_in){/*{{{*/

	this->enum_type=in_enum_type;

	/*Copy path to emulator*/
	this->pt_path = xNew<char>(strlen(pt_path_in)+1);
	xMemCpy<char>(this->pt_path, pt_path_in,(strlen(pt_path_in)+1));

	/*Activate interpretor*/
	_error_("not finished yet");
}
/*}}}*/
EmulatorParam::~EmulatorParam(){/*{{{*/
	xDelete<char>(this->pt_path);
}
/*}}}*/

/*Object virtual functions definitions:*/
Param* EmulatorParam::copy() {/*{{{*/

	_error_("not implemented");

}
/*}}}*/
void EmulatorParam::DeepEcho(void){/*{{{*/

	_printf_(setw(22)<<"   EmulatoParam "<<setw(35)<<left<<EnumToStringx(this->enum_type)<<" guard: "<<this->guard<<", mod: "<< this->mod <<\n");
}
/*}}}*/
void EmulatorParam::Echo(void){/*{{{*/
	this->DeepEcho();
}
/*}}}*/
int  EmulatorParam::Id(void){ return -1; }/*{{{*/
/*}}}*/
void EmulatorParam::Marshall(MarshallHandle* marshallhandle){ /*{{{*/

	_error_("Not implemented yet");

}
/*}}}*/
int EmulatorParam::ObjectEnum(void){/*{{{*/

	return EmulatorParamEnum;

}
/*}}}*/

/*EmulatorParam virtual functions definitions: */

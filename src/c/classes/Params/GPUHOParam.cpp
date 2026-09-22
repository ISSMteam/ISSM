/*!\file GPUHOParam.c
 * \brief: implementation of the GPUHOParam object
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

/*GPUHOParam constructors and destructor*/
GPUHOParam::GPUHOParam(){/*{{{*/
	this->Kff_dofarraysize = 0;
	this->Pf_valuesarraysize = 0;
	this->Kff_valuesarraysize = 0;
	this->PDStressHO_valuesarraysize = 0;
	this->PDStressHO_basisarraysize = 0;
	this->PDStressHO_factorarraysize = 0;
	return;
}
/*}}}*/
GPUHOParam::~GPUHOParam(){/*{{{*/
}
/*}}}*/

/*Object virtual functions definitions:*/
Param* GPUHOParam::copy() {/*{{{*/

	_error_("not implemented");

}
/*}}}*/
void GPUHOParam::DeepEcho(void){/*{{{*/

	_error_("not implemented");

}
/*}}}*/
void GPUHOParam::Echo(void){/*{{{*/
	this->DeepEcho();
}
/*}}}*/
void GPUHOParam::Marshall(MarshallHandle* marshallhandle){ /*{{{*/

	_error_("Not implemented yet");

}
/*}}}*/

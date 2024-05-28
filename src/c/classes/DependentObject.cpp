/*!\file DependentObject.c
 * \brief: implementation of the DependentObject object
 */

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./classes.h"
#include "shared/shared.h"
#include "../modules/modules.h"

/*DependentObject constructors and destructor*/
DependentObject::DependentObject(){/*{{{*/
	this->name=NULL;
	this->response_value=0.;
}
/*}}}*/
DependentObject::DependentObject(char* in_name){/*{{{*/

	this->name=xNew<char>(strlen(in_name)+1); xMemCpy<char>(this->name,in_name,strlen(in_name)+1);
	this->response_value=0.;

}/*}}}*/
DependentObject::DependentObject(char* in_name,IssmDouble in_response){/*{{{*/

	this->name=xNew<char>(strlen(in_name)+1); xMemCpy<char>(this->name,in_name,strlen(in_name)+1);
	this->response_value=in_response;

}/*}}}*/
DependentObject::~DependentObject(){ //destructor/*{{{*/
	xDelete<char>(this->name);
}/*}}}*/

/*Object virtual functions definitions:*/
Object* DependentObject::copy(void) { /*{{{*/
	return new DependentObject(name,response_value);
} /*}}}*/
void DependentObject::DeepEcho(void){/*{{{*/
	this->Echo();
}
/*}}}*/
void DependentObject::Echo(void){/*{{{*/

	_printf_("DependentObject:\n");
	_printf_("   name: " << this->name << "\n");
	_printf_("   response_value: " << this->response_value<< "\n");
}
/*}}}*/
int  DependentObject::Id(void){ return -1; }/*{{{*/
/*}}}*/
int  DependentObject::ObjectEnum(void){/*{{{*/

	return DependentObjectEnum;

}
/*}}}*/
void DependentObject::Marshall(MarshallHandle* marshallhandle){/*{{{*/

	int object_enum = DependentObjectEnum;
	marshallhandle->call(object_enum);

	/*Marshall name (tricky)*/
	marshallhandle->call(this->name);

	marshallhandle->call(this->response_value);
}/*}}}*/

/*DependentObject methods: */
void  DependentObject::RecordResponsex(FemModel* femmodel){/*{{{*/

	int ierr;

	/*Is this some special type of response for which we need to go in the output definitions? :*/
	if(StringToEnumx(this->name,false)==-1){
		ierr = OutputDefinitionsResponsex(&this->response_value, femmodel, this->name);
		if(ierr) _error_("Could not find response "<<this->name);
	}
	else{
		femmodel->Responsex(&this->response_value, this->name);
	}
}
/*}}}*/
IssmDouble DependentObject::GetValue(void){/*{{{*/
	return this->response_value;
}
/*}}}*/
void DependentObject::ResetResponseValue(){/*{{{*/
	this->response_value=0.;
}
/*}}}*/

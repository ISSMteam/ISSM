/*!\file ControlInput.c
 * \brief: implementation of the ControlInput object
 */

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./ControlInput.h"
#include "./ElementInput.h"
#include "./TriaInput.h"
#include "./PentaInput.h"
#include "./TransientInput.h"
//#include "../../toolkits/objects/Vector.h"

/*ControlInput constructors and destructor*/
ControlInput::ControlInput(){/*{{{*/
	control_id  = 0;
	values      = NULL;
	minvalues   = NULL;
	maxvalues   = NULL;
}
/*}}}*/
ControlInput::ControlInput(int nbe, int nbv,int input_layout_enum,int interp,int id){/*{{{*/

	this->control_id  = id;
	this->layout_enum = input_layout_enum;

	switch(this->layout_enum){
		case TriaInputEnum:
			this->values     =new TriaInput(nbe,nbv,interp);
			this->minvalues  =new TriaInput(nbe,nbv,interp);
			this->maxvalues  =new TriaInput(nbe,nbv,interp);
			break;
		case PentaInputEnum:
			this->values     =new PentaInput(nbe,nbv,interp);
			this->minvalues  =new PentaInput(nbe,nbv,interp);
			this->maxvalues  =new PentaInput(nbe,nbv,interp);
			break;
		default:
			_error_("Input of Enum \"" << EnumToStringx(input_layout_enum) << "\" not supported yet by ControlInput");
	}
}
/*}}}*/
ControlInput::ControlInput(int enum_in,int nbe, int nbv,int id,IssmDouble* times, int numtimes){/*{{{*/

	this->enum_type   = enum_in;
	this->control_id  = id;
	this->layout_enum = TransientInputEnum; /*Tria or Penta?*/

	this->values     =new TransientInput(enum_in,nbe,nbv,times,numtimes);
	this->minvalues  =new TransientInput(enum_in,nbe,nbv,times,numtimes);
	this->maxvalues  =new TransientInput(enum_in,nbe,nbv,times,numtimes);
}
/*}}}*/
ControlInput::~ControlInput(){/*{{{*/
	delete values;
	delete minvalues;
	delete maxvalues;
}
/*}}}*/

/*Object virtual functions definitions:*/
Input* ControlInput::copy() {/*{{{*/

	ControlInput* output=NULL;

	output = new ControlInput();
	output->enum_type=this->enum_type;
	output->control_id=this->control_id;
	output->layout_enum = this->layout_enum;

	if(values)      output->values      = this->values->copy();
	if(minvalues)   output->minvalues   = this->minvalues->copy();
	if(maxvalues)   output->maxvalues   = this->maxvalues->copy();

	return output;
}
/*}}}*/
void ControlInput::Configure(Parameters* params){/*{{{*/
	this->values->Configure(params);
	this->minvalues->Configure(params);
	this->maxvalues->Configure(params);
}
/*}}}*/
void ControlInput::DeepEcho(void){/*{{{*/

	_printf_("ControlInput:\n");
	_printf_(setw(15)<<"   ControlInput "<<setw(25)<<left<<EnumToStringx(this->enum_type)<<"\n");
	_printf_(setw(15)<<"   Layout       "<<setw(25)<<left<<EnumToStringx(this->layout_enum)<<"\n");
	_printf_("---values: \n");     if (values)      values->Echo();
	_printf_("---minvalues: \n");  if (minvalues)   minvalues->Echo();
	_printf_("---maxvalues: \n");  if (maxvalues)   maxvalues->Echo();
}
/*}}}*/
void ControlInput::Echo(void){/*{{{*/
	this->DeepEcho();
}
/*}}}*/
void ControlInput::Marshall(MarshallHandle* marshallhandle){ /*{{{*/

	int object_enum = ControlInputEnum;
   marshallhandle->call(object_enum);

   marshallhandle->call(this->control_id);
	marshallhandle->call(this->enum_type);
	marshallhandle->call(this->layout_enum);

	/*Allocate memory*/
	if(marshallhandle->OperationNumber()==MARSHALLING_LOAD){
		switch(this->layout_enum){
			case TriaInputEnum:
				this->values     =new TriaInput();
				this->minvalues  =new TriaInput();
				this->maxvalues  =new TriaInput();
				break;
			case PentaInputEnum:
				this->values     =new PentaInput();
				this->minvalues  =new PentaInput();
				this->maxvalues  =new PentaInput();
				break;
			case TransientInputEnum:
				this->values     =new TransientInput();
				this->minvalues  =new TransientInput();
				this->maxvalues  =new TransientInput();
				break;
			default:
				_error_("Input of Enum \"" << EnumToStringx(this->layout_enum) << "\" not supported yet");
		}
	}

	this->values->Marshall(marshallhandle);
	this->minvalues->Marshall(marshallhandle);
	this->maxvalues->Marshall(marshallhandle);
}
/*}}}*/

void ControlInput::SetControl(int interp,int numindices,int* indices,IssmDouble* values_in,IssmDouble* values_min,IssmDouble* values_max){/*{{{*/

	_assert_(this);

	/*Set input*/
	if(this->values->ObjectEnum()==TriaInputEnum || this->values->ObjectEnum()==PentaInputEnum){
		xDynamicCast<ElementInput*>(this->values)->SetInput(interp,numindices,indices,values_in);
		xDynamicCast<ElementInput*>(this->minvalues)->SetInput(interp,numindices,indices,values_min);
		xDynamicCast<ElementInput*>(this->maxvalues)->SetInput(interp,numindices,indices,values_max);
	}
	else{
		_error_("not supported");
	}
}
/*}}}*/
TriaInput* ControlInput::GetTriaInput(){/*{{{*/

	/*Cast and return*/
	if(this->values->ObjectEnum()==TriaInputEnum){
		return xDynamicCast<TriaInput*>(this->values);
	}
	else if(this->values->ObjectEnum()==TransientInputEnum){
		return xDynamicCast<TransientInput*>(this->values)->GetTriaInput();
	}
	else{
		_error_("Cannot return a TriaInput");
	}

}
/*}}}*/
PentaInput* ControlInput::GetPentaInput(){/*{{{*/

	/*Cast and return*/
	if(this->values->ObjectEnum()!=PentaInputEnum){
		_error_("Cannot return a PentaInput");
	}
	return xDynamicCast<PentaInput*>(this->values);

}
/*}}}*/
ElementInput* ControlInput::GetInput(const char* data){/*{{{*/

	if(strcmp(data,"value")==0){
		_assert_(values);
		return xDynamicCast<ElementInput*>(values);
	}
	else if (strcmp(data,"lowerbound")==0){
		_assert_(minvalues);
		return xDynamicCast<ElementInput*>(minvalues);
	}
	else if (strcmp(data,"upperbound")==0){
		_assert_(maxvalues);
		return xDynamicCast<ElementInput*>(maxvalues);
	}
	else{
		_error_("Data " << data << " not supported yet");
	}
}
/*}}}*/
TransientInput* ControlInput::GetTransientInput(const char* data){/*{{{*/

	if(this->values->ObjectEnum()==TransientInputEnum){
		if(strcmp(data,"value")==0){
			_assert_(values);
			return xDynamicCast<TransientInput*>(values);
		}
		else if (strcmp(data,"lowerbound")==0){
			_assert_(minvalues);
			return xDynamicCast<TransientInput*>(minvalues);
		}
		else if (strcmp(data,"upperbound")==0){
			_assert_(maxvalues);
			return xDynamicCast<TransientInput*>(maxvalues);
		}
		else{
			_error_("Data " << data << " not supported yet");
		}
	}
	else{
		_error_("not supported");
	}

}
/*}}}*/

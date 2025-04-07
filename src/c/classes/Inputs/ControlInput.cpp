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
	savedvalues = NULL;
	minvalues   = NULL;
	maxvalues   = NULL;
	gradient    = NULL;
}
/*}}}*/
ControlInput::ControlInput(int nbe, int nbv,int input_layout_enum,int interp,int id){/*{{{*/

	this->control_id  = id;
	this->layout_enum = input_layout_enum;

	switch(this->layout_enum){
		case TriaInputEnum:
			this->values     =new TriaInput(nbe,nbv,interp);
			this->savedvalues=new TriaInput(nbe,nbv,interp);
			this->minvalues  =new TriaInput(nbe,nbv,interp);
			this->maxvalues  =new TriaInput(nbe,nbv,interp);
			this->gradient   =new TriaInput(nbe,nbv,interp);
			break;
		case PentaInputEnum:
			this->values     =new PentaInput(nbe,nbv,interp);
			this->savedvalues=new PentaInput(nbe,nbv,interp);
			this->minvalues  =new PentaInput(nbe,nbv,interp);
			this->maxvalues  =new PentaInput(nbe,nbv,interp);
			this->gradient   =new PentaInput(nbe,nbv,interp);
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
	this->savedvalues=new TransientInput(enum_in,nbe,nbv,times,numtimes);
	this->minvalues  =new TransientInput(enum_in,nbe,nbv,times,numtimes);
	this->maxvalues  =new TransientInput(enum_in,nbe,nbv,times,numtimes);
	this->gradient   =new TransientInput(enum_in,nbe,nbv,times,numtimes);
}
/*}}}*/
ControlInput::~ControlInput(){/*{{{*/
	delete values;
	delete savedvalues;
	delete minvalues;
	delete maxvalues;
	delete gradient;
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
	if(savedvalues) output->savedvalues = this->savedvalues->copy();
	if(minvalues)   output->minvalues   = this->minvalues->copy();
	if(maxvalues)   output->maxvalues   = this->maxvalues->copy();
	if(gradient)    output->gradient    = this->gradient->copy();

	return output;
}
/*}}}*/
void ControlInput::Configure(Parameters* params){/*{{{*/
	this->values->Configure(params);
	this->savedvalues->Configure(params);
	this->minvalues->Configure(params);
	this->maxvalues->Configure(params);
	this->gradient->Configure(params);
}
/*}}}*/
void ControlInput::DeepEcho(void){/*{{{*/

	_printf_("ControlInput:\n");
	_printf_(setw(15)<<"   ControlInput "<<setw(25)<<left<<EnumToStringx(this->enum_type)<<"\n");
	_printf_(setw(15)<<"   Layout       "<<setw(25)<<left<<EnumToStringx(this->layout_enum)<<"\n");
	_printf_("---values: \n");     if (values)      values->Echo();
	_printf_("---savedvalues: \n");if (savedvalues) savedvalues->Echo();
	_printf_("---minvalues: \n");  if (minvalues)   minvalues->Echo();
	_printf_("---maxvalues: \n");  if (maxvalues)   maxvalues->Echo();
	_printf_("---gradient: \n");   if (gradient){    gradient->Echo();} else{_printf_("     Not set yet\n");}
}
/*}}}*/
void ControlInput::Echo(void){/*{{{*/
	this->DeepEcho();
}
/*}}}*/
int  ControlInput::Id(void){ return -1; }/*{{{*/
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
				this->savedvalues=new TriaInput();
				this->minvalues  =new TriaInput();
				this->maxvalues  =new TriaInput();
				this->gradient   =new TriaInput();
				break;
			case PentaInputEnum:
				this->values     =new PentaInput();
				this->savedvalues=new PentaInput();
				this->minvalues  =new PentaInput();
				this->maxvalues  =new PentaInput();
				this->gradient   =new PentaInput();
				break;
			case TransientInputEnum:
				this->values     =new TransientInput();
				this->savedvalues=new TransientInput();
				this->minvalues  =new TransientInput();
				this->maxvalues  =new TransientInput();
				this->gradient   =new TransientInput();
				break;
			default:
				_error_("Input of Enum \"" << EnumToStringx(this->layout_enum) << "\" not supported yet");
		}
	}

	this->values->Marshall(marshallhandle);
	this->savedvalues->Marshall(marshallhandle);
	this->minvalues->Marshall(marshallhandle);
	this->maxvalues->Marshall(marshallhandle);
	this->gradient->Marshall(marshallhandle);
}
/*}}}*/
int  ControlInput::ObjectEnum(void){/*{{{*/

	return ControlInputEnum;

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
void ControlInput::SetGradient(int interp,int numindices,int* indices,IssmDouble* values_in){/*{{{*/

	_assert_(this);
	_assert_(this->gradient);
	if(this->gradient->ObjectEnum()==TriaInputEnum || this->gradient->ObjectEnum()==PentaInputEnum){
		xDynamicCast<ElementInput*>(this->gradient)->SetInput(interp,numindices,indices,values_in);
	}
	else{
		_error_("not supported");
	}
}
/*}}}*/
void ControlInput::SetGradient(int interp,int numindices,int* indices,IssmDouble* values_in,int n){/*{{{*/

	if(this->values->ObjectEnum()!=TransientInputEnum)_error_("you are in the wrong place, go home");
	_assert_(this);
	_assert_(this->gradient);
	_error_("S");

	//NEW??
	//this->gradient->SetInput(interp,numindices,indices,values_in);
}
/*}}}*/
void ControlInput::AverageAndReplace(void){/*{{{*/
	this->values->AverageAndReplace();
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

	if(this->gradient->ObjectEnum()==TriaInputEnum || this->gradient->ObjectEnum()==PentaInputEnum){
		if(strcmp(data,"value")==0){
			_assert_(values);
			return xDynamicCast<ElementInput*>(values);
		}
		else if(strcmp(data,"savedvalues")==0){
			_assert_(savedvalues);
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
		else if (strcmp(data,"gradient")==0){
			_assert_(gradient);
			return xDynamicCast<ElementInput*>(gradient);
		}
		else{
			_error_("Data " << data << " not supported yet");
		}
	}
	else{
		int* temp = xNew<int>(3);
		_error_("Gradient is of type "<<EnumToStringx(this->gradient->ObjectEnum()) <<", which is not supported yet");
	}

}
/*}}}*/
TransientInput* ControlInput::GetTransientInput(const char* data){/*{{{*/

	if(this->values->ObjectEnum()==TransientInputEnum){
		if(strcmp(data,"value")==0){
			_assert_(values);
			return xDynamicCast<TransientInput*>(values);
		}
		else if(strcmp(data,"savedvalues")==0){
			_assert_(savedvalues);
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
		else if (strcmp(data,"gradient")==0){
			_assert_(gradient);
			return xDynamicCast<TransientInput*>(gradient);
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

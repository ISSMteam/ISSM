/*!\file TransientFileInput.c
 * \brief: implementation of the TransientFileInput object
 */
/*Headers*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include <numeric>
#include "./TransientFileInput.h"
#include "./TriaInput.h"
#include "./PentaInput.h"
#include "../../shared/shared.h"
#include "../Params/Parameters.h"

/*TransientFileInput constructors and destructor*/
TransientFileInput::TransientFileInput(){/*{{{*/

	this->enum_type=UNDEF;
	this->inputs=NULL;
	this->numtimesteps=0;
	this->parameters=NULL;
	this->timesteps=NULL;

	this->current_input=NULL;
	this->current_step=-1;

}
/*}}}*/
TransientFileInput::TransientFileInput(int in_enum_type,int nbe,int nbv,char* filename_in,IssmPDouble period_in){/*{{{*/

	/*Set Enum*/
	this->enum_type=in_enum_type;
	this->numberofelements_local = nbe;
	this->numberofvertices_local = nbv;
   this->filename= xNew<char>(strlen(filename_in)+1);
   xMemCpy<char>(this->filename,filename_in,(strlen(filename_in)+1));
   this->loading_period = period_in; _assert_(period_in>0.);

	/*Set current inputs to empty for now*/
	this->numtimesteps = 0;
	this->timesteps    = NULL;
	this->inputs       = NULL;

	this->parameters    = NULL;
	this->current_input = NULL;
	this->current_step  = -1;
}
/*}}}*/
TransientFileInput::~TransientFileInput(){/*{{{*/

	for(int i=0;i<this->numtimesteps;i++){
		delete this->inputs[i];
	}
	xDelete<Input*>(this->inputs);
	xDelete<IssmDouble>(this->timesteps);
   xDelete<char>(this->filename);

	if(this->current_input) delete this->current_input;
}
/*}}}*/

/*Object virtual functions definitions:*/
Input* TransientFileInput::copy() {/*{{{*/

	TransientFileInput* output=NULL;

	output = new TransientFileInput();
	output->enum_type=this->enum_type;
	output->numtimesteps=this->numtimesteps;
	if(this->numtimesteps>0){
		output->timesteps=xNew<IssmDouble>(this->numtimesteps);
		xMemCpy(output->timesteps,this->timesteps,this->numtimesteps);
		output->inputs = xNew<Input*>(this->numtimesteps);
		for(int i=0;i<this->numtimesteps;i++){
			if(this->inputs[i]){
				output->inputs[i] = this->inputs[i]->copy();
			}
			else{
				output->inputs[i] = NULL;
			}
		}
	}
	output->parameters=this->parameters;

   output->filename= xNew<char>(strlen(this->filename)+1);
   xMemCpy<char>(output->filename,this->filename,(strlen(this->filename)+1));
   output->loading_period = this->loading_period;

	return output;
}/*}}}*/
void TransientFileInput::DeepEcho(void){/*{{{*/

	_printf_("TransientFileInput:\n");
	_printf_("   enum: " << this->enum_type << " (" << EnumToStringx(this->enum_type) << ")\n");
   _printf_("   filename:       " << this->filename <<"\n");
   _printf_("   loading_period: " << this->loading_period<<"\n");
	_printf_("   numtimesteps: " << this->numtimesteps << "\n");
	_printf_("---inputs: \n");
	for(int i=0;i<this->numtimesteps;i++){
		_printf_("   time: " << this->timesteps[i]<<"  ");
		if(this->inputs[i]) this->inputs[i]->DeepEcho();
		else                _printf_(" NOT SET! \n");
	}
}
/*}}}*/
void TransientFileInput::Configure(Parameters* params){/*{{{*/
	this->parameters=params;
}
/*}}}*/
void TransientFileInput::Echo(void){/*{{{*/
	_printf_("TransientFileInput:\n");
   _printf_("   enum: " << this->enum_type << " (" << EnumToStringx(this->enum_type) << ")\n");
   _printf_("   filename:       " << this->filename <<"\n");
   _printf_("   loading_period: " << this->loading_period<<"\n");
	_printf_("   numtimesteps: " << this->numtimesteps << "\n");
	_printf_("---inputs: \n");
	for(int i=0;i<this->numtimesteps;i++){
		_printf_("   time: " << this->timesteps[i]<<"  ");
		if(this->inputs[i]) this->inputs[i]->Echo();
		else                _printf_(" NOT SET! \n");
	}
}
/*}}}*/
int  TransientFileInput::Id(void){ return -1; }/*{{{*/
/*}}}*/
void TransientFileInput::Marshall(MarshallHandle* marshallhandle){ /*{{{*/

	bool isnull;

	int object_enum = TransientFileInputEnum;
   marshallhandle->call(object_enum);

	marshallhandle->call(this->numberofelements_local);
	marshallhandle->call(this->numberofvertices_local);
	marshallhandle->call(this->enum_type);
   marshallhandle->call(this->filename);
   marshallhandle->call(this->loading_period);
	marshallhandle->call(this->numtimesteps);
	marshallhandle->call(this->timesteps,numtimesteps);

	/*Allocate memory if need be*/
	if(marshallhandle->OperationNumber() == MARSHALLING_LOAD){
		int N = this->numtimesteps; _assert_(N>=0 && N<1e6);
		if(N){
			this->inputs = xNew<Input*>(N);
			for(int i=0;i<N;i++) this->inputs[i] = NULL;
		}
		else{
			this->inputs = NULL;
		}
	}

	/*Marshall!*/
	if(marshallhandle->OperationNumber()!=MARSHALLING_LOAD){
		for(int i=0;i<this->numtimesteps;i++){

			//_assert_(this->inputs[i]);
			isnull = false;
			if(!this->inputs[i]) isnull = true;
			marshallhandle->call(isnull);

			if(!isnull){
				object_enum = this->inputs[i]->ObjectEnum();
				marshallhandle->call(object_enum);
				this->inputs[i]->Marshall(marshallhandle);
			}
		}
	}
	else{
		for(int i=0;i<this->numtimesteps;i++){
			marshallhandle->call(isnull);
			if(!isnull){
				marshallhandle->call(object_enum);

				if(object_enum==TriaInputEnum){
					TriaInput* triainput2=new TriaInput();
					triainput2->Marshall(marshallhandle);
					this->inputs[i]=triainput2;
				}
				else if(object_enum==PentaInputEnum){
					PentaInput* pentainput2=new PentaInput();
					pentainput2->Marshall(marshallhandle);
					this->inputs[i]=pentainput2;
				}
				else{
					_error_("input "<<EnumToStringx(object_enum)<<" not supported");
				}
			}
		}
	}

}
/*}}}*/
int  TransientFileInput::ObjectEnum(void){/*{{{*/

	return TransientFileInputEnum;

}
/*}}}*/

/*Intermediary*/
bool TransientFileInput::IsFileInputUpdate(IssmDouble time){/*{{{*/


	/*Do we have anything loaded right now?*/
	if(this->numtimesteps==0) return true;

	/*Check if we are less than */
	_assert_(this->timesteps);
	_assert_(this->loading_period>0.);
	if(time - this->timesteps[this->numtimesteps-1] < this->loading_period){
		return true;
	}

	return false;
}/*}}}*/

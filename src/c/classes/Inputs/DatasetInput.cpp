/*!\file DatasetInput.c
 * \brief: implementation of the datasetinput object
 */
/*Headers*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./DatasetInput.h"
#include "./TriaInput.h"
#include "./PentaInput.h"
#include "./TransientInput.h"

/*DatasetInput constructors and destructor*/
DatasetInput::DatasetInput(){/*{{{*/
	this->inputs    = NULL;
	this->numids    = 0;
	this->ids       = NULL;
	this->numberofelements_local = -1;
	this->numberofvertices_local = -1;
}
/*}}}*/
DatasetInput::DatasetInput(int nbe, int nbv){/*{{{*/
	this->inputs    = NULL;
	this->numids    = 0;
	this->ids       = NULL;
	this->numberofelements_local = nbe;
	this->numberofvertices_local = nbv;
}
/*}}}*/
DatasetInput::~DatasetInput(){/*{{{*/
	xDelete<int>(this->ids);
	for(int i=0;i<this->numids;i++){
		delete this->inputs[i];
	}
	xDelete<Input*>(this->inputs);
}
/*}}}*/

/*Object virtual functions definitions:*/
Input* DatasetInput::copy() {/*{{{*/

	DatasetInput* output=NULL;

	output = new DatasetInput();
	output->numids=this->numids;
	if(this->numids>0){
		output->ids=xNew<int>(output->numids);
		xMemCpy(output->ids,this->ids,output->numids);
		output->inputs = xNew<Input*>(this->numids);
		for(int i=0;i<this->numids;i++){
			output->inputs[i] = this->inputs[i]->copy();
		}
	}

	return output;
}
/*}}}*/
void DatasetInput::Configure(Parameters* params){/*{{{*/
	for(int i=0;i<this->numids;i++){
		this->inputs[i]->Configure(params);
	}
}
/*}}}*/
void DatasetInput::DeepEcho(void){/*{{{*/

	_printf_("DatasetInput:\n");
	_printf_("   numids:"<< this->numids<< "\n");
	_printf_("      ids: ");
	for(int i=0;i<this->numids;i++) _printf_(this->ids[i]<<" ("<<EnumToStringx(this->ids[i])<<") ");
	_printf_("\n");
	//_printf_("   inputs: \n"); inputs->Echo();
}
/*}}}*/
void DatasetInput::Echo(void){/*{{{*/
	this->DeepEcho();
}
/*}}}*/
int  DatasetInput::Id(void){ return -1; }/*{{{*/
/*}}}*/
void DatasetInput::Marshall(MarshallHandle* marshallhandle){ /*{{{*/

	int object_enum = DatasetInputEnum;
   marshallhandle->call(object_enum);

	marshallhandle->call(this->numids);
	marshallhandle->call(this->numberofelements_local);
	marshallhandle->call(this->numberofvertices_local);
	marshallhandle->call(this->ids,numids);

	/*Allocate memory if need be*/
	int N = this->numids; _assert_(N>=0 && N<1e6);
	if(marshallhandle->OperationNumber() == MARSHALLING_LOAD){
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
		for(int i=0;i<N;i++){
			_assert_(this->inputs[i]);
			object_enum = this->inputs[i]->ObjectEnum();
			marshallhandle->call(object_enum);
			this->inputs[i]->Marshall(marshallhandle);
		}
	}
	else{
		for(int i=0;i<N;i++){
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
			else if(object_enum==TransientInputEnum){
				TransientInput* transinput2=new TransientInput();
				transinput2->Marshall(marshallhandle);
				this->inputs[i]=transinput2;
			}
			else{
				_error_("input "<<EnumToStringx(object_enum)<<" not supported");
			}
		}
	}

}
/*}}}*/
int  DatasetInput::ObjectEnum(void){/*{{{*/
	return DatasetInputEnum;
}/*}}}*/

void DatasetInput::SetTriaInput(int id,int interp_in,int numinds,int* rows,IssmDouble* values_in){ /*{{{*/

	int  index = -1;
	for(int i=0;i<this->numids;i++){
		if(this->ids[i] == id) index = i;
	}

	/*Create new input if not found*/
	if(index == -1){
		int* new_ids = xNew<int>(this->numids+1);
		if(this->numids) xMemCpy(new_ids,this->ids,this->numids);
		new_ids[this->numids] = id;

		Input** new_inputs = xNew<Input*>(this->numids+1);
		if(this->numids) xMemCpy(new_inputs,this->inputs,this->numids);
		new_inputs[this->numids] = new TriaInput(this->numberofelements_local,this->numberofvertices_local,interp_in);
		index = this->numids;

		xDelete<int>(this->ids);
		this->ids = new_ids;
		xDelete<Input*>(this->inputs);
		this->inputs = new_inputs;

		this->numids ++;
	}

	/*Set input*/
	if(this->inputs[index]->ObjectEnum()!=TriaInputEnum) _error_("cannot add Element values to a "<<EnumToStringx(this->inputs[index]->ObjectEnum()));
	TriaInput* input = xDynamicCast<TriaInput*>(this->inputs[index]);
	input->SetInput(interp_in,numinds,rows,values_in);

}
/*}}}*/
void DatasetInput::SetPentaInput(int id,int interp_in,int numinds,int* rows,IssmDouble* values_in){ /*{{{*/

	int  index = -1;
	for(int i=0;i<this->numids;i++){
		if(this->ids[i] == id) index = i;
	}

	/*Create new input if not found*/
	if(index == -1){
		int* new_ids = xNew<int>(this->numids+1);
		if(this->numids) xMemCpy(new_ids,this->ids,this->numids);
		new_ids[this->numids] = id;

		Input** new_inputs = xNew<Input*>(this->numids+1);
		if(this->numids) xMemCpy(new_inputs,this->inputs,this->numids);
		new_inputs[this->numids] = new PentaInput(this->numberofelements_local,this->numberofvertices_local,interp_in);
		index = this->numids;

		xDelete<int>(this->ids);
		this->ids = new_ids;
		xDelete<Input*>(this->inputs);
		this->inputs = new_inputs;

		this->numids ++;
	}

	/*Set input*/
	if(this->inputs[index]->ObjectEnum()!=PentaInputEnum) _error_("cannot add Element values to a "<<EnumToStringx(this->inputs[index]->ObjectEnum()));
	PentaInput* input = xDynamicCast<PentaInput*>(this->inputs[index]);
	input->SetInput(interp_in,numinds,rows,values_in);

}
/*}}}*/
TransientInput* DatasetInput::SetTransientInput(int id,IssmDouble* times,int numtimes){ /*{{{*/

	int  index = -1;
	for(int i=0;i<this->numids;i++){
		if(this->ids[i] == id) index = i;
	}

	/*Create new input if not found*/
	if(index == -1){
		int* new_ids = xNew<int>(this->numids+1);
		if(this->numids) xMemCpy(new_ids,this->ids,this->numids);
		new_ids[this->numids] = id;

		Input** new_inputs = xNew<Input*>(this->numids+1);
		if(this->numids) xMemCpy(new_inputs,this->inputs,this->numids);
		new_inputs[this->numids] = new TransientInput(NoneEnum,this->numberofelements_local,this->numberofvertices_local,times,numtimes);
		index = this->numids;

		xDelete<int>(this->ids);
		this->ids = new_ids;
		xDelete<Input*>(this->inputs);
		this->inputs = new_inputs;

		this->numids ++;
	}

	/*Set input*/
	if(this->inputs[index]->ObjectEnum()!=TransientInputEnum) _error_("cannot add values to a "<<EnumToStringx(this->inputs[index]->ObjectEnum()));
	TransientInput* input = xDynamicCast<TransientInput*>(this->inputs[index]);
	return input;
}
/*}}}*/
void DatasetInput::GetInputValue(IssmDouble* pvalue,Gauss* gauss,int id){ /*{{{*/

	int  index = -1;
	for(int i=0;i<this->numids;i++){
		if(this->ids[i] == id) index = i;
	}

	/*Create new input if not found*/
	if(index == -1){
		this->Echo();
		_error_("Could not find input "<<id<<" ("<<EnumToStringx(id)<<"?) in DatasetInput");
	}

	Input* input = this->inputs[index];

	if(this->inputs[index]->ObjectEnum()==TransientInputEnum){
		input = xDynamicCast<TransientInput*>(this->inputs[index])->current_input;
	}

	input->GetInputValue(pvalue,gauss);

}
/*}}}*/
IssmDouble DatasetInput::GetInputMin(void){ /*{{{*/

	IssmDouble minvalue,newminvalue;
	for(int i=0;i<this->numids;i++){

		Input* input = this->inputs[i];

		if(this->inputs[i]->ObjectEnum()==TransientInputEnum){
			input = xDynamicCast<TransientInput*>(this->inputs[i])->current_input;
		}
		newminvalue=input->GetInputMin();
		if(i==0)minvalue=newminvalue;
		else minvalue=min(minvalue,newminvalue);
	}
	return minvalue;

}
/*}}}*/
TransientInput* DatasetInput::GetTransientInputByOffset(int offset){/*{{{*/

	_assert_(offset>=0 && offset<this->numids);
	_assert_(this->inputs[offset]);

	/*Cast and return*/
	if(this->inputs[offset]->ObjectEnum()==TransientInputEnum){
		return xDynamicCast<TransientInput*>(this->inputs[offset]);
	}
	else{
		_error_("Cannot return a TransientInput");
	}
}/*}}}*/
TriaInput* DatasetInput::GetTriaInput(void){/*{{{*/

	return this->GetTriaInputByOffset(0);

}/*}}}*/
TriaInput* DatasetInput::GetTriaInputByOffset(int offset){/*{{{*/

	_assert_(offset>=0 && offset<this->numids);
	_assert_(this->inputs[offset]);

	/*Cast and return*/
	if(this->inputs[offset]->ObjectEnum()==TransientInputEnum){
		return xDynamicCast<TransientInput*>(this->inputs[offset])->GetTriaInput();
	}
	if(this->inputs[offset]->ObjectEnum()!=TriaInputEnum){
		_error_("Cannot return a TriaInput");
	}
	return xDynamicCast<TriaInput*>(this->inputs[offset]);

}/*}}}*/
PentaInput* DatasetInput::GetPentaInputByOffset(int offset){/*{{{*/

	_assert_(offset>=0 && offset<this->numids);
	_assert_(this->inputs[offset]);

	/*Cast and return*/
	if(this->inputs[offset]->ObjectEnum()==TransientInputEnum){
		return xDynamicCast<TransientInput*>(this->inputs[offset])->GetPentaInput();
	}
	if(this->inputs[offset]->ObjectEnum()!=PentaInputEnum){
		_error_("Cannot return a PentaInput");
	}
	return xDynamicCast<PentaInput*>(this->inputs[offset]);

}/*}}}*/

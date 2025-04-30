/*\file Inputs.cpp
 * \brief: Implementation of the Inputs class, derived from DataSet class.
 */

/*Headers: {{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./Input.h"
#include "./Inputs.h"

#include "./BoolInput.h"
#include "./IntInput.h"
#include "./DoubleInput.h"
#include "./ElementInput.h"
#include "./SegInput.h"
#include "./TriaInput.h"
#include "./PentaInput.h"
#include "./TransientInput.h"
#include "./TransientFileInput.h"
#include "./ControlInput.h"
#include "./DatasetInput.h"
#include "./ArrayInput.h"
#include "./IntArrayInput.h"
using namespace std;
/*}}}*/

/*Object constructors and destructor*/
Inputs::Inputs(void){/*{{{*/

	this->numberofelements_local = 0;
	this->numberofvertices_local = 0;

	/*Initialize pointers*/
	for(int i=0;i<NUMINPUTS;i++) this->inputs[i] = NULL;
}
/*}}}*/
Inputs::Inputs(int nbe,int nbv){/*{{{*/

	this->numberofelements_local = nbe;
	this->numberofvertices_local = nbv;

	/*Initialize pointers*/
	for(int i=0;i<NUMINPUTS;i++) this->inputs[i] = NULL;
}
/*}}}*/
Inputs::~Inputs(){/*{{{*/
	for(int i=0;i<NUMINPUTS;i++){
		if(this->inputs[i]) delete this->inputs[i];
	}
	return;
}
/*}}}*/

Inputs* Inputs::Copy(void){/*{{{*/

	Inputs* output = new Inputs(this->numberofelements_local,this->numberofvertices_local);

	for(int i=0;i<NUMINPUTS;i++){
		if(this->inputs[i]) output->inputs[i]=this->inputs[i]->copy();
	}

	return output;
}/*}}}*/
void Inputs::DeepEcho(void){/*{{{*/
	for(int i=0;i<NUMINPUTS;i++) {
		if(this->inputs[i]) this->inputs[i]->DeepEcho();
	}
	return;
}
/*}}}*/
void Inputs::Echo(void){/*{{{*/
	_printf_("Inputs Echo:\n");
	for(int i=0;i<NUMINPUTS;i++) {
		if(this->inputs[i]) _printf_(setw(25)<<EnumToStringx(i+InputsSTARTEnum+1)<<": set as "<<EnumToStringx(this->inputs[i]->ObjectEnum())<<"\n");
	}
	return;
}
/*}}}*/
void Inputs::DeepEcho(int enum_in){/*{{{*/
	int index= EnumToIndex(enum_in);
	if(this->inputs[index])this->inputs[index]->DeepEcho();
	return;
}
/*}}}*/
void Inputs::Echo(int enum_in){/*{{{*/
	_printf_("Inputs Echo:\n");
	int index= EnumToIndex(enum_in);
	if(this->inputs[index])this->inputs[index]->Echo();
	return;
}
/*}}}*/
void Inputs::Marshall(MarshallHandle* marshallhandle){/*{{{*/

	int num_inputs=0;
	int index;

	int object_enum = InputsEnum;
   marshallhandle->call(object_enum);
	marshallhandle->call(this->numberofelements_local);
	marshallhandle->call(this->numberofvertices_local);

	if(marshallhandle->OperationNumber()!=MARSHALLING_LOAD){

		/*Marshall num_inputs first*/
		for(int i=0;i<NUMINPUTS;i++){
			if(this->inputs[i]) num_inputs++;
		}
		marshallhandle->call(num_inputs);

		/*Marshall Parameters one by one now*/
		for(int i=0;i<NUMINPUTS;i++){
			if(this->inputs[i]){
				object_enum = this->inputs[i]->ObjectEnum();
				marshallhandle->call(i);
				marshallhandle->call(object_enum);
				this->inputs[i]->Marshall(marshallhandle);
			}
		}
	}
	else{

		/*Get number of inputs marshalled*/
		marshallhandle->call(num_inputs);

		/*Recover input2eters one by one*/
		for(int i=0;i<num_inputs;i++){

			/*Recover enum of object first: */
			marshallhandle->call(index);
			marshallhandle->call(object_enum);

			if(object_enum==BoolInputEnum){
				BoolInput* boolinput2=new BoolInput();
				boolinput2->Marshall(marshallhandle);
				this->inputs[index]=boolinput2;
			}
			else if(object_enum==IntInputEnum){
				IntInput* intinput2=new IntInput();
				intinput2->Marshall(marshallhandle);
				this->inputs[index]=intinput2;
			}
			else if(object_enum==TriaInputEnum){
				TriaInput* triainput2=new TriaInput();
				triainput2->Marshall(marshallhandle);
				this->inputs[index]=triainput2;
			}
			else if(object_enum==PentaInputEnum){
				PentaInput* pentainput2=new PentaInput();
				pentainput2->Marshall(marshallhandle);
				this->inputs[index]=pentainput2;
			}
			else if(object_enum==ControlInputEnum){
				ControlInput* input=new ControlInput();
				input->Marshall(marshallhandle);
				this->inputs[index]=input;
			}
			else if(object_enum==TransientInputEnum){
				TransientInput* input=new TransientInput();
				input->Marshall(marshallhandle);
				this->inputs[index]=input;
			}
			else if(object_enum==DatasetInputEnum){
				DatasetInput* input=new DatasetInput();
				input->Marshall(marshallhandle);
				this->inputs[index]=input;
			}
			else if(object_enum==ArrayInputEnum){
				ArrayInput* input=new ArrayInput();
				input->Marshall(marshallhandle);
				this->inputs[index]=input;
			}
			else{
				_error_("input "<<EnumToStringx(object_enum)<<" not supported");
			}
		}
	}
}
/*}}}*/

void Inputs::AddInput(Input* newinput){/*{{{*/

	/*Get Enum from Param*/
	_assert_(newinput);
	int input_enum = newinput->InstanceEnum();

	/*Get index in array*/
	#ifdef _ISSM_DEBUG_
	if(input_enum<=InputsSTARTEnum) _error_("Cannot add input: Enum "<<EnumToStringx(input_enum)<<" should appear after InputsSTARTEnum");
	if(input_enum>=InputsENDEnum)   _error_("Cannot add input: Enum "<<EnumToStringx(input_enum)<<" should appear before InputsENDEnum");
	#endif
	int index = input_enum - InputsSTARTEnum -1;

	/*Delete input if it already exists*/
	if(this->inputs[index]){
		delete this->inputs[index];
		this->inputs[index] = NULL;
	}

	/*Add input to array*/
	this->inputs[index] = newinput;
}
/*}}}*/
void Inputs::ChangeEnum(int oldenumtype,int newenumtype){/*{{{*/

	/*Get indices from enums*/
	int index_old = EnumToIndex(oldenumtype);
	int index_new = EnumToIndex(newenumtype);

	/*Delete input if it already exists*/
	if(this->inputs[index_new]) delete this->inputs[index_new];

	/*Make sure that old one exists*/
	if(!this->inputs[index_old]){
		_error_("Input "<<EnumToStringx(oldenumtype)<<" not found");
	}

	/*Replace Enums*/
	this->inputs[index_old]->ChangeEnum(newenumtype);
	this->inputs[index_new] = this->inputs[index_old];
	this->inputs[index_old] = NULL;
}/*}}}*/
void Inputs::Configure(Parameters* parameters){/*{{{*/
	for(int i=0;i<NUMINPUTS;i++){
		if(this->inputs[i]) this->inputs[i]->Configure(parameters);
	}
}
/*}}}*/
int  Inputs::DeleteInput(int input_enum){/*{{{*/

	int index = EnumToIndex(input_enum);
	if(this->inputs[index]){
		delete this->inputs[index];
		this->inputs[index] = NULL;
	}

	return 1;
}
/*}}}*/
void Inputs::DuplicateInput(int original_enum,int new_enum){/*{{{*/

	_assert_(this);

	/*Get indices from enums*/
	int index_ori = EnumToIndex(original_enum);
	int index_new = EnumToIndex(new_enum);

	/*Delete input if it already exists*/
	if(this->inputs[index_new]) delete this->inputs[index_new];

	/*Make sure that old one exists*/
	if(!this->inputs[index_ori]){
		_error_("Input "<<EnumToStringx(original_enum)<<" not found");
	}

	/*Make a copy*/
	Input* copy=this->inputs[index_ori]->copy();

	/*Add copy*/
	this->inputs[index_new] = copy;
}
/*}}}*/
void Inputs::ZAXPY(IssmDouble alpha, int xenum, int yenum, int zenum){/*{{{*/

	_assert_(this);

	/*Get indices from enums*/
	int index_x = EnumToIndex(xenum);
	int index_y = EnumToIndex(yenum);
	int index_z = EnumToIndex(zenum);

	/*Delete output if it already exists*/
	if(this->inputs[index_z]) delete this->inputs[index_z];

	/*Make sure that old one exists*/
	if(!this->inputs[index_x]) _error_("Input "<<EnumToStringx(xenum)<<" not found");
	if(!this->inputs[index_y]) _error_("Input "<<EnumToStringx(yenum)<<" not found");

	/*Make a copy*/
	this->inputs[index_z]=this->inputs[index_y]->copy();

	/*AXPY: */
	this->inputs[index_z]->AXPY(this->inputs[index_x],alpha);
}
/*}}}*/
void Inputs::AXPY(IssmDouble alpha, int xenum, int yenum ){/*{{{*/

	_assert_(this);

	/*Get indices from enums*/
	int index_x = EnumToIndex(xenum);
	int index_y = EnumToIndex(yenum);

	/*Make sure that old one exists*/
	if(!this->inputs[index_x]) _error_("Input "<<EnumToStringx(xenum)<<" not found");
	if(!this->inputs[index_y]) _error_("Input "<<EnumToStringx(yenum)<<" not found");

	/*AXPY: */
	this->inputs[index_y]->AXPY(this->inputs[index_x],alpha);
}
/*}}}*/
void Inputs::Shift(int xenum, IssmDouble alpha){/*{{{*/

	_assert_(this);

	/*Get indices from enums*/
	int index_x = EnumToIndex(xenum);

	/*Make sure that x exists*/
	if(!this->inputs[index_x]) _error_("Input "<<EnumToStringx(xenum)<<" not found");

	/*Shift: */
	this->inputs[index_x]->Shift(alpha);
}
/*}}}*/
void Inputs::AverageAndReplace(int inputenum){/*{{{*/

	_assert_(this);

	/*Get indices from enums*/
	int index = EnumToIndex(inputenum);
	if(!this->inputs[index]) _error_("Input "<<EnumToStringx(inputenum)<<" not found");

	this->inputs[index]->AverageAndReplace();
}
/*}}}*/
int  Inputs::EnumToIndex(int enum_in){/*{{{*/

	_assert_(this);

	/*Make sure this parameter is at the right place*/
	#ifdef _ISSM_DEBUG_
	if(enum_in<=InputsSTARTEnum){
		//int* temp = xNew<int>(3);
		_error_("Enum "<<EnumToStringx(enum_in)<<" should appear after InputsSTARTEnum");
	}
	if(enum_in>=InputsENDEnum){
		_error_("Enum "<<EnumToStringx(enum_in)<<" should appear before InputsENDEnum");
	}
	#endif
	return enum_in - InputsSTARTEnum -1;
}/*}}}*/
bool Inputs::Exist(int enum_in){/*{{{*/

	_assert_(this);

	int index = EnumToIndex(enum_in);
	if(this->inputs[index]) return true;
	return false;
}
/*}}}*/
int Inputs::GetInputObjectEnum(int enum_in){/*{{{*/

	_assert_(this);

	int index = EnumToIndex(enum_in);
	if(!this->inputs[index]) _error_("Input "<<EnumToStringx(enum_in)<<" not found");
	return this->inputs[index]->ObjectEnum();
}
/*}}}*/
void Inputs::GetInputsInterpolations(int* pnuminputs,int** pinterpolations,int** pinputenums){/*{{{*/

	/*First count number of inputs*/
	int count = 0;
	for(int i=0;i<NUMINPUTS;i++){
		if(this->inputs[i]) count++;
	}
	int numinputs = count;

	/*Allocate output*/
	int* interpolations = xNew<int>(count);
	int* enumlist       = xNew<int>(count);

	/*Go through all inputs and assign interpolation in vector*/
	count = 0;
	for(int i=0;i<NUMINPUTS;i++){

		Input* input=this->inputs[i];
		if(!input) continue;

		enumlist[count] = i+InputsSTARTEnum+1;
		switch(input->ObjectEnum()){
			case BoolInputEnum:
			case IntInputEnum:
				interpolations[count] = input->ObjectEnum();
				break;
			case TriaInputEnum:
				interpolations[count] = input->GetResultInterpolation();
				break;
			default:
				_error_("Input "<<EnumToStringx(input->ObjectEnum())<<" not supported yet");
		}
		count++;
	}
	_assert_(count == numinputs);

	/*Return pointer*/
	*pnuminputs = numinputs;
	*pinterpolations = interpolations;
	*pinputenums = enumlist;

}/*}}}*/
SegInput* Inputs::GetSegInput(int enum_in){/*{{{*/

	_assert_(this);

	/*Get input id*/
	int id = EnumToIndex(enum_in);

	/*Check that it has the right format*/
	Input* input = this->inputs[id];
	if(!input) return NULL;

	return input->GetSegInput();
}/*}}}*/
TriaInput* Inputs::GetTriaInput(int enum_in){/*{{{*/

	_assert_(this);

	/*Get input id*/
	int id = EnumToIndex(enum_in);

	/*Check that it has the right format*/
	Input* input = this->inputs[id];
	if(!input) return NULL;

	return input->GetTriaInput();
}/*}}}*/
TriaInput* Inputs::GetTriaInput(int enum_in,IssmDouble time){/*{{{*/

	/*Get input id*/
	int id = EnumToIndex(enum_in);

	/*Check that it has the right format*/
	Input* input = this->inputs[id];
	if(!input) return NULL;

	if(input->ObjectEnum()==TransientInputEnum){
		return xDynamicCast<TransientInput*>(input)->GetTriaInput(time);
	}
	else{
		return input->GetTriaInput();
	}
}/*}}}*/
TriaInput* Inputs::GetTriaInput(int enum_in,IssmDouble start_time,IssmDouble end_time,int averaging_method){/*{{{*/

	/*Get input id*/
	int id = EnumToIndex(enum_in);

	/*Check that it has the right format*/
	Input* input = this->inputs[id];
	if(!input) return NULL;

	if(input->ObjectEnum()==TransientInputEnum){
		return xDynamicCast<TransientInput*>(input)->GetTriaInput(start_time,end_time,averaging_method);
	}
	else{
		_error_("Input "<<EnumToStringx(enum_in)<<" is not an TransientInput");
	}
}/*}}}*/
PentaInput* Inputs::GetPentaInput(int enum_in){/*{{{*/

	/*Get input id*/
	int id = EnumToIndex(enum_in);

	/*Check that it has the right format*/
	Input* input = this->inputs[id];
	if(!input) return NULL;

	return input->GetPentaInput();
}/*}}}*/
PentaInput* Inputs::GetPentaInput(int enum_in,IssmDouble time){/*{{{*/

	/*Get input id*/
	int id = EnumToIndex(enum_in);

	/*Check that it has the right format*/
	Input* input = this->inputs[id];
	if(!input) return NULL;

	if(input->ObjectEnum()==TransientInputEnum){
		return xDynamicCast<TransientInput*>(input)->GetPentaInput(time);
	}
	else{
		return input->GetPentaInput();
	}
}/*}}}*/
PentaInput* Inputs::GetPentaInput(int enum_in,IssmDouble start_time,IssmDouble end_time,int averaging_method){/*{{{*/

	/*Get input id*/
	int id = EnumToIndex(enum_in);

	/*Check that it has the right format*/
	Input* input = this->inputs[id];
	if(!input) return NULL;

	if(input->ObjectEnum()==TransientInputEnum){
		return xDynamicCast<TransientInput*>(input)->GetPentaInput(start_time,end_time,averaging_method);
	}
	else{
		_error_("Input "<<EnumToStringx(enum_in)<<" is not an TransientInput");
	}
}/*}}}*/
TransientInput* Inputs::GetTransientInput(int enum_in){/*{{{*/

	/*Get input id*/
	int id = EnumToIndex(enum_in);

	/*Check that it has the right format*/
	Input* input = this->inputs[id];
	if(!input) return NULL;

	if(input->ObjectEnum() != TransientInputEnum){
		_error_("Input "<<EnumToStringx(enum_in)<<" is not an TransientInput");
	}

	/*Cast and return*/
	TransientInput* output = xDynamicCast<TransientInput*>(input);
	return output;
}/*}}}*/
ElementInput* Inputs::GetControlInputData(int enum_in,const char* data){/*{{{*/

	/*Get input id*/
	int id = EnumToIndex(enum_in);

	/*Check that it has the right format*/
	Input* input = this->inputs[id];
	if(!input) return NULL;
	if(input->ObjectEnum() != ControlInputEnum){
		_error_("Input "<<EnumToStringx(enum_in)<<" is not an ControlInput");
	}

	/*Cast and return*/
	return xDynamicCast<ControlInput*>(input)->GetInput(data);
}/*}}}*/
DatasetInput* Inputs::GetDatasetInput(int enum_in){/*{{{*/

	/*Get input id*/
	int id = EnumToIndex(enum_in);

	/*Check that it has the right format*/
	Input* input = this->inputs[id];
	if(!input) return NULL;
	if(input->ObjectEnum() != DatasetInputEnum){
		_error_("Input "<<EnumToStringx(enum_in)<<" is not an DatasetInput");
	}

	/*Cast and return*/
	return xDynamicCast<DatasetInput*>(input);
}/*}}}*/
ControlInput* Inputs::GetControlInput(int enum_in){/*{{{*/

	/*Get input id*/
	int id = EnumToIndex(enum_in);

	/*Check that it has the right format*/
	Input* input = this->inputs[id];
	if(!input) return NULL;
	if(input->ObjectEnum() != ControlInputEnum){
		_error_("Input "<<EnumToStringx(enum_in)<<" is not an ControlInput");
	}

	/*Cast and return*/
	return xDynamicCast<ControlInput*>(input);
}/*}}}*/
void Inputs::GetArrayPtr(int enum_in,int row,IssmDouble** pvalues,int* pN){/*{{{*/

	/*Get input id*/
	int id = EnumToIndex(enum_in);

	/*Create it if necessary*/
	if(this->inputs[id]){
		if(this->inputs[id]->ObjectEnum()!=ArrayInputEnum) _error_(EnumToStringx(this->inputs[id]->ObjectEnum())<<" cannot return an array");
	}
	else{
		_error_("Input "<<EnumToStringx(enum_in)<<" not found");
	}

	/*Set input*/
	ArrayInput* input = xDynamicCast<ArrayInput*>(this->inputs[id]);
	input->GetArrayPtr(row,pvalues,pN);
}/*}}}*/
void Inputs::GetIntArrayPtr(int enum_in,int row,int** pvalues,int* pN){/*{{{*/

	/*Get input id*/
	int id = EnumToIndex(enum_in);

	/*Create it if necessary*/
	if(this->inputs[id]){
		if(this->inputs[id]->ObjectEnum()!=IntArrayInputEnum) _error_(EnumToStringx(this->inputs[id]->ObjectEnum())<<" cannot return an int array");
	}
	else{
		_error_("Input "<<EnumToStringx(enum_in)<<" not found");
	}

	/*Set input*/
	IntArrayInput* input = xDynamicCast<IntArrayInput*>(this->inputs[id]);
	input->GetArrayPtr(row,pvalues,pN);
}/*}}}*/
void Inputs::GetArray(int enum_in,int row,IssmDouble** pvalues,int* pN){/*{{{*/

	/*Get input id*/
	int id = EnumToIndex(enum_in);

	/*Create it if necessary*/
	if(this->inputs[id]){
		if(this->inputs[id]->ObjectEnum()!=ArrayInputEnum) _error_(EnumToStringx(this->inputs[id]->ObjectEnum())<<" cannot return an array");
	}
	else{
		_error_("Input "<<EnumToStringx(enum_in)<<" not found");
	}

	/*Set input*/
	ArrayInput* input = xDynamicCast<ArrayInput*>(this->inputs[id]);
	input->GetArray(row,pvalues,pN);
}/*}}}*/
void Inputs::GetIntArray(int enum_in,int row,int** pvalues,int* pN){/*{{{*/

	/*Get input id*/
	int id = EnumToIndex(enum_in);

	/*Create it if necessary*/
	if(this->inputs[id]){
		if(this->inputs[id]->ObjectEnum()!=IntArrayInputEnum) _error_(EnumToStringx(this->inputs[id]->ObjectEnum())<<" cannot return an int array");
	}
	else{
		_error_("Input "<<EnumToStringx(enum_in)<<" not found");
	}

	/*Set input*/
	IntArrayInput* input = xDynamicCast<IntArrayInput*>(this->inputs[id]);
	input->GetArray(row,pvalues,pN);
}/*}}}*/
void Inputs::GetInputValue(bool* pvalue,int enum_in,int index){/*{{{*/

	/*Get input id*/
	int id = EnumToIndex(enum_in);

	/*Create it if necessary*/
	if(this->inputs[id]){
		if(this->inputs[id]->ObjectEnum()!=BoolInputEnum) _error_(EnumToStringx(this->inputs[id]->ObjectEnum())<<" cannot return a bool");
	}
	else{
		_error_("Input "<<EnumToStringx(enum_in)<<" not found");
	}

	/*Set input*/
	BoolInput* input = xDynamicCast<BoolInput*>(this->inputs[id]);
	input->GetInput(pvalue,index);
}/*}}}*/
void Inputs::GetInputValue(int* pvalue,int enum_in,int index){/*{{{*/

	/*Get input id*/
	int id = EnumToIndex(enum_in);

	/*Create it if necessary*/
	if(this->inputs[id]){
		if(this->inputs[id]->ObjectEnum()!=IntInputEnum) _error_(EnumToStringx(this->inputs[id]->ObjectEnum())<<" cannot return a int");
	}
	else{
		int* temp = xNew<int>(3);
		_error_("Input "<<EnumToStringx(enum_in)<<" not found");
	}

	/*Set input*/
	IntInput* input = xDynamicCast<IntInput*>(this->inputs[id]);
	input->GetInput(pvalue,index);
}/*}}}*/
void Inputs::GetInputValue(IssmDouble* pvalue,int enum_in,int index){/*{{{*/

	/*Get input id*/
	int id = EnumToIndex(enum_in);

	/*Create it if necessary*/
	if(this->inputs[id]){
		if(this->inputs[id]->ObjectEnum()!=DoubleInputEnum) _error_(EnumToStringx(this->inputs[id]->ObjectEnum())<<" cannot return a double!");
	}
	else{
		int* temp = xNew<int>(3);
		_error_("Input "<<EnumToStringx(enum_in)<<" not found");
	}

	/*Set input*/
	DoubleInput* input = xDynamicCast<DoubleInput*>(this->inputs[id]);
	input->GetInput(pvalue,index);

}/*}}}*/
bool Inputs::IsFileInputUpdate(IssmDouble time){/*{{{*/

	for(int i=0;i<NUMINPUTS;i++){
		if(this->inputs[i]){
			if(this->inputs[i]->ObjectEnum()==TransientFileInputEnum){
				TransientFileInput* input = xDynamicCast<TransientFileInput*>(this->inputs[i]);
				if(input->IsFileInputUpdate(time)) return true;
			}
		}
	}

	return false;
}/*}}}*/
void Inputs::ResultInterpolation(int* pinterpolation,int* pnodesperelement,int* parray_size, int output_enum){/*{{{*/

	/*Get input */
	int     index = EnumToIndex(output_enum);
	Input* input = this->inputs[index];

	/*Check that it is found*/
	if(!input){
		_error_("Input "<<EnumToStringx(output_enum)<<" not found and cannot be added to model results");
	}

	/*Assign output pointer*/
	*pinterpolation   = input->GetResultInterpolation();
	*pnodesperelement = input->GetResultNumberOfNodes();
	*parray_size      = input->GetResultArraySize();
}/*}}}*/
void Inputs::SetInput(int enum_in,int index,bool value){/*{{{*/

	/*Get input id*/
	int id = EnumToIndex(enum_in);

	/*Create it if necessary*/
	if(this->inputs[id]){
		if(this->inputs[id]->ObjectEnum()!=BoolInputEnum) _error_("cannot add a bool to a "<<EnumToStringx(this->inputs[id]->ObjectEnum()));
	}
	else{
		this->inputs[id] = new BoolInput(this->numberofelements_local);
	}

	/*Set input*/
	BoolInput* input = xDynamicCast<BoolInput*>(this->inputs[id]);
	input->SetInput(index,value);
}/*}}}*/
void Inputs::SetInput(int enum_in,int index,int value){/*{{{*/

	/*Get input id*/
	int id = EnumToIndex(enum_in);

	/*Create it if necessary*/
	if(this->inputs[id]){
		if(this->inputs[id]->ObjectEnum()!=IntInputEnum) _error_("cannot add an int to a "<<EnumToStringx(this->inputs[id]->ObjectEnum()));
	}
	else{
		this->inputs[id] = new IntInput(this->numberofelements_local);
	}

	/*Set input*/
	IntInput* input = xDynamicCast<IntInput*>(this->inputs[id]);
	input->SetInput(index,value);
}/*}}}*/
void Inputs::SetDoubleInput(int enum_in,int index,IssmDouble value){/*{{{*/

	/*Get input id*/
	int id = EnumToIndex(enum_in);

	/*Create it if necessary*/
	if(this->inputs[id]){
		if(this->inputs[id]->ObjectEnum()!=DoubleInputEnum) _error_("cannot add a double to a "<<EnumToStringx(this->inputs[id]->ObjectEnum()));
	}
	else{
		this->inputs[id] = new DoubleInput(this->numberofelements_local);
	}

	/*Set input*/
	DoubleInput* input = xDynamicCast<DoubleInput*>(this->inputs[id]);
	input->SetInput(index,value);
}/*}}}*/
void Inputs::SetArrayInput(int enum_in,int row,IssmDouble* values,int numlayers){/*{{{*/

	bool recreate = false;

	/*Get input id*/
	int id = EnumToIndex(enum_in);

	/*Create it if necessary*/
	if(this->inputs[id]){
		if(this->inputs[id]->ObjectEnum()!=ArrayInputEnum){
			delete this->inputs[id];
			recreate = true;
		}
	}
	else{
		recreate = true;
	}

	if(recreate){
		this->inputs[id] = new ArrayInput(this->numberofelements_local);
	}

	/*Set input*/
	ArrayInput* input = xDynamicCast<ArrayInput*>(this->inputs[id]);
	input->SetInput(row,numlayers,values);
}/*}}}*/
void Inputs::SetIntArrayInput(int enum_in,int row,int* values,int numlayers){/*{{{*/

	bool recreate = false;

	/*Get input id*/
	int id = EnumToIndex(enum_in);

	/*Create it if necessary*/
	if(this->inputs[id]){
		if(this->inputs[id]->ObjectEnum()!=IntArrayInputEnum){
			delete this->inputs[id];
			recreate = true;
		}
	}
	else{
		recreate = true;
	}

	if(recreate){
		this->inputs[id] = new IntArrayInput(this->numberofelements_local);
	}

	/*Set input*/
	IntArrayInput* input = xDynamicCast<IntArrayInput*>(this->inputs[id]);
	input->SetInput(row,numlayers,values);
}/*}}}*/
TransientInput* Inputs::SetDatasetTransientInput(int enum_in,int dataset_id,IssmDouble* times,int numtimes){/*{{{*/

	bool recreate = false;
	/*Get input id*/
	int id = EnumToIndex(enum_in);

	/*Create it if necessary*/
	if(this->inputs[id]){
		if(this->inputs[id]->ObjectEnum()!=DatasetInputEnum){
			delete this->inputs[id];
			recreate = true;
		}
	}
	else{
		recreate = true;
	}

	if(recreate){
		this->inputs[id] = new DatasetInput(this->numberofelements_local,this->numberofvertices_local);
	}

	/*Get Dataset Input now*/
	DatasetInput* input = xDynamicCast<DatasetInput*>(this->inputs[id]);

	/*Create and return transient input*/
	return input->SetTransientInput(dataset_id,times,numtimes);
}/*}}}*/
void Inputs::SetTransientInput(int enum_in,IssmDouble* times,int numtimes){/*{{{*/

	/*Get input id*/
	int id = EnumToIndex(enum_in);

	/*Create it if necessary*/
	if(this->inputs[id]){
		/*Input already there, make sure it is the right type*/
		if(this->inputs[id]->ObjectEnum()!=TransientInputEnum){
			_error_("cannot add a TransientInput to a "<<EnumToStringx(this->inputs[id]->ObjectEnum()));
		}
	}
	else{
		this->inputs[id] = new TransientInput(enum_in,this->numberofelements_local,this->numberofvertices_local,times,numtimes);
	}
}/*}}}*/
void Inputs::SetControlInput(int enum_in,int layout,int interpolation,int control_id){/*{{{*/

	bool recreate = false;

	/*Get input id*/
	int id = EnumToIndex(enum_in);

	/*Create it if necessary*/
	if(this->inputs[id]){
		if(this->inputs[id]->ObjectEnum()!=ControlInputEnum){
			delete this->inputs[id];
			recreate = true;
		}
	}
	else{
		recreate = true;
	}

	if(recreate){
		this->inputs[id] = new ControlInput(this->numberofelements_local,this->numberofvertices_local,layout,interpolation,control_id);
	}

}/*}}}*/
void Inputs::SetTransientControlInput(int enum_in,int control_id,IssmDouble* times,int numtimes){/*{{{*/

	bool recreate = false;

	/*Get input id*/
	int id = EnumToIndex(enum_in);

	/*Create it if necessary*/
	if(this->inputs[id]){
		if(this->inputs[id]->ObjectEnum()!=ControlInputEnum){
			delete this->inputs[id];
			recreate = true;
		}
	}
	else{
		recreate = true;
	}

	if(recreate){
		this->inputs[id] = new ControlInput(enum_in,this->numberofelements_local,this->numberofvertices_local,control_id,times,numtimes);
	}

}/*}}}*/
void Inputs::SetTriaControlInputGradient(int enum_in,int interpolation,int numindices,int* indices,IssmDouble* values){/*{{{*/

	/*Get input id*/
	int id = EnumToIndex(enum_in);

	/*Create it if necessary*/
	if(!this->inputs[id]) _error_("could not find Input "<<EnumToStringx(enum_in));
	if( this->inputs[id]->ObjectEnum()!=ControlInputEnum) _error_("Input "<<EnumToStringx(enum_in)<<" is not a ControlInput");

	/*Set input*/
	ControlInput* input = xDynamicCast<ControlInput*>(this->inputs[id]);
	input->SetGradient(interpolation,numindices,indices,values);
}/*}}}*/
void Inputs::SetTriaControlInputGradient(int enum_in,int interpolation,int numindices,int* indices,IssmDouble* values,int n){/*{{{*/

	/*Get input id*/
	int id = EnumToIndex(enum_in);

	/*Create it if necessary*/
	if(!this->inputs[id]) _error_("could not find Input "<<EnumToStringx(enum_in));
	if( this->inputs[id]->ObjectEnum()!=ControlInputEnum) _error_("Input "<<EnumToStringx(enum_in)<<" is not a ControlInput");

	/*Set input*/
	ControlInput* input = xDynamicCast<ControlInput*>(this->inputs[id]);
	input->SetGradient(interpolation,numindices,indices,values,n);
}/*}}}*/
void Inputs::SetTriaDatasetInput(int enum_in,int id_in,int interpolation,int numindices,int* indices,IssmDouble* values){/*{{{*/

	bool recreate = false;
	/*Get input id*/
	int id = EnumToIndex(enum_in);

	/*Create it if necessary*/
	if(this->inputs[id]){
		if(this->inputs[id]->ObjectEnum()!=DatasetInputEnum){
			delete this->inputs[id];
			recreate = true;
		}
	}
	else{
		recreate = true;
	}

	if(recreate){
		this->inputs[id] = new DatasetInput(this->numberofelements_local,this->numberofvertices_local);
	}

	/*Set input*/
	DatasetInput* input = xDynamicCast<DatasetInput*>(this->inputs[id]);
	input->SetTriaInput(id_in,P1Enum,numindices,indices,values);
}/*}}}*/
void Inputs::SetTriaInput(int enum_in,int interpolation,int row,IssmDouble value){/*{{{*/

	/*This one only supports P0 and P1 because it assumes col=0*/
	_assert_(interpolation==P0Enum || interpolation==P1Enum);

	/*Get input id*/
	int id = EnumToIndex(enum_in);

	/*Create it if necessary*/
	if(this->inputs[id]){
		if(this->inputs[id]->ObjectEnum()!=TriaInputEnum) _error_("cannot add a bool to a "<<EnumToStringx(this->inputs[id]->ObjectEnum()));
	}
	else{
		this->inputs[id] = new TriaInput(this->numberofelements_local,this->numberofvertices_local,interpolation);
	}

	/*Set input*/
	TriaInput* input = xDynamicCast<TriaInput*>(this->inputs[id]);
	input->SetInput(interpolation,row,value);
}/*}}}*/
void Inputs::SetTriaInput(int enum_in,int interpolation,int numindices,int* indices,IssmDouble* values){/*{{{*/

	/*Get input id*/
	int id = EnumToIndex(enum_in);

	/*Create it if necessary*/
	if(this->inputs[id]){
		if(this->inputs[id]->ObjectEnum()!=TriaInputEnum){
			_error_("cannot add Element values to a "<<EnumToStringx(this->inputs[id]->ObjectEnum())<<" while trying to set "<<EnumToStringx(enum_in));
		}
	}
	else{
		this->inputs[id] = new TriaInput(this->numberofelements_local,this->numberofvertices_local,interpolation);
	}

	/*Set input*/
	TriaInput* input = xDynamicCast<TriaInput*>(this->inputs[id]);
	input->SetInput(interpolation,numindices,indices,values);
}/*}}}*/
void Inputs::SetTriaInput(int enum_in,int interpolation,int row,int numindices,IssmDouble* values){/*{{{*/

	/*Get input id*/
	int id = EnumToIndex(enum_in);

	/*Create it if necessary*/
	if(this->inputs[id]){
		if(this->inputs[id]->ObjectEnum()!=TriaInputEnum) _error_("cannot add Element values to a "<<EnumToStringx(this->inputs[id]->ObjectEnum()));
	}
	else{
		this->inputs[id] = new TriaInput(this->numberofelements_local,this->numberofvertices_local,interpolation);
	}

	/*Set input*/
	TriaInput* input = xDynamicCast<TriaInput*>(this->inputs[id]);
	input->SetInput(interpolation,row,numindices,values);
}/*}}}*/
void Inputs::SetPentaControlInputGradient(int enum_in,int interpolation,int numindices,int* indices,IssmDouble* values){/*{{{*/

	/*Get input id*/
	int id = EnumToIndex(enum_in);

	/*Create it if necessary*/
	if(!this->inputs[id]) _error_("could not find Input "<<EnumToStringx(enum_in));
	if( this->inputs[id]->ObjectEnum()!=ControlInputEnum) _error_("Input "<<EnumToStringx(enum_in)<<" is not a ControlInput");

	/*Set input*/
	ControlInput* input = xDynamicCast<ControlInput*>(this->inputs[id]);
	input->SetGradient(interpolation,numindices,indices,values);
}/*}}}*/
void Inputs::SetPentaDatasetInput(int enum_in,int id_in,int interpolation,int numindices,int* indices,IssmDouble* values){/*{{{*/

	bool recreate = false;
	/*Get input id*/
	int id = EnumToIndex(enum_in);

	/*Create it if necessary*/
	if(this->inputs[id]){
		if(this->inputs[id]->ObjectEnum()!=DatasetInputEnum){
			delete this->inputs[id];
			recreate = true;
		}
	}
	else{
		recreate = true;
	}

	if(recreate){
		this->inputs[id] = new DatasetInput(this->numberofelements_local,this->numberofvertices_local);
	}

	/*Set input*/
	DatasetInput* input = xDynamicCast<DatasetInput*>(this->inputs[id]);
	input->SetPentaInput(id_in,P1Enum,numindices,indices,values);
}/*}}}*/
void Inputs::SetPentaInput(int enum_in,int interpolation,int row,IssmDouble value){/*{{{*/

	/*This one only supports P0 and P1 because it assumes col=0*/
	_assert_(interpolation==P0Enum || interpolation==P1Enum);

	/*Get input id*/
	int id = EnumToIndex(enum_in);

	/*Create it if necessary*/
	if(this->inputs[id]){
		if(this->inputs[id]->ObjectEnum()!=PentaInputEnum) _error_("cannot add a bool to a "<<EnumToStringx(this->inputs[id]->ObjectEnum()));
	}
	else{
		this->inputs[id] = new PentaInput(this->numberofelements_local,this->numberofvertices_local,interpolation);
	}

	/*Set input*/
	PentaInput* input = xDynamicCast<PentaInput*>(this->inputs[id]);
	input->SetInput(interpolation,row,value);
}/*}}}*/
void Inputs::SetPentaInput(int enum_in,int interpolation,int numindices,int* indices,IssmDouble* values){/*{{{*/

	/*Get input id*/
	int id = EnumToIndex(enum_in);

	/*Create it if necessary*/
	if(this->inputs[id]){
		if(this->inputs[id]->ObjectEnum()!=PentaInputEnum) _error_("cannot add Element values to a "<<EnumToStringx(this->inputs[id]->ObjectEnum()));
	}
	else{
		this->inputs[id] = new PentaInput(this->numberofelements_local,this->numberofvertices_local,interpolation);
	}

	/*Set input*/
	PentaInput* input = xDynamicCast<PentaInput*>(this->inputs[id]);
	input->SetInput(interpolation,numindices,indices,values);
}/*}}}*/
void Inputs::SetPentaInput(int enum_in,int interpolation,int row,int numindices,IssmDouble* values){/*{{{*/

	/*Get input id*/
	int id = EnumToIndex(enum_in);

	/*Create it if necessary*/
	if(this->inputs[id]){
		if(this->inputs[id]->ObjectEnum()!=PentaInputEnum) _error_("cannot add Element values to a "<<EnumToStringx(this->inputs[id]->ObjectEnum()));
	}
	else{
		this->inputs[id] = new PentaInput(this->numberofelements_local,this->numberofvertices_local,interpolation);
	}

	/*Set input*/
	PentaInput* input = xDynamicCast<PentaInput*>(this->inputs[id]);
	input->SetInput(interpolation,row,numindices,values);
}/*}}}*/

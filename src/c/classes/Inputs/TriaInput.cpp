/*!\file TriaInput.c
 * \brief: implementation of the TriaInput object
 */

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "../classes.h"
#include "../../shared/shared.h"
#include "./TriaInput.h"

/*TriaInput constructors and destructor*/
TriaInput::TriaInput(void){/*{{{*/

	this->numberofelements_local = -1;
	this->numberofvertices_local = -1;
	this->isserved       = false;
	this->isserved_collapsed= 0;
	this->M = -1;
	this->N = -1;
	this->values         = NULL;
	this->element_values = NULL;

}/*}}}*/
TriaInput::TriaInput(int nbe_in,int nbv_in,int interp_in){/*{{{*/

	_assert_(nbe_in>0);
	_assert_(nbe_in<1e11);
	_assert_(nbv_in>0);
	_assert_(nbv_in<1e11);
	this->numberofelements_local = nbe_in;
	this->numberofvertices_local = nbv_in;
	this->isserved       = false;
	this->isserved_collapsed = 0;

	/*Reset takes care of the rest*/
	this->Reset(interp_in);
}/*}}}*/
TriaInput::~TriaInput(){/*{{{*/
	if(this->element_values) xDelete<IssmDouble>(this->element_values);
	if(this->values)         xDelete<IssmDouble>(this->values);
}
/*}}}*/
void TriaInput::Reset(int interp_in){/*{{{*/

	/*Clean up*/
	if(this->values)         xDelete<IssmDouble>(this->values);
	if(this->element_values) xDelete<IssmDouble>(this->element_values);

	/*Set interpolation*/
	this->interpolation  = interp_in;

	/*Create Sizes*/
	if(this->interpolation==P1Enum){
		this->M = this->numberofvertices_local;
		this->N = 1;
	}
	else{
		this->M = this->numberofelements_local;
		this->N = TriaRef::NumberofNodes(interp_in);
	}

	/*Allocate Pointers*/
	this->values         = xNewZeroInit<IssmDouble>(this->M*this->N);
	this->element_values = xNewZeroInit<IssmDouble>(TriaRef::NumberofNodes(interp_in));
}/*}}}*/

/*Object virtual functions definitions:*/
Input* TriaInput::copy() {/*{{{*/

	TriaInput* output = new TriaInput(this->numberofelements_local,this->numberofvertices_local,this->interpolation);

	xMemCpy<IssmDouble>(output->values,this->values,this->M*this->N);
	xMemCpy<IssmDouble>(output->element_values,this->element_values,TriaRef::NumberofNodes(this->interpolation));

	return output;
}
/*}}}*/
void TriaInput::DeepEcho(void){/*{{{*/
	_printf_("TriaInput Echo:\n");
	_printf_("   interpolation: "<<EnumToStringx(this->interpolation)<<"\n");
	_printf_("   Size:          "<<M<<"x"<<N<<"\n");
	_printf_("   isserved:      "<<(isserved?"true":"false") << "\n");
	_printf_("   isserved_collapsed: "<<isserved_collapsed << "\n");
	if(isserved){
		_printf_("   current values:      ");
		for(int i=0;i<3;i++) _printf_(" "<<this->element_values[i]);
		_printf_("] ("<<EnumToStringx(this->interpolation)<<")\n");
	}
	printarray(this->values,this->M,this->N);
	//_printf_(setw(15)<<"   TriaInput "<<setw(25)<<left<<EnumToStringx(this->enum_type)<<" "<<(value?"true":"false") << "\n");
}
/*}}}*/
void TriaInput::Echo(void){/*{{{*/
	_printf_("TriaInput Echo:\n");
	_printf_("   interpolation: "<<EnumToStringx(this->interpolation)<<"\n");
	_printf_("   Size:          "<<M<<"x"<<N<<"\n");
	_printf_("   isserved:      "<<(isserved?"true":"false") << "\n");
	_printf_("   isserved_collapsed: "<<isserved_collapsed << "\n");
	if(isserved){
		_printf_("   current values:      ");
		_printf_("[ ");
		for(int i=0;i<TriaRef::NumberofNodes(this->interpolation);i++) _printf_(" "<<this->element_values[i]);
		_printf_("] ("<<EnumToStringx(this->interpolation)<<")\n");
	}
}
/*}}}*/
int  TriaInput::Id(void){/*{{{*/
	return -1;
}/*}}}*/
void TriaInput::Marshall(MarshallHandle* marshallhandle){ /*{{{*/

	int object_enum = TriaInputEnum;
   marshallhandle->call(object_enum);

	marshallhandle->call(this->numberofelements_local);
	marshallhandle->call(this->numberofvertices_local);
	marshallhandle->call(this->interpolation);
	marshallhandle->call(this->M);
	marshallhandle->call(this->N);
	this->isserved = false;
	this->isserved_collapsed = 0;
	if(this->M*this->N){
		marshallhandle->call(this->values,this->M*this->N);
	}
	else this->values = NULL;

	if(marshallhandle->OperationNumber() == MARSHALLING_LOAD){
		this->element_values = xNewZeroInit<IssmDouble>(TriaRef::NumberofNodes(this->interpolation));
	}

}
/*}}}*/
int  TriaInput::ObjectEnum(void){/*{{{*/
	return TriaInputEnum;
}
/*}}}*/

/*TriaInput management*/
void TriaInput::SetInput(int interp_in,int row,IssmDouble value_in){/*{{{*/

	_assert_(this);
	_assert_(row>=0);
	_assert_(row<this->M);
	_assert_(this->N==1);

	this->values[row] = value_in;
	this->isserved = false;
}
/*}}}*/
void TriaInput::SetInput(int interp_in,int numindices,int* indices,IssmDouble* values_in){/*{{{*/

	_assert_(this);
	if(interp_in==P1Enum && this->interpolation==P1Enum){
		_assert_(this->N==1);
		for(int i=0;i<numindices;i++){
			int row = indices[i];
			_assert_(row>=0);
			_assert_(row<this->M);
			this->values[row] = values_in[i];
		}
	}
	else if(interp_in==P0Enum && this->interpolation==P0Enum){
		_assert_(this->N==1);
		for(int i=0;i<numindices;i++){
			int row = indices[i];
			_assert_(row>=0);
			_assert_(row<this->M);
			this->values[row] = values_in[i];
		}
	}
	else if(this->interpolation!=P1Enum && interp_in==P1Enum){
		this->Reset(interp_in);
		for(int i=0;i<numindices;i++){
			int row = indices[i];
			_assert_(row>=0);
			_assert_(row<this->M);
			this->values[row] = values_in[i];
		}
	}
	else{
		_error_("Cannot convert "<<EnumToStringx(this->interpolation)<<" to "<<EnumToStringx(interp_in));
	}
	this->isserved = false;
}
/*}}}*/
void TriaInput::SetInput(int interp_in,int row,int numindices,IssmDouble* values_in){/*{{{*/

	_assert_(this);
	if(interp_in==this->interpolation){
		_assert_(this->N==numindices);
	}
	else{
		this->Reset(interp_in);
		_assert_(this->N==numindices);
	}
	for(int i=0;i<numindices;i++) this->values[row*this->N+i] = values_in[i];
	this->isserved = false;
}
/*}}}*/
void TriaInput::Serve(int numindices,int* indices){/*{{{*/

	_assert_(this);
	_assert_(this->N==1);

	for(int i=0;i<numindices;i++){
		int row = indices[i];
		_assert_(row>=0);
		_assert_(row<this->M);
		this->element_values[i] = this->values[row];
	}

	/*Set input as served*/
	this->isserved = true;
	this->isserved_collapsed = 0;
}
/*}}}*/
void TriaInput::Serve(int row,int numindices){/*{{{*/

	_assert_(this);
	_assert_(this->N==numindices);
	_assert_(row<this->M);
	_assert_(row>=0);

	for(int i=0;i<numindices;i++){
		this->element_values[i] = this->values[row*this->N+i];
	}

	/*Set input as served*/
	this->isserved = true;
	this->isserved_collapsed = 0;
} /*}}}*/
void TriaInput::ServeCollapsed(int row,int id1,int id2){/*{{{*/

	_assert_(this);
	_assert_(this->N>=3);
	_assert_(row<this->M);
	_assert_(row>=0);
	_assert_(id1>=0 && id1<3);
	_assert_(id2>=0 && id2<3);

	this->element_values[0] = this->values[row*this->N+id1];
	this->element_values[1] = this->values[row*this->N+id2];

	/*Set input as served*/
	this->isserved = true;
	this->isserved_collapsed = 1;
}/*}}}*/
void TriaInput::SetServeCollapsed(bool status){/*{{{*/
	this->isserved_collapsed = 1;
}/*}}}*/
int  TriaInput::GetInterpolation(){/*{{{*/
	return this->interpolation;
}/*}}}*/
void TriaInput::GetInputAverage(IssmDouble* pvalue){/*{{{*/
	_assert_(this);
	_assert_(this->isserved);

	int        numnodes  = this->NumberofNodes(this->interpolation);
	if(this->isserved_collapsed) numnodes = 2;
	IssmDouble numnodesd = reCast<int,IssmDouble>(numnodes);
	IssmDouble value     = 0.;

	for(int i=0;i<numnodes;i++) value+=this->element_values[i];
	value = value/numnodesd;

	*pvalue=value;
}/*}}}*/
IssmDouble TriaInput::GetInputMin(void){/*{{{*/
	_assert_(this);
	_assert_(this->isserved);

	int        numnodes  = this->NumberofNodes(this->interpolation);
	if(this->isserved_collapsed) numnodes = 2;
	IssmDouble min=this->element_values[0];

	for(int i=1;i<numnodes;i++){
		if(this->element_values[i]<min) min=this->element_values[i];
	}
	return min;
}/*}}}*/
IssmDouble TriaInput::GetInputMax(void){/*{{{*/
	_assert_(this);
	_assert_(this->isserved);

	int        numnodes  = this->NumberofNodes(this->interpolation);
	if(this->isserved_collapsed) numnodes = 2;
	IssmDouble max=this->element_values[0];

	for(int i=1;i<numnodes;i++){
		if(this->element_values[i]>max) max=this->element_values[i];
	}
	return max;
}/*}}}*/
IssmDouble TriaInput::GetInputMaxAbs(void){/*{{{*/
	_assert_(this);
	_assert_(this->isserved);

	int        numnodes  = this->NumberofNodes(this->interpolation);
	if(this->isserved_collapsed) numnodes = 2;
	IssmDouble maxabs=fabs(this->element_values[0]);

	for(int i=1;i<numnodes;i++){
		if(fabs(this->element_values[i])>maxabs) maxabs=fabs(this->element_values[i]);
	}
	return maxabs;
}/*}}}*/
void TriaInput::GetInputDerivativeValue(IssmDouble* derivativevalues, IssmDouble* xyz_list, Gauss* gauss){/*{{{*/
	_assert_(this);
	_assert_(this->isserved);

	if(this->isserved_collapsed){
		_assert_(gauss->Enum()==GaussSegEnum);
		if(this->interpolation==P0Enum){
			derivativevalues[0] = 0.;
		}
		else{
			SegRef temp;
			temp.GetInputDerivativeValue(derivativevalues,this->element_values,xyz_list,(GaussSeg*)gauss,P1Enum);
		}
	}
	else{
		_assert_(gauss->Enum()==GaussTriaEnum);
		TriaRef::GetInputDerivativeValue(derivativevalues,this->element_values,xyz_list,(GaussTria*)gauss,this->interpolation);
	}
}/*}}}*/
void TriaInput::GetInputValue(IssmDouble* pvalue,Gauss* gauss){/*{{{*/
	_assert_(this);
	_assert_(this->isserved);
	if(this->isserved_collapsed){
		_assert_(gauss->Enum()==GaussSegEnum);
		if(this->interpolation==P0Enum){
			*pvalue = this->element_values[0];
		}
		else{
			SegRef temp;
			temp.GetInputValue(pvalue,this->element_values,(GaussSeg*)gauss,P1Enum);
		}
	}
	else{
		_assert_(gauss->Enum()==GaussTriaEnum);
		TriaRef::GetInputValue(pvalue,this->element_values,(GaussTria*)gauss,this->interpolation);
	}
}/*}}}*/
int  TriaInput::GetResultArraySize(void){/*{{{*/
	return 1;
}
/*}}}*/
int  TriaInput::GetResultInterpolation(void){/*{{{*/
	if(this->interpolation==P0Enum || this->interpolation==P0DGEnum){
		return P0Enum;
	}
	return P1Enum;
}/*}}}*/
int  TriaInput::GetResultNumberOfNodes(void){/*{{{*/
	return TriaRef::NumberofNodes(this->interpolation);
}
/*}}}*/
void TriaInput::Scale(IssmDouble alpha){/*{{{*/

	for(int i=0;i<this->M*this->N;i++) this->values[i] = alpha*this->values[i];
	for(int i=0;i<TriaRef::NumberofNodes(this->interpolation);i++) this->element_values[i] = alpha*this->element_values[i];
}
/*}}}*/
void TriaInput::Pow(IssmDouble alpha){/*{{{*/

	for(int i=0;i<this->M*this->N;i++) this->values[i] = pow(this->values[i],alpha);
	for(int i=0;i<TriaRef::NumberofNodes(this->interpolation);i++) this->element_values[i] = pow(this->element_values[i],alpha);
}
/*}}}*/
void TriaInput::AXPY(Input* xinput,IssmDouble alpha){/*{{{*/

	/*xinput is of the same type, so cast it: */
	if(xinput->ObjectEnum()!=TriaInputEnum) _error_("Operation not permitted because xinput is of type " << EnumToStringx(xinput->ObjectEnum()));
	TriaInput* xtriainput=xDynamicCast<TriaInput*>(xinput);
	if(xtriainput->GetInterpolation()!=this->interpolation) _error_("Operation not permitted because xinput is of type " << EnumToStringx(xinput->ObjectEnum()));

	/*Carry out the AXPY operation depending on type:*/
	for(int i=0;i<this->M*this->N;i++) this->values[i] = alpha*xtriainput->values[i] + this->values[i];
	for(int i=0;i<TriaRef::NumberofNodes(this->interpolation);i++) this->element_values[i] = alpha*xtriainput->element_values[i] + this->element_values[i];
}
/*}}}*/
void TriaInput::Shift(IssmDouble alpha){/*{{{*/

	/*Carry out the shift operation:*/
	for(int i=0;i<this->M*this->N;i++) this->values[i] +=alpha;
	for(int i=0;i<TriaRef::NumberofNodes(this->interpolation);i++) this->element_values[i] += alpha;
}
/*}}}*/
void TriaInput::PointWiseMult(Input* xinput){/*{{{*/

	/*xinput is of the same type, so cast it: */
	if(xinput->ObjectEnum()!=TriaInputEnum) _error_("Operation not permitted because xinput is of type " << EnumToStringx(xinput->ObjectEnum()));
	TriaInput* xtriainput=xDynamicCast<TriaInput*>(xinput);
	if(xtriainput->GetInterpolation()!=this->interpolation) _error_("Operation not permitted because xinput is of type " << EnumToStringx(xinput->ObjectEnum()));

	/* we need to check that the vector sizes are identical*/
	if(xtriainput->M!=this->M||xtriainput->N!=this->N) _error_("Operation not permitted because the inputs have different sizes");

	/*Carry out the AXPY operation depending on type:*/
	for(int i=0;i<this->M*this->N;i++) this->values[i] = xtriainput->values[i] * this->values[i];
	for(int i=0;i<TriaRef::NumberofNodes(this->interpolation);i++) this->element_values[i] = xtriainput->element_values[i] * this->element_values[i];
}
/*}}}*/
void TriaInput::AverageAndReplace(void){/*{{{*/

	if(this->M!=this->numberofelements_local) _error_("not implemented for P1");

	/*Get local sum and local size*/
	IssmDouble sum  = 0.;
	int        weight;
	for(int i=0;i<this->M*this->N;i++) sum += this->values[i];
	weight = this->M*this->N;

	/*Get sum across all procs*/
	IssmDouble all_sum;
	int        all_weight;
	ISSM_MPI_Allreduce((void*)&sum,(void*)&all_sum,1,ISSM_MPI_DOUBLE,ISSM_MPI_SUM,IssmComm::GetComm());
	ISSM_MPI_Allreduce((void*)&weight,(void*)&all_weight,1,ISSM_MPI_INT,ISSM_MPI_SUM,IssmComm::GetComm());

	/*Divide by number of procs*/
	IssmDouble newvalue = all_sum/reCast<IssmPDouble>(all_weight);

	/*Now replace existing input*/
	this->Reset(P0Enum);
	for(int i=0;i<this->M*this->N;i++) this->values[i] = newvalue;
}
/*}}}*/

/*Object functions*/

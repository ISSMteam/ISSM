/*!\file SegInput.c
 * \brief: implementation of the SegInput object
 */

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "../classes.h"
#include "../../shared/shared.h"
#include "./SegInput.h"

/*SegInput constructors and destructor*/
SegInput::SegInput(void){/*{{{*/

	this->numberofelements_local = -1;
	this->numberofvertices_local = -1;
	this->isserved       = false;
	this->M = -1;
	this->N = -1;
	this->values         = NULL;
	this->element_values = NULL;

}/*}}}*/
SegInput::SegInput(int nbe_in,int nbv_in,int interp_in){/*{{{*/

	_assert_(nbe_in>0);
	_assert_(nbe_in<1e11);
	_assert_(nbv_in>0);
	_assert_(nbv_in<1e11);
	this->numberofelements_local = nbe_in;
	this->numberofvertices_local = nbv_in;
	this->isserved       = false;

	/*Reset takes care of the rest*/
	this->Reset(interp_in);
}/*}}}*/
SegInput::~SegInput(){/*{{{*/
	if(this->element_values) xDelete<IssmDouble>(this->element_values);
	if(this->values)         xDelete<IssmDouble>(this->values);
}
/*}}}*/
void SegInput::Reset(int interp_in){/*{{{*/

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
		this->N = SegRef::NumberofNodes(interp_in);
	}

	/*Allocate Pointers*/
	this->values         = xNewZeroInit<IssmDouble>(this->M*this->N);
	this->element_values = xNewZeroInit<IssmDouble>(SegRef::NumberofNodes(interp_in));
}/*}}}*/

/*Object virtual functions definitions:*/
Input* SegInput::copy() {/*{{{*/

	SegInput* output = new SegInput(this->numberofelements_local,this->numberofvertices_local,this->interpolation);

	xMemCpy<IssmDouble>(output->values,this->values,this->M*this->N);
	xMemCpy<IssmDouble>(output->element_values,this->element_values,SegRef::NumberofNodes(this->interpolation));

	return output;
}
/*}}}*/
void SegInput::DeepEcho(void){/*{{{*/
	_printf_("SegInput Echo:\n");
	_printf_("   interpolation: "<<EnumToStringx(this->interpolation)<<"\n");
	_printf_("   Size:          "<<M<<"x"<<N<<"\n");
	_printf_("   isserved:      "<<(isserved?"true":"false") << "\n");
	if(isserved){
		_printf_("   current values:      ");
		for(int i=0;i<3;i++) _printf_(" "<<this->element_values[i]);
		_printf_("] ("<<EnumToStringx(this->interpolation)<<")\n");
	}
	printarray(this->values,this->M,this->N);
	//_printf_(setw(15)<<"   SegInput "<<setw(25)<<left<<EnumToStringx(this->enum_type)<<" "<<(value?"true":"false") << "\n");
}
/*}}}*/
void SegInput::Echo(void){/*{{{*/
	_printf_("SegInput Echo:\n");
	_printf_("   interpolation: "<<EnumToStringx(this->interpolation)<<"\n");
	_printf_("   Size:          "<<M<<"x"<<N<<"\n");
	_printf_("   isserved:      "<<(isserved?"true":"false") << "\n");
	if(isserved){
		_printf_("   current values:      ");
		_printf_("[ ");
		for(int i=0;i<3;i++) _printf_(" "<<this->element_values[i]);
		_printf_("] ("<<EnumToStringx(this->interpolation)<<")\n");
	}
}
/*}}}*/
int  SegInput::Id(void){/*{{{*/
	return -1;
}/*}}}*/
void SegInput::Marshall(MarshallHandle* marshallhandle){ /*{{{*/

	int object_enum = SegInputEnum;
   marshallhandle->call(object_enum);

	marshallhandle->call(this->numberofelements_local);
	marshallhandle->call(this->numberofvertices_local);
	marshallhandle->call(this->interpolation);
	marshallhandle->call(this->M);
	marshallhandle->call(this->N);
	this->isserved = false;
	if(this->M*this->N){
		marshallhandle->call(this->values,this->M*this->N);
	}
	else this->values = NULL;

	if(marshallhandle->OperationNumber() == MARSHALLING_LOAD){
		this->element_values = xNewZeroInit<IssmDouble>(SegRef::NumberofNodes(this->interpolation));
	}

}
/*}}}*/
int  SegInput::ObjectEnum(void){/*{{{*/
	return SegInputEnum;
}
/*}}}*/

/*SegInput management*/
void SegInput::SetInput(int interp_in,int row,IssmDouble value_in){/*{{{*/

	_assert_(this);
	_assert_(row>=0);
	_assert_(row<this->M);
	_assert_(this->N==1);

	this->values[row] = value_in;
	this->isserved = false;
}
/*}}}*/
void SegInput::SetInput(int interp_in,int numindices,int* indices,IssmDouble* values_in){/*{{{*/

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
void SegInput::SetInput(int interp_in,int row,int numindices,IssmDouble* values_in){/*{{{*/

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
void SegInput::Serve(int numindices,int* indices){/*{{{*/

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
}
/*}}}*/
void SegInput::Serve(int row,int numindices){/*{{{*/

	_assert_(this);
	_assert_(this->N==numindices);
	_assert_(row<this->M);
	_assert_(row>=0);

	for(int i=0;i<numindices;i++){
		this->element_values[i] = this->values[row*this->N+i];
	}

	/*Set input as served*/
	this->isserved = true;
} /*}}}*/
int  SegInput::GetInterpolation(){/*{{{*/
	return this->interpolation;
}/*}}}*/
void SegInput::GetInputAverage(IssmDouble* pvalue){/*{{{*/
	_assert_(this);
	_assert_(this->isserved);

	int        numnodes  = this->NumberofNodes(this->interpolation);
	IssmDouble numnodesd = reCast<int,IssmDouble>(numnodes);
	IssmDouble value     = 0.;

	for(int i=0;i<numnodes;i++) value+=this->element_values[i];
	value = value/numnodesd;

	*pvalue=value;
}/*}}}*/
IssmDouble SegInput::GetInputMin(void){/*{{{*/
	_assert_(this);
	_assert_(this->isserved);

	int        numnodes  = this->NumberofNodes(this->interpolation);
	IssmDouble min=this->element_values[0];

	for(int i=1;i<numnodes;i++){
		if(this->element_values[i]<min) min=this->element_values[i];
	}
	return min;
}/*}}}*/
IssmDouble SegInput::GetInputMax(void){/*{{{*/
	_assert_(this);
	_assert_(this->isserved);

	int        numnodes  = this->NumberofNodes(this->interpolation);
	IssmDouble max=this->element_values[0];

	for(int i=1;i<numnodes;i++){
		if(this->element_values[i]>max) max=this->element_values[i];
	}
	return max;
}/*}}}*/
IssmDouble SegInput::GetInputMaxAbs(void){/*{{{*/
	_assert_(this);
	_assert_(this->isserved);

	int        numnodes  = this->NumberofNodes(this->interpolation);
	IssmDouble maxabs=fabs(this->element_values[0]);

	for(int i=1;i<numnodes;i++){
		if(fabs(this->element_values[i])>maxabs) maxabs=fabs(this->element_values[i]);
	}
	return maxabs;
}/*}}}*/
void SegInput::GetInputDerivativeValue(IssmDouble* derivativevalues, IssmDouble* xyz_list, Gauss* gauss){/*{{{*/
	_assert_(this);
	_assert_(this->isserved);
	_assert_(gauss->Enum()==GaussSegEnum);
	SegRef::GetInputDerivativeValue(derivativevalues,this->element_values,xyz_list,(GaussSeg*)gauss,this->interpolation);
}/*}}}*/
void SegInput::GetInputValue(IssmDouble* pvalue,Gauss* gauss){/*{{{*/
	_assert_(this);
	_assert_(this->isserved);
	_assert_(gauss->Enum()==GaussSegEnum);
	SegRef::GetInputValue(pvalue,this->element_values,(GaussSeg*)gauss,this->interpolation);
}/*}}}*/
int  SegInput::GetResultArraySize(void){/*{{{*/
	return 1;
}
/*}}}*/
int  SegInput::GetResultInterpolation(void){/*{{{*/
	if(this->interpolation==P0Enum || this->interpolation==P0DGEnum){
		return P0Enum;
	}
	return P1Enum;
}/*}}}*/
int  SegInput::GetResultNumberOfNodes(void){/*{{{*/
	return SegRef::NumberofNodes(this->interpolation);
}
/*}}}*/
void SegInput::Scale(IssmDouble alpha){/*{{{*/

	for(int i=0;i<this->M*this->N;i++) this->values[i] = alpha*this->values[i];
	for(int i=0;i<SegRef::NumberofNodes(this->interpolation);i++) this->element_values[i] = alpha*this->element_values[i];
}
/*}}}*/
void SegInput::Pow(IssmDouble alpha){/*{{{*/

	for(int i=0;i<this->M*this->N;i++) this->values[i] = pow(this->values[i],alpha);
	for(int i=0;i<SegRef::NumberofNodes(this->interpolation);i++) this->element_values[i] = pow(this->element_values[i],alpha);
}
/*}}}*/
void SegInput::AXPY(Input* xinput,IssmDouble alpha){/*{{{*/

	/*xinput is of the same type, so cast it: */
	if(xinput->ObjectEnum()!=SegInputEnum) _error_("Operation not permitted because xinput is of type " << EnumToStringx(xinput->ObjectEnum()));
	SegInput* xseginput=xDynamicCast<SegInput*>(xinput);
	if(xseginput->GetInterpolation()!=this->interpolation) _error_("Operation not permitted because xinput is of type " << EnumToStringx(xinput->ObjectEnum()));

	/*Carry out the AXPY operation depending on type:*/
	for(int i=0;i<this->M*this->N;i++) this->values[i] = alpha*xseginput->values[i] + this->values[i];
	for(int i=0;i<SegRef::NumberofNodes(this->interpolation);i++) this->element_values[i] = alpha*xseginput->element_values[i] + this->element_values[i];
}
/*}}}*/
void SegInput::PointWiseMult(Input* xinput){/*{{{*/

	/*xinput is of the same type, so cast it: */
	if(xinput->ObjectEnum()!=SegInputEnum) _error_("Operation not permitted because xinput is of type " << EnumToStringx(xinput->ObjectEnum()));
	SegInput* xseginput=xDynamicCast<SegInput*>(xinput);
	if(xseginput->GetInterpolation()!=this->interpolation) _error_("Operation not permitted because xinput is of type " << EnumToStringx(xinput->ObjectEnum()));

	/* we need to check that the vector sizes are identical*/
	if(xseginput->M!=this->M||xseginput->N!=this->N) _error_("Operation not permitted because the inputs have different sizes");

	/*Carry out the AXPY operation depending on type:*/
	for(int i=0;i<this->M*this->N;i++) this->values[i] = xseginput->values[i] * this->values[i];
	for(int i=0;i<SegRef::NumberofNodes(this->interpolation);i++) this->element_values[i] = xseginput->element_values[i] * this->element_values[i];
}
/*}}}*/

/*Object functions*/

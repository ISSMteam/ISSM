/*!\file PentaInput.c
 * \brief: implementation of the PentaInput object
 */

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "../classes.h"
#include "../../shared/shared.h"
#include "./PentaInput.h"

/*PentaInput constructors and destructor*/
PentaInput::PentaInput(void){/*{{{*/

	this->numberofelements_local = -1;
	this->numberofvertices_local = -1;
	this->isserved       = false;
	this->isserved_collapsed= 0;
	this->M = -1;
	this->N = -1;
	this->values         = NULL;
	this->element_values = NULL;

}/*}}}*/
PentaInput::PentaInput(int nbe_in,int nbv_in,int interp_in){/*{{{*/

	_assert_(nbe_in>0);
	_assert_(nbe_in<1e11);
	_assert_(nbv_in>0);
	_assert_(nbv_in<1e11);
	this->numberofelements_local = nbe_in;
	this->numberofvertices_local = nbv_in;
	this->isserved           = false;
	this->isserved_collapsed = 0;

	/*Reset takes care of the rest*/
	this->Reset(interp_in);

}/*}}}*/
PentaInput::~PentaInput(){/*{{{*/
	if(this->element_values) xDelete<IssmDouble>(this->element_values);
	if(this->values)         xDelete<IssmDouble>(this->values);
}
/*}}}*/
void PentaInput::Reset(int interp_in){/*{{{*/

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
		this->N = PentaRef::NumberofNodes(interp_in);
	}

	/*Allocate Pointers*/
	this->values         = xNewZeroInit<IssmDouble>(this->M*this->N);
	this->element_values = xNewZeroInit<IssmDouble>(PentaRef::NumberofNodes(interp_in));
}/*}}}*/

/*Object virtual functions definitions:*/
Input* PentaInput::copy() {/*{{{*/

	/*Create output*/
	PentaInput* output = new PentaInput(this->numberofelements_local,this->numberofvertices_local,this->interpolation);

	/*Copy values*/
	xMemCpy<IssmDouble>(output->values,this->values,this->M*this->N);

	/*Return output*/
	return output;

}
/*}}}*/
void PentaInput::DeepEcho(void){/*{{{*/
	_printf_("PentaInput Echo:\n");
	_printf_("   interpolation:      "<<EnumToStringx(this->interpolation)<<"\n");
	_printf_("   nbe_local:          "<<this->numberofvertices_local<<"\n");
	_printf_("   nbv_local:          "<<this->numberofelements_local<<"\n");
	_printf_("   Size:               "<<M<<"x"<<N<<"\n");
	_printf_("   isserved:           "<<(isserved?"true":"false") << "\n");
	_printf_("   isserved_collapsed: "<<isserved_collapsed << "\n");
	if(isserved){
		_printf_("   current values:      ");
		if(isserved_collapsed){
			_printf_("[ ");
			for(int i=0;i<3;i++) _printf_(" "<<this->element_values[i]);
			_printf_("] ("<<EnumToStringx(this->interpolation)<<")\n");
		}
		else{
			_printf_("[ ");
			for(int i=0;i<PentaRef::NumberofNodes(this->interpolation);i++) _printf_(" "<<this->element_values[i]);
			_printf_("] ("<<EnumToStringx(this->interpolation)<<")\n");
		}
	}
}
/*}}}*/
void PentaInput::Echo(void){/*{{{*/
	_printf_(setw(15)<<"   PentaInput "<<setw(25)<<left<<EnumToStringx(-1));
	if(isserved){
		if(isserved_collapsed){
			_printf_("[ ");
			for(int i=0;i<3;i++) _printf_(" "<<this->element_values[i]);
			_printf_("] ("<<EnumToStringx(this->interpolation)<<")\n");
		}
		else{
			_printf_("[ ");
			for(int i=0;i<PentaRef::NumberofNodes(this->interpolation);i++) _printf_(" "<<this->element_values[i]);
			_printf_("] ("<<EnumToStringx(this->interpolation)<<")\n");
		}
	}
}
/*}}}*/
int  PentaInput::Id(void){/*{{{*/
	return -1;
}/*}}}*/
void PentaInput::Marshall(MarshallHandle* marshallhandle){ /*{{{*/

	int object_enum = PentaInputEnum;
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
		this->element_values = xNewZeroInit<IssmDouble>(PentaRef::NumberofNodes(this->interpolation));
	}
}
/*}}}*/
int  PentaInput::ObjectEnum(void){/*{{{*/
	return PentaInputEnum;
}
/*}}}*/

/*PentaInput management*/
void PentaInput::SetInput(int interp_in,int row,IssmDouble value_in){/*{{{*/

	_assert_(this);
	_assert_(row>=0);
	_assert_(row<this->M);
	_assert_(this->N==1);

	this->values[row] = value_in;
	this->isserved    = false;
}
/*}}}*/
void PentaInput::SetInput(int interp_in,int numindices,int* indices,IssmDouble* values_in){/*{{{*/

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
		_error_("not supported");
	}

	this->isserved    = false;
}
/*}}}*/
void PentaInput::SetInput(int interp_in,int row,int numindices,IssmDouble* values_in){/*{{{*/

	_assert_(this);
	if(interp_in==this->interpolation){
		_assert_(this->N==numindices);
	}
	else{
		this->Reset(interp_in);
		_assert_(this->N==numindices);
	}
	for(int i=0;i<numindices;i++) this->values[row*this->N+i] = values_in[i];

	this->isserved    = false;
}
/*}}}*/
void PentaInput::Serve(int numindices,int* indices){/*{{{*/

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
void PentaInput::Serve(int row,int numindices){/*{{{*/

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
}/*}}}*/
void PentaInput::ServeCollapsed(int row,int state){/*{{{*/

	_assert_(this);
	_assert_(this->N>=3);
	_assert_(row<this->M);
	_assert_(row>=0);

	if(state==1){
		for(int i=0;i<3;i++) this->element_values[i] = this->values[row*this->N+i];
		for(int i=3;i<6;i++) this->element_values[i] = 0.;
	}
	else if(state==2){
		for(int i=0;i<3;i++) this->element_values[i] = this->values[row*this->N+3+i];
		for(int i=3;i<6;i++) this->element_values[i] = 0.;
	}
	else{
		_error_("not supported");
	}

	/*Set input as served*/
	this->isserved = true;
	this->isserved_collapsed = state;
}/*}}}*/
void PentaInput::SetServeCollapsed(int state){/*{{{*/
	this->isserved_collapsed = state;
}/*}}}*/
int  PentaInput::GetInterpolation(){/*{{{*/
	return this->interpolation;
}/*}}}*/
void PentaInput::GetInputAverage(IssmDouble* pvalue){/*{{{*/
	_assert_(this);
	_assert_(this->isserved);

	/*Output*/
	IssmDouble value = 0.;

	if(this->isserved_collapsed){
		if(this->interpolation==P0Enum){
			value = this->element_values[0];
		}
		else{
			/*Assume P1...*/
			value = 1./3.*(this->element_values[0] +  this->element_values[1] +  this->element_values[2]);
		}
	}
	else{
		int        numnodes  = this->NumberofNodes(this->interpolation);
		IssmDouble numnodesd = reCast<int,IssmDouble>(numnodes);

		for(int i=0;i<numnodes;i++) value+=this->element_values[i];
		value = value/numnodesd;
	}

	*pvalue=value;
}/*}}}*/
IssmDouble PentaInput::GetInputMin(void){/*{{{*/
	_assert_(this);
	_assert_(this->isserved);

	int  numnodes  = this->NumberofNodes(this->interpolation);
	if(this->isserved_collapsed) numnodes = 3;
	IssmDouble min=this->element_values[0];

	for(int i=1;i<numnodes;i++){
		if(this->element_values[i]<min) min=this->element_values[i];
	}
	return min;
}/*}}}*/
IssmDouble PentaInput::GetInputMax(void){/*{{{*/
	_assert_(this);
	_assert_(this->isserved);

	int numnodes  = this->NumberofNodes(this->interpolation);
	if(this->isserved_collapsed) numnodes = 3;
	IssmDouble max=this->element_values[0];

	for(int i=1;i<numnodes;i++){
		if(this->element_values[i]>max) max=this->element_values[i];
	}
	return max;
}/*}}}*/
IssmDouble PentaInput::GetInputMaxAbs(void){/*{{{*/
	_assert_(this);
	_assert_(this->isserved);

	int numnodes  = this->NumberofNodes(this->interpolation);
	if(this->isserved_collapsed) numnodes = 3;
	IssmDouble maxabs=fabs(this->element_values[0]);

	for(int i=1;i<numnodes;i++){
		if(fabs(this->element_values[i])>maxabs) maxabs=fabs(this->element_values[i]);
	}
	return maxabs;
}/*}}}*/
void PentaInput::GetInputDerivativeValue(IssmDouble* derivativevalues, IssmDouble* xyz_list, Gauss* gauss){/*{{{*/
	_assert_(this);
	_assert_(this->isserved);
	if(this->isserved_collapsed){
		_assert_(gauss->Enum()==GaussTriaEnum);
		if(this->interpolation==P0Enum){
			derivativevalues[0] = 0.;
			derivativevalues[1] = 0.;
		}
		else{
			TriaRef temp;
			temp.GetInputDerivativeValue(derivativevalues,this->element_values,xyz_list,(GaussTria*)gauss,P1Enum);
		}
	}
	else{
		_assert_(gauss->Enum()==GaussPentaEnum);
		PentaRef::GetInputDerivativeValue(derivativevalues,this->element_values,xyz_list,(GaussPenta*)gauss,this->interpolation);
	}
}/*}}}*/
void PentaInput::GetInputValue(IssmDouble* pvalue,Gauss* gauss){/*{{{*/
	_assert_(this);
	_assert_(this->isserved);
	if(this->isserved_collapsed){
		_assert_(gauss->Enum()==GaussTriaEnum);
		if(this->interpolation==P0Enum){
			*pvalue = this->element_values[0];
		}
		else{
			TriaRef temp;
			temp.GetInputValue(pvalue,this->element_values,(GaussTria*)gauss,P1Enum);
		}
	}
	else{
		_assert_(gauss->Enum()==GaussPentaEnum);
		PentaRef::GetInputValue(pvalue,this->element_values,(GaussPenta*)gauss,this->interpolation);
	}
}/*}}}*/
int  PentaInput::GetResultArraySize(void){/*{{{*/
	return 1;
}
/*}}}*/
int  PentaInput::GetResultInterpolation(void){/*{{{*/
	if(this->interpolation==P0Enum || this->interpolation==P0DGEnum){
		return P0Enum;
	}
	return P1Enum;
}/*}}}*/
int  PentaInput::GetResultNumberOfNodes(void){/*{{{*/
	return PentaRef::NumberofNodes(this->interpolation);
}
/*}}}*/
void PentaInput::Scale(IssmDouble alpha){/*{{{*/

	for(int i=0;i<this->M*this->N;i++) this->values[i] = alpha*this->values[i];
	for(int i=0;i<PentaRef::NumberofNodes(this->interpolation);i++) this->element_values[i] = alpha*this->element_values[i];
}
/*}}}*/
void PentaInput::Pow(IssmDouble alpha){/*{{{*/

	for(int i=0;i<this->M*this->N;i++) this->values[i] = pow(this->values[i],alpha);
	for(int i=0;i<PentaRef::NumberofNodes(this->interpolation);i++) this->element_values[i] = pow(this->element_values[i],alpha);
}
/*}}}*/
void PentaInput::AXPY(Input* xinput,IssmDouble alpha){/*{{{*/

	/*xinput is of the same type, so cast it: */
	if(xinput->ObjectEnum()!=PentaInputEnum) _error_("Operation not permitted because xinput is of type " << EnumToStringx(xinput->ObjectEnum()));
	PentaInput* xpentainput=xDynamicCast<PentaInput*>(xinput);
	if(xpentainput->GetInterpolation()!=this->interpolation) _error_("Operation not permitted because xinput is of type " << EnumToStringx(xinput->ObjectEnum()));

	/*Carry out the AXPY operation depending on type:*/
	for(int i=0;i<this->M*this->N;i++) this->values[i] = alpha*xpentainput->values[i] + this->values[i];
	for(int i=0;i<PentaRef::NumberofNodes(this->interpolation);i++) this->element_values[i] = alpha*xpentainput->element_values[i] + this->element_values[i];
}
/*}}}*/
void PentaInput::PointWiseMult(Input* xinput){/*{{{*/

	/*xinput is of the same type, so cast it: */
	if(xinput->ObjectEnum()!=PentaInputEnum) _error_("Operation not permitted because xinput is of type " << EnumToStringx(xinput->ObjectEnum()));
	PentaInput* xpentainput=xDynamicCast<PentaInput*>(xinput);
	if(xpentainput->GetInterpolation()!=this->interpolation) _error_("Operation not permitted because xinput is of type " << EnumToStringx(xinput->ObjectEnum()));

	/* we need to check that the vector sizes are identical*/
	if(xpentainput->M!=this->M||xpentainput->N!=this->N) _error_("Operation not permitted because the inputs have different sizes");

	/*Carry out the pointwise operation depending on type:*/
	for(int i=0;i<this->M*this->N;i++) this->values[i] = xpentainput->values[i] * this->values[i];
	for(int i=0;i<PentaRef::NumberofNodes(this->interpolation);i++) this->element_values[i] = xpentainput->element_values[i] * this->element_values[i];
}
/*}}}*/

/*Object functions*/

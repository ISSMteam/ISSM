/*!\file ControlParam.c
 * \brief: implementation of the ControlParam object
 */

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "../classes.h"
#include "shared/shared.h"

/*ControlParam constructors and destructor*/
ControlParam::ControlParam(){/*{{{*/
	return;
}
/*}}}*/
ControlParam::ControlParam(IssmDouble* in_value, IssmDouble* in_minvalue, IssmDouble* in_maxvalue, int in_enum_type, int in_M,int in_N){/*{{{*/

	this->enum_type=in_enum_type;
	this->M=in_M;
	this->N=in_N;

	/*Sanity check, can't hurt*/
	if(this->N<1) _error_("Parameter is empty");
	if(this->M<1) _error_("Parameter is empty");
	if(this->M>2) _error_("Cannot handle more than 2 rows (as a TransientParam)");

	/*Assign value*/
	this->value=xNew<IssmDouble>(M*N);
	xMemCpy<IssmDouble>(value,in_value,M*N);

	/*Assign other fields*/
	this->minvalue=xNew<IssmDouble>(N);
	xMemCpy<IssmDouble>(minvalue,in_minvalue,N);
	this->maxvalue=xNew<IssmDouble>(N);
	xMemCpy<IssmDouble>(maxvalue,in_maxvalue,N);
	this->gradient=xNewZeroInit<IssmDouble>(N);
}
/*}}}*/
ControlParam::~ControlParam(){/*{{{*/
	xDelete<IssmDouble>(value);
	xDelete<IssmDouble>(minvalue);
	xDelete<IssmDouble>(maxvalue);
	xDelete<IssmDouble>(gradient);
	return;
}
/*}}}*/

/*Object virtual functions definitions:*/
Param* ControlParam::copy() {/*{{{*/

	ControlParam* output = new ControlParam();
	output->enum_type=this->enum_type;
	output->M=this->M;
	output->N=this->N;
	if(value){
		output->value=xNew<IssmDouble>(this->M*this->N);
		xMemCpy<IssmDouble>(output->value,this->value,this->M*this->N);
	}
	if(minvalue){
		output->minvalue=xNew<IssmDouble>(this->N);
		xMemCpy<IssmDouble>(output->minvalue,this->minvalue,this->N);
	}
	if(maxvalue){
		output->maxvalue=xNew<IssmDouble>(this->N);
		xMemCpy<IssmDouble>(output->maxvalue,this->maxvalue,this->N);
	}
	if(gradient){
		output->gradient=xNew<IssmDouble>(this->N);
		xMemCpy<IssmDouble>(output->gradient,this->gradient,this->N);
	}
	return output;

}
/*}}}*/
void ControlParam::DeepEcho(void){/*{{{*/

	_printf_(setw(22)<<"   ControlParam "<<setw(35)<<left<<EnumToStringx(this->enum_type)<<"\n ");
	if (value) _printf_("---value: ");
	for(int i=0;i<this->M;i++) _printf_(" "<< this->value[i]);
	_printf_("]\n");
	if (minvalue) _printf_("---minvalue: ");
	for(int i=0;i<this->M;i++) _printf_(" "<< this->minvalue[i]);
	_printf_("]\n");
	if (maxvalue) _printf_("---maxvalue: ");
	for(int i=0;i<this->M;i++) _printf_(" "<< this->maxvalue[i]);
	_printf_("]\n");
	if (gradient) _printf_("---gradient: " << this->gradient << "\n");
}
/*}}}*/
void ControlParam::Echo(void){/*{{{*/
	this->DeepEcho();
}
/*}}}*/
int  ControlParam::Id(void){ return -1; }/*{{{*/
/*}}}*/
void ControlParam::Marshall(MarshallHandle* marshallhandle){ /*{{{*/

	int object_enum = ControlParamEnum;
   marshallhandle->call(object_enum);
	marshallhandle->call(this->enum_type);
	marshallhandle->call(this->M);
	marshallhandle->call(this->N);
	marshallhandle->call(this->value,this->M*this->N);
	marshallhandle->call(this->minvalue,this->N);
	marshallhandle->call(this->maxvalue,this->N);
	marshallhandle->call(this->gradient,this->N);

}
/*}}}*/
int  ControlParam::ObjectEnum(void){/*{{{*/

	return ControlParamEnum;

}
/*}}}*/

void  ControlParam::GetParameterValue(IssmDouble** poutput,int* pN, const char* data){/*{{{*/

	IssmDouble* output=xNew<IssmDouble>(N);

	if(strcmp(data,"value")==0){
		xMemCpy<IssmDouble>(output,value,N);
	}
	else if (strcmp(data,"lowerbound")==0){
		xMemCpy<IssmDouble>(output,minvalue,N);
	}
	else if (strcmp(data,"upperbound")==0){
		xMemCpy<IssmDouble>(output,maxvalue,N);
	}
	else if (strcmp(data,"gradient")==0){
		xMemCpy<IssmDouble>(output,gradient,N);
	}
	else{
		_error_("Data " << data << " not supported yet");
	}

	/*Assign output pointers:*/
	if(pN) *pN=N;
	*poutput=output;
}
/*}}}*/
void  ControlParam::GetParameterValue(IssmDouble* poutput){/*{{{*/

	/*Copy entire vector if M==1, or first row if M==2*/
	if(M==1){
		xMemCpy<IssmDouble>(poutput,value,N);
		return;
	}

	_error_("STOP");

}
/*}}}*/
void  ControlParam::GetParameterValue(IssmDouble* poutput, IssmDouble time){/*{{{*/

	if(M==1){
		*poutput = value[0];
		return;
	}

	IssmDouble *timesteps = &this->value[1*this->N+0];
	IssmDouble output;
	bool       found;

	/*Ok, we have the time, go through the timesteps, and figure out which interval we 
	 *fall within. Then interpolate the values on this interval: */
	if(time<timesteps[0]){
		/*get values for the first time: */
		output=this->value[0];
		found=true;
	}
	else if(time>timesteps[this->N-1]){
		/*get values for the last time: */
		output=this->value[this->N-1];
		found=true;
	}
	else{
		/*Find which interval we fall within: */
		for(int i=0;i<this->N;i++){
			if(time==timesteps[i]){
				/*We are right on one step time: */
				output=this->value[i];
				found=true;
				break; //we are done with the time interpolation.
			}
			else{
				if(timesteps[i]<time && time<timesteps[i+1]){
					/*ok, we have the interval ]i:i+1[. Interpolate linearly for now: */
					IssmDouble deltat=timesteps[i+1]-timesteps[i];
					IssmDouble alpha=(time-timesteps[i])/deltat;
					output=(1.0-alpha)*this->value[i] + alpha*this->value[i+1];
					found=true;
					break;
				}
				else continue; //keep looking on the next interval
			}
		}
	}
	if(!found)_error_("did not find time interval on which to interpolate values");
	//_printf_("for time = "<<time/31536000.<<" yr, melt = "<<output*31536000.<<" m/yr\n");

	*poutput=output;
}
/*}}}*/
void  ControlParam::GetParameterValue(IssmDouble** poutput, int* pN){/*{{{*/

	/*This method should be specific to VectorParams, only one tow required*/
	_assert_(N>0);
	_assert_(M==1);
	IssmDouble* output=xNew<IssmDouble>(N);
	xMemCpy<IssmDouble>(output,value,N);

	/*Assign output pointers:*/
	if(pN) *pN=N;
	*poutput=output;
}
/*}}}*/
void  ControlParam::SetValue(IssmDouble* poutput,int in_M, int in_N){/*{{{*/

	_assert_(in_N==this->N);
	_assert_(in_M==1);
	xMemCpy<IssmDouble>(this->value,poutput,in_N);
}
/*}}}*/
void  ControlParam::SetGradient(IssmDouble* poutput,int in_M, int in_N){/*{{{*/

	_assert_(in_M==1);
	xMemCpy<IssmDouble>(this->gradient,poutput,in_N);
}
/*}}}*/
void  ControlParam::GetVectorFromControl(Vector<IssmDouble>* vector,int control_index,int in_N,const char* data,int offset){/*{{{*/

	/*Get list of ids for this element and this control*/
	_assert_(in_N==this->N);
	int* idlist = xNew<int>(this->N);
	for(int i=0;i<this->N;i++) idlist[i] = offset+i;

	/*Get data*/
	IssmDouble* values = NULL;
	GetParameterValue(&values, NULL, data);

	/*Enter data in vector*/
	vector->SetValues(this->N,idlist,values,INS_VAL);

	/*Clean up*/
	xDelete<int>(idlist);
	xDelete<IssmDouble>(values);

}/*}}}*/

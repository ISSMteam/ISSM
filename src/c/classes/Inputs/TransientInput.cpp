/*!\file TransientInput.c
 * \brief: implementation of the TransientInput object
 */
/*Headers*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include <numeric>
#include "./TransientInput.h"
#include "./TriaInput.h"
#include "./PentaInput.h"
#include "../../shared/shared.h"
#include "../Params/Parameters.h"

/*TransientInput constructors and destructor*/
TransientInput::TransientInput(){/*{{{*/

	this->enum_type=UNDEF;
	this->inputs=NULL;
	this->numtimesteps=0;
	this->parameters=NULL;
	this->timesteps=NULL;

	this->current_input=NULL;
	this->current_step=-1;

}
/*}}}*/
TransientInput::TransientInput(int in_enum_type,int nbe,int nbv,IssmDouble* timesin,int N){/*{{{*/

	/*Set Enum*/
	this->enum_type=in_enum_type;
	this->numberofelements_local = nbe;
	this->numberofvertices_local = nbv;

	/*Allocate values and timesteps, and copy: */
	_assert_(N>=0 && N<1e6);
	this->numtimesteps=N;
	if(N>0){
		this->timesteps=xNew<IssmDouble>(N);
		xMemCpy(this->timesteps,timesin,N);

		this->inputs     = xNew<Input*>(N);
		for(int i=0;i<N;i++) this->inputs[i] = NULL;
	}
	else{
		this->timesteps=0;
		this->inputs   =0;
	}
	this->parameters = NULL;
	this->current_input=NULL;
	this->current_step=-1;
}
/*}}}*/
TransientInput::~TransientInput(){/*{{{*/

	for(int i=0;i<this->numtimesteps;i++){
		delete this->inputs[i];
	}
	xDelete<Input*>(this->inputs);
	xDelete<IssmDouble>(this->timesteps);

	if(this->current_input) delete this->current_input;
}
/*}}}*/

/*Object virtual functions definitions:*/
Input* TransientInput::copy() {/*{{{*/

	TransientInput* output=NULL;

	output = new TransientInput();
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

	return output;
}/*}}}*/
void TransientInput::DeepEcho(void){/*{{{*/

	_printf_("TransientInput:\n");
	_printf_("   enum: " << this->enum_type << " (" << EnumToStringx(this->enum_type) << ")\n");
	_printf_("   numtimesteps: " << this->numtimesteps << "\n");
	_printf_("---inputs: \n");
	for(int i=0;i<this->numtimesteps;i++){
		_printf_("   time: " << this->timesteps[i]<<"  ");
		if(this->inputs[i]) this->inputs[i]->DeepEcho();
		else                _printf_(" NOT SET! \n");
	}
}
/*}}}*/
void TransientInput::Configure(Parameters* params){/*{{{*/
	this->parameters=params;
}
/*}}}*/
void TransientInput::Echo(void){/*{{{*/
	_printf_("TransientInput:\n");
	_printf_("   enum: " << this->enum_type << " (" << EnumToStringx(this->enum_type) << ")\n");
	_printf_("   numtimesteps: " << this->numtimesteps << "\n");
	_printf_("---inputs: \n");
	for(int i=0;i<this->numtimesteps;i++){
		_printf_("   time: " << this->timesteps[i]<<"  ");
		if(this->inputs[i]) this->inputs[i]->Echo();
		else                _printf_(" NOT SET! \n");
	}
}
/*}}}*/
int  TransientInput::Id(void){ return -1; }/*{{{*/
/*}}}*/
void TransientInput::Marshall(MarshallHandle* marshallhandle){ /*{{{*/

	bool isnull;

	int object_enum = TransientInputEnum;
   marshallhandle->call(object_enum);

	marshallhandle->call(this->numberofelements_local);
	marshallhandle->call(this->numberofvertices_local);
	marshallhandle->call(this->enum_type);
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
int  TransientInput::ObjectEnum(void){/*{{{*/

	return TransientInputEnum;

}
/*}}}*/

/*Intermediary*/
void TransientInput::AddTriaTimeInput(IssmDouble time,int numindices,int* indices,IssmDouble* values_in,int interp_in){/*{{{*/

	/*Check whether this is the last time step that we have*/
	if(this->numtimesteps){
		if(fabs(this->timesteps[this->numtimesteps-1]-time)<1.0e-5){
			this->AddTriaTimeInput(this->numtimesteps-1,numindices,indices,values_in,interp_in);
			return;
		}
	}

	/*This is a new time step! we need to add it to the list*/
	if(this->numtimesteps>0 && time<this->timesteps[this->numtimesteps-1]) _error_("timestep values must increase sequentially, here " << this->timesteps[this->numtimesteps-1] <<" is the last step but smaller than the preceding "<< time<<"\n");

	IssmDouble *old_timesteps = NULL;
	Input    **old_inputs    = NULL;
	if (this->numtimesteps > 0){
		old_timesteps=xNew<IssmDouble>(this->numtimesteps);
		xMemCpy(old_timesteps,this->timesteps,this->numtimesteps);
		xDelete<IssmDouble>(this->timesteps);
		old_inputs=xNew<Input*>(this->numtimesteps);
		xMemCpy(old_inputs,this->inputs,this->numtimesteps);
		xDelete<Input*>(this->inputs);
	}

	this->numtimesteps=this->numtimesteps+1;
	this->timesteps=xNew<IssmDouble>(this->numtimesteps);
	this->inputs   = xNew<Input*>(this->numtimesteps);

	if (this->numtimesteps > 1){
		xMemCpy(this->inputs,old_inputs,this->numtimesteps-1);
		xMemCpy(this->timesteps,old_timesteps,this->numtimesteps-1);
		xDelete(old_timesteps);
		xDelete<Input*>(old_inputs);
	}

	/*go ahead and plug: */
	this->timesteps[this->numtimesteps-1] = time;
	this->inputs[this->numtimesteps-1]    = NULL;
	this->AddTriaTimeInput(this->numtimesteps-1,numindices,indices,values_in,interp_in);

}
/*}}}*/
void TransientInput::AddPentaTimeInput(IssmDouble time,int numindices,int* indices,IssmDouble* values_in,int interp_in){/*{{{*/

	/*Check whether this is the last time step that we have*/
	if(this->numtimesteps){
		if(fabs(this->timesteps[this->numtimesteps-1]-time)<1.0e-5){
			this->AddPentaTimeInput(this->numtimesteps-1,numindices,indices,values_in,interp_in);
			return;
		}
	}

	/*This is a new time step! we need to add it to the list*/
	if(this->numtimesteps>0 && time<this->timesteps[this->numtimesteps-1]) _error_("timestep values must increase sequentially");

	IssmDouble *old_timesteps = NULL;
	Input    **old_inputs    = NULL;
	if (this->numtimesteps > 0){
		old_timesteps=xNew<IssmDouble>(this->numtimesteps);
		xMemCpy(old_timesteps,this->timesteps,this->numtimesteps);
		xDelete<IssmDouble>(this->timesteps);
		old_inputs=xNew<Input*>(this->numtimesteps);
		xMemCpy(old_inputs,this->inputs,this->numtimesteps);
		xDelete<Input*>(this->inputs);
	}

	this->numtimesteps=this->numtimesteps+1;
	this->timesteps=xNew<IssmDouble>(this->numtimesteps);
	this->inputs   = xNew<Input*>(this->numtimesteps);

	if (this->numtimesteps > 1){
		xMemCpy(this->inputs,old_inputs,this->numtimesteps-1);
		xMemCpy(this->timesteps,old_timesteps,this->numtimesteps-1);
		xDelete(old_timesteps);
		xDelete<Input*>(old_inputs);
	}

	/*go ahead and plug: */
	this->timesteps[this->numtimesteps-1] = time;
	this->inputs[this->numtimesteps-1]    = NULL;
	this->AddPentaTimeInput(this->numtimesteps-1,numindices,indices,values_in,interp_in);

}
/*}}}*/
void TransientInput::AddTriaTimeInput(int step,int numindices,int* indices,IssmDouble* values_in,int interp_in){/*{{{*/

	_assert_(step>=0 && step<this->numtimesteps);

	/*Create it if necessary*/
	if(this->inputs[step]){
		if(this->inputs[step]->ObjectEnum()!=TriaInputEnum) _error_("cannot add Element values to a "<<EnumToStringx(this->inputs[step]->ObjectEnum()));
	}
	else{
		this->inputs[step] = new TriaInput(this->numberofelements_local,this->numberofvertices_local,interp_in);
	}

	/*Set input*/
	TriaInput* input = xDynamicCast<TriaInput*>(this->inputs[step]);
	input->SetInput(interp_in,numindices,indices,values_in);

}
/*}}}*/
void TransientInput::AddPentaTimeInput(int step,int numindices,int* indices,IssmDouble* values_in,int interp_in){/*{{{*/

	_assert_(step>=0 && step<this->numtimesteps);

	/*Create it if necessary*/
	if(this->inputs[step]){
		if(this->inputs[step]->ObjectEnum()!=PentaInputEnum) _error_("cannot add Element values to a "<<EnumToStringx(this->inputs[step]->ObjectEnum()));
	}
	else{
		this->inputs[step] = new PentaInput(this->numberofelements_local,this->numberofvertices_local,interp_in);
	}

	/*Set input*/
	PentaInput* input = xDynamicCast<PentaInput*>(this->inputs[step]);
	input->SetInput(interp_in,numindices,indices,values_in);

}
/*}}}*/
void TransientInput::GetAllTimes(IssmDouble** ptimesteps,int* pnumtimesteps){/*{{{*/

	if(ptimesteps){
		*ptimesteps=xNew<IssmDouble>(this->numtimesteps);
		xMemCpy(*ptimesteps,this->timesteps,this->numtimesteps);
	}
	if(pnumtimesteps){
		*pnumtimesteps = this->numtimesteps;
	}

}
/*}}}*/
TriaInput* TransientInput::GetTriaInput(){/*{{{*/

	IssmDouble time;
	this->parameters->FindParam(&time,TimeEnum);
	return this->GetTriaInput(time);

}
/*}}}*/
TriaInput* TransientInput::GetTriaInput(IssmDouble time){/*{{{*/

	/*Set current time input*/
	this->SetCurrentTimeInput(time);
	_assert_(this->current_input);

	/*Cast and return*/
	if(this->current_input->ObjectEnum()!=TriaInputEnum){
		_error_("Cannot return a TriaInput");
	}
	return xDynamicCast<TriaInput*>(this->current_input);

}
/*}}}*/
TriaInput* TransientInput::GetTriaInput(IssmDouble start_time, IssmDouble end_time, int averaging_method){/*{{{*/

	/*Set current time input*/
	this->SetAverageAsCurrentTimeInput(start_time,end_time,averaging_method);
	_assert_(this->current_input);

	/*Cast and return*/
	if(this->current_input->ObjectEnum()!=TriaInputEnum){
		_error_("Cannot return a TriaInput");
	}
	return xDynamicCast<TriaInput*>(this->current_input);

}
/*}}}*/
TriaInput* TransientInput::GetTriaInput(int offset){/*{{{*/

	/*Check offset*/
	if(offset<0 || offset>this->numtimesteps-1){
		_error_("Cannot return input for offset "<<offset);
	}
	Input* input = this->inputs[offset];

	/*Cast and return*/
	_assert_(input);
	if(input->ObjectEnum()!=TriaInputEnum) _error_("Cannot return a TriaInput");
	return xDynamicCast<TriaInput*>(input);

}
/*}}}*/
PentaInput* TransientInput::GetPentaInput(){/*{{{*/

	IssmDouble time;
	this->parameters->FindParam(&time,TimeEnum);
	return this->GetPentaInput(time);
}
/*}}}*/
PentaInput* TransientInput::GetPentaInput(IssmDouble time){/*{{{*/

	/*Set current time input*/
	this->SetCurrentTimeInput(time);
	_assert_(this->current_input);

	/*Cast and return*/
	if(this->current_input->ObjectEnum()!=PentaInputEnum){
		_error_("Cannot return a PentaInput");
	}
	return xDynamicCast<PentaInput*>(this->current_input);

}
/*}}}*/
PentaInput* TransientInput::GetPentaInput(int offset){/*{{{*/

	/*Check offset*/
	if(offset<0 || offset>this->numtimesteps-1){
		_error_("Cannot return input for offset "<<offset);
	}
	Input* input = this->inputs[offset];

	/*Cast and return*/
	if(input->ObjectEnum()!=PentaInputEnum) _error_("Cannot return a PentaInput");
	return xDynamicCast<PentaInput*>(input);

}
/*}}}*/
PentaInput* TransientInput::GetPentaInput(IssmDouble start_time, IssmDouble end_time, int averaging_method){/*{{{*/

	/*Set current time input*/
	this->SetAverageAsCurrentTimeInput(start_time,end_time,averaging_method);
	_assert_(this->current_input);

	/*Cast and return*/
	if(this->current_input->ObjectEnum()!=PentaInputEnum){
		_error_("Cannot return a PentaInput");
	}
	return xDynamicCast<PentaInput*>(this->current_input);

}
/*}}}*/

void TransientInput::SetCurrentTimeInput(IssmDouble time){/*{{{*/

	/*First, recover current time from parameters: */
	bool linear_interp,average,cycle;
	int  timestepping;
	IssmDouble dt;
	this->parameters->FindParam(&linear_interp,TimesteppingInterpForcingEnum);
	this->parameters->FindParam(&average,TimesteppingAverageForcingEnum);
	this->parameters->FindParam(&cycle,TimesteppingCycleForcingEnum);
	this->parameters->FindParam(&dt,TimesteppingTimeStepEnum);          /*transient core time step*/
	this->parameters->FindParam(&timestepping,TimesteppingTypeEnum);

	if(cycle){

		/*Change input time if we cycle through the forcing*/
		IssmDouble time0 = this->timesteps[0];
		IssmDouble time1 = this->timesteps[this->numtimesteps - 1];

		if(timestepping!=AdaptiveTimesteppingEnum){
			/*We need the end time to be the last timestep that would be taken*/
			/* i.e., the case where GEMB has time stamps (finer timestep) after the last timestep */
			/* warning: this assumes dt = constant!!*/
			IssmDouble nsteps = reCast<int,IssmDouble>(time1/dt);
			if (reCast<IssmDouble>(nsteps)<time1/dt) nsteps=nsteps+1;
			time1 = nsteps*dt;
		}

		/*See by how many intervals we have to offset time*/
		IssmDouble deltat = time1-time0;

		//int num_intervals = floor((time-time0)/deltat); //Cannot do that because of AD!
		int num_intervals = reCast<int,IssmDouble>(fabs(time-time0)/deltat);

		/*Uncomment following line if you would like to apply a cycle BEFORE the time series starts*/
		if(time<time0) num_intervals = -num_intervals-1;

		if(fabs(time-time0)/deltat == reCast<IssmDouble>(num_intervals)){
			/*Hack to make sure we always cover the last value of the series (discussion with Nicole)*/
			time = time1;
		}
		else{
			/*Now offset time so that we do the right interpolation below*/
			time = time - num_intervals*deltat;
		}
	}

	/*Figure step out*/
	int offset, prevoffset;
	if(!binary_search(&offset,time,this->timesteps,this->numtimesteps)){
		_error_("Input not found (is TransientInput sorted ?)");
	}
	if(!binary_search(&prevoffset,reCast<IssmDouble>(time-dt),this->timesteps,this->numtimesteps)){
		_error_("Input not found (is TransientInput sorted ?)");
	}

	if (offset==-1){

		/*get values for the first time: */
		_assert_(time<this->timesteps[0]);

		/*If already processed return*/
		if(this->current_step==0.) return;

		/*Prepare input*/
		if(this->current_input) delete this->current_input;
		this->current_step = 0.;
		this->current_input = this->inputs[0]->copy();

	}
	else if(offset-prevoffset>1 && prevoffset >=0 && average){
		/*get values for the last time: */
		_assert_(time>=this->timesteps[offset]);

		/*If already processed return*/
		if(this->current_step==reCast<IssmDouble>(offset)) return;

		/*Prepare input*/
		if(this->current_input) delete this->current_input;
		this->current_step  = reCast<IssmDouble>(offset);

		this->current_input = this->inputs[prevoffset]->copy();
		for(int i=prevoffset+1;i<offset;i++) {
			this->current_input->AXPY(this->inputs[i],+1.0);
		}
		this->current_input->Scale(1./(offset-prevoffset));

	}
	else if(offset==(this->numtimesteps-1) || !linear_interp){

		/*get values for the last time: */
		_assert_(time>=this->timesteps[offset]);

		/*If already processed return*/
		if(this->current_step==reCast<IssmDouble>(offset)) return;

		/*Prepare input*/
		if(this->current_input) delete this->current_input;
		this->current_step  = reCast<IssmDouble>(offset);
		this->current_input = this->inputs[offset]->copy();
	}
	else {

		/*Interpolate */
		_assert_(time>=this->timesteps[offset] && time<this->timesteps[offset+1]);

		/*get values between two times [offset:offset+1[, Interpolate linearly*/
		IssmDouble deltat=this->timesteps[offset+1]-this->timesteps[offset];
		IssmDouble this_step = reCast<IssmDouble>(offset) + (time - this->timesteps[offset])/deltat;

		/*If already processed return*/
		if(fabs(this->current_step-this_step)<1.e-5) return;

		/*Prepare input*/
		if(this->current_input) delete this->current_input;
		this->current_step = this_step;
		IssmDouble alpha2=(time-this->timesteps[offset])/deltat;
		IssmDouble alpha1=(1.0-alpha2);

		Input* input1=this->inputs[offset];
		Input* input2=this->inputs[offset+1];

		this->current_input = input1->copy();
		this->current_input->Scale(alpha1);
		this->current_input->AXPY(input2,alpha2);
	}
}/*}}}*/
void TransientInput::SetAverageAsCurrentTimeInput(IssmDouble start_time,IssmDouble end_time, int averaging_method){/*{{{*/

	IssmDouble  dt,durinv;
	IssmDouble  dtsum=0;
	IssmDouble  timespan,mid_step;
	int         found,start_offset,end_offset,input_offset;

	/*go through the timesteps, and grab offset for start and end*/
	found=binary_search(&start_offset,start_time,this->timesteps,this->numtimesteps);
	if(!found) _error_("Input not found (is TransientInput sorted ?)");
	found=binary_search(&end_offset,end_time,this->timesteps,this->numtimesteps);
	if(!found) _error_("Input not found (is TransientInput sorted ?)");

	if(start_offset==-1){
		timespan=this->timesteps[end_offset]-start_time;
	}
	else{
		timespan=this->timesteps[end_offset]-this->timesteps[start_offset];
	}
	mid_step=reCast<IssmDouble>(start_offset)+0.5*timespan;
	/*If already processed return, we set step in the middle of the interval*/
	if(fabs(this->current_step-mid_step)<1.e-5) return;
	/*If not processed set current_step*/
	if(this->current_input) delete this->current_input;
	this->current_step = mid_step;

	int offset=start_offset;
	while(offset < end_offset){
		if(offset==start_offset){
			dt=this->timesteps[offset+1]-start_time;
			_assert_(dt>0.);
			if(offset==end_offset-1){
				dt=end_time-start_time;
				_assert_(dt>0.);
			}
		}
		else if(offset==end_offset-1){
			dt=end_time-this->timesteps[offset];
			_assert_(dt>0.);
		}
		else{
			dt=this->timesteps[offset+1]-this->timesteps[offset];
			_assert_(dt>0.);
		}
		Input* stepinput=this->inputs[offset+1]->copy();
		switch(averaging_method){
			case 0: /*Arithmetic mean*/
				if(offset==start_offset){
					this->current_input=stepinput->copy();
					this->current_input->Scale(dt);
				}
				else{
					this->current_input->AXPY(stepinput,dt);
				}
				break;
			case 1: /*Geometric mean*/
				if(offset==start_offset){
					this->current_input = stepinput->copy();
					this->current_input->Scale(dt);
				}
				else{
					stepinput->Scale(dt);
					this->current_input->PointWiseMult(stepinput);
				}
				break;
			case 2: /*Harmonic mean*/
				if(offset==start_offset){
					this->current_input = stepinput->copy();
					this->current_input->Pow(-1);
					this->current_input->Scale(dt);
				}
				else{
					stepinput->Pow(-1);
					this->current_input->AXPY(stepinput,dt);
				}
			default:
				_error_("averaging method is not recognised");
		}
		delete stepinput;
		dtsum+=dt;
		offset+=1;
	}
	_assert_(dtsum>0);
	durinv=1./dtsum;
	/*Integration done, now normalize*/
	switch(averaging_method){
		case 0: //Arithmetic mean
			this->current_input->Scale(durinv);
			break;
		case 1: /*Geometric mean*/
			this->current_input->Pow(durinv);
			break;
		case 2: /*Harmonic mean*/
			this->current_input->Scale(durinv);
			this->current_input->Pow(-1);
		default:
			_error_("averaging method is not recognised");
	}
}/*}}}*/
IssmDouble  TransientInput::GetTimeByOffset(int offset){/*{{{*/
	if(offset<0) offset=0;
	_assert_(offset<this->numtimesteps);
	return this->timesteps[offset];
}
/*}}}*/
int  TransientInput::GetTimeInputOffset(IssmDouble time){/*{{{*/

	int offset;

	/*go through the timesteps, and figure out which interval we
	 *     *fall within. Then interpolate the values on this interval: */
	int found=binary_search(&offset,time,this->timesteps,this->numtimesteps);
	if(!found) _error_("Input not found (is TransientInput sorted ?)");

	return offset;
}
/*}}}*/

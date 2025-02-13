/*!\file TransientParam.c
 * \brief: implementation of the TransientParam object
 */

/*header files: */
/*{{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "../classes.h"
#include "../../shared/shared.h"
/*}}}*/

/*TransientParam constructors and destructor*/
TransientParam::TransientParam(){/*{{{*/
	return;
}
/*}}}*/
TransientParam::TransientParam(int in_enum_type,IssmDouble* in_values,IssmDouble* in_time,bool interpolation_on,bool cycle_in,int in_N){/*{{{*/

	_assert_(in_values && in_time);

	this->enum_type=in_enum_type;
	this->N=in_N;
	this->interpolation=interpolation_on;
	this->cycle=cycle_in;

	this->values=xNew<IssmDouble>(N);
	xMemCpy<IssmDouble>(this->values,in_values,N);

	this->timesteps=xNew<IssmDouble>(N);
	xMemCpy<IssmDouble>(this->timesteps,in_time,N);
}
/*}}}*/
TransientParam::~TransientParam(){/*{{{*/
	xDelete<IssmDouble>(values);
	xDelete<IssmDouble>(timesteps);
}
/*}}}*/

/*Object virtual functions definitions:*/
Param* TransientParam::copy() {/*{{{*/

	return new TransientParam(this->enum_type,this->values,this->timesteps,this->interpolation,this->cycle,this->N);

}
/*}}}*/
void TransientParam::DeepEcho(void){/*{{{*/

	_printf_("TransientParam:\n");
	_printf_("   enum: " << this->enum_type << " (" << EnumToStringx(this->enum_type) << ")\n");
	_printf_("   size: " << this->N << "\n");
	for(int i=0;i<this->N;i++){
		_printf_(   "time: " << this->timesteps[i] << " value: " << this->values[i] << "\n");
	}
}
/*}}}*/
void TransientParam::Echo(void){/*{{{*/

	_printf_("TransientParam:\n");
	_printf_("   enum: " << this->enum_type << " (" << EnumToStringx(this->enum_type) << ")\n");
	_printf_("   size: " << this->N << "\n");

}
/*}}}*/
int  TransientParam::Id(void){ return -1; }/*{{{*/
/*}}}*/
void TransientParam::Marshall(MarshallHandle* marshallhandle){ /*{{{*/

	int object_enum = TransientParamEnum;
   marshallhandle->call(object_enum);

	marshallhandle->call(this->enum_type);
	marshallhandle->call(this->interpolation);
	marshallhandle->call(this->cycle);
	marshallhandle->call(this->N);
	if(marshallhandle->OperationNumber()==MARSHALLING_LOAD){
		values=xNew<IssmDouble>(N);
		timesteps=xNew<IssmDouble>(N);
	}
	marshallhandle->call(this->values,N);
	marshallhandle->call(this->timesteps,N);

}
/*}}}*/
int  TransientParam::ObjectEnum(void){/*{{{*/

	return TransientParamEnum;

}
/*}}}*/

/*TransientParam virtual functions definitions: */
void  TransientParam::GetParameterValue(IssmDouble* pdouble,IssmDouble time){/*{{{*/

	IssmDouble output;
	bool   found;

	if(this->cycle) _error_("not implemented yet");

	/*Ok, we have the time, go through the timesteps, and figure out which interval we 
	 *fall within. Then interpolate the values on this interval: */
	if(time<this->timesteps[0]){
		/*get values for the first time: */
		output=this->values[0];
		found=true;
	}
	else if(time>this->timesteps[this->N-1]){
		/*get values for the last time: */
		output=this->values[this->N-1];
		found=true;
	}
	else{
		/*Find which interval we fall within: */
		for(int i=0;i<this->N;i++){
			if(time==this->timesteps[i]){
				/*We are right on one step time: */
				output=this->values[i];
				found=true;
				break; //we are done with the time interpolation.
			}
			else{
				if(this->timesteps[i]<time && time<this->timesteps[i+1]){
					/*ok, we have the interval ]i:i+1[. Interpolate linearly for now: */
					IssmDouble deltat=this->timesteps[i+1]-this->timesteps[i];
					IssmDouble alpha=(time-this->timesteps[i])/deltat;
					if(this->interpolation==true) output=(1.0-alpha)*this->values[i] + alpha*this->values[i+1];
					else output=this->values[i];
					found=true;
					break;
				}
				else continue; //keep looking on the next interval
			}
		}
	}
	if(!found)_error_("did not find time interval on which to interpolate values");
	*pdouble=output;
}
/*}}}*/
void  TransientParam::GetParameterValue(IssmDouble* pdouble,IssmDouble time, int timestepping, IssmDouble dt){/*{{{*/

	bool cycle = this->cycle;

	if(cycle){

		/*Change input time if we cycle through the forcing*/
		IssmDouble time0 = this->timesteps[0];
		IssmDouble time1 = this->timesteps[this->N - 1];

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

		this->cycle = false;
	}

	this->GetParameterValue(pdouble,time);

	this->cycle = cycle;
}
/*}}}*/


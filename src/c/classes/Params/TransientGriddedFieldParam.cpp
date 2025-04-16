/*!\file TransientGriddedFieldParam.c
 * \brief: implementation of the TransientGriddedFieldParam object
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

/*TransientGriddedFieldParam constructors and destructor*/
TransientGriddedFieldParam::TransientGriddedFieldParam(){/*{{{*/
	return;
}
/*}}}*/
TransientGriddedFieldParam::TransientGriddedFieldParam(int in_enum_type,IssmDouble* in_values,IssmDouble* in_time,bool interpolation_on,bool cycle_in,int in_N,int in_M, int in_T){/*{{{*/

	_assert_(in_values && in_time);

	this->enum_type=in_enum_type;
	this->M=in_M; //Number of rows (LAT)
	this->N=in_N; //Number of columns (LON)
	this->MN=this->M*this->N;; //Number of data per time step
	this->T=in_T; //Number of time steps
	this->interpolation=interpolation_on;
	this->cycle=cycle_in;

	this->values=xNew<IssmDouble>(M*N*T);
	xMemCpy<IssmDouble>(this->values,in_values,M*N*T);

	this->timesteps=xNew<IssmDouble>(T);
	xMemCpy<IssmDouble>(this->timesteps,in_time,T);
}
/*}}}*/
TransientGriddedFieldParam::~TransientGriddedFieldParam(){/*{{{*/
	xDelete<IssmDouble>(values);
	xDelete<IssmDouble>(timesteps);
}/*}}}*/

/*Object virtual functions definitions:*/
Param* TransientGriddedFieldParam::copy() {/*{{{*/

	return new TransientGriddedFieldParam(this->enum_type,this->values,this->timesteps,this->interpolation,this->cycle,this->M,this->N, this->T);

}
/*}}}*/
void TransientGriddedFieldParam::DeepEcho(void){/*{{{*/

	_printf_("TransientGriddedFieldParam:\n");
	_printf_("   enum: " << this->enum_type << " (" << EnumToStringx(this->enum_type) << ")\n");
	_printf_("   number of time steps: " << this->T << "\n");
	_printf_("   number of rows (LAT): " << this->M << "\n");
	_printf_("   number of columns (LON): " << this->N << "\n");
	for(int t=0;t<this->T;t++){
		_printf_("	time: " << this->timesteps[t] << "\n");
		for(int i=0;i<this->N;i++){
			for(int k=0;k<this->M;k++){
				_printf_("		values: " << this->values[t*this->MN+k*this->N+i] << "\n");
			}
			_printf_("\n");
		}
	}
}
/*}}}*/
void TransientGriddedFieldParam::Echo(void){/*{{{*/

	_printf_("TransientGriddedFieldParam:\n");
	_printf_("   enum: " << this->enum_type << " (" << EnumToStringx(this->enum_type) << ")\n");
	_printf_("   size: " << this->M << " by " << this->N <<  " by " << this->T << "\n");

}
/*}}}*/
int  TransientGriddedFieldParam::Id(void){ return -1; }/*{{{*/
/*}}}*/
void TransientGriddedFieldParam::Marshall(MarshallHandle* marshallhandle){ /*{{{*/

	int object_enum = TransientGriddedFieldParamEnum;
   marshallhandle->call(object_enum);

	marshallhandle->call(this->enum_type);
	marshallhandle->call(this->interpolation);
	marshallhandle->call(this->cycle);
	marshallhandle->call(this->M);
	marshallhandle->call(this->N);
	marshallhandle->call(this->MN);
	marshallhandle->call(this->T);
	if(marshallhandle->OperationNumber()==MARSHALLING_LOAD){
		values    = xNew<IssmDouble>(MN*T);
		timesteps = xNew<IssmDouble>(T);
	}
	marshallhandle->call(this->values,MN*T);
	marshallhandle->call(this->timesteps,T);

}/*}}}*/
int  TransientGriddedFieldParam::ObjectEnum(void){/*{{{*/

	return TransientGriddedFieldParamEnum;

}/*}}}*/

/*TransientGriddedFieldParam virtual functions definitions: */
void  TransientGriddedFieldParam::GetParameterValue(IssmDouble* pdouble,int row,int column,IssmDouble time){/*{{{*/

	IssmDouble output;
	bool       found;
	_assert_(row>=0 && row<this->M); 
	_assert_(column>=0 && column<this->N); 

	if(this->cycle) _error_("not implemented yet");

	/*Ok, we have the time and row, go through the timesteps, and figure out which interval we 
	 *fall within. Then interpolate the values on this interval: */
	if(time<this->timesteps[0]){
		/*get values for the first time: */
		output=this->values[row*this->N+column];
		found=true;
	}
	else if(time>this->timesteps[this->T-1]){
		/*get values for the last time: */
		output=this->values[row*this->N+column+(this->T-1)*this->MN];
		found=true;
	}
	else{
		/*Find which interval we fall within: */
		for(int i=0;i<this->T;i++){
			if(time==this->timesteps[i]){
				/*We are right on one step time: */
				output = this->values[i*this->MN+row*this->N+column];
				found=true;
				break; //we are done with the time interpolation.
			}
			else{
				if(this->timesteps[i]<time && time<this->timesteps[i+1]){
					/*ok, we have the interval [i:i+1]. Interpolate linearly for now: */
					IssmDouble deltat = this->timesteps[i+1]-this->timesteps[i];
					IssmDouble alpha  = (time-this->timesteps[i])/deltat;
					if(this->interpolation==true) output=(1.0-alpha)*this->values[i*this->MN+row*this->N+column] + alpha*this->values[(i+1)*this->MN+row*this->N+column];
					else output=this->values[i*this->MN+row*this->N+column];
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
void  TransientGriddedFieldParam::GetParameterValue(IssmDouble* pdouble,int row,int column,IssmDouble time, int timestepping, IssmDouble dt){/*{{{*/

	bool cycle = this->cycle;

	_assert_(row>=0 && row<this->M); 
	_assert_(column>=0 && column<this->N); 

	if(cycle){

		/*Change input time if we cycle through the forcing*/
		IssmDouble time0 = this->timesteps[0];
		IssmDouble time1 = this->timesteps[this->T - 1];

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

		this->cycle=false;
	}

	this->GetParameterValue(pdouble,row,column,time);

	this->cycle = cycle;

}
/*}}}*/


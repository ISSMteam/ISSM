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
	this->MN=in_M*in_N;; //Number of data per time step
	this->T=in_T; //Number of time steps
	this->interpolation=interpolation_on;
	this->cycle=cycle_in;

	this->values=xNew<IssmDouble>(in_M*in_N*in_T);
	xMemCpy<IssmDouble>(this->values,in_values,in_M*in_N*in_T);

	this->timesteps=xNew<IssmDouble>(in_T);
	xMemCpy<IssmDouble>(this->timesteps,in_time,in_T);
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
	this->GetParameterValue(pdouble,NULL,row,column,time);
}
/*}}}*/
void  TransientGriddedFieldParam::GetParameterValue(IssmDouble* pdouble,int* index,int row,int column,IssmDouble time){/*{{{*/

	IssmDouble output;
	int output_id;
	bool       found;
	_assert_(row>=0 && row<this->M); 
	_assert_(column>=0 && column<this->N); 

	if(this->cycle) _error_("not implemented yet");

	/*Ok, we have the time and row, go through the timesteps, and figure out which interval we 
	 *fall within. Then interpolate the values on this interval: */
	if(time<this->timesteps[0]){
		/*get values for the first time: */
		output=this->values[(row*this->N+column)*this->T];
		output_id = -1;
		found=true;
	}
	else if(time>=this->timesteps[this->T-1]){
		/*get values for the last time: */
		output=this->values[(row*this->N+column)*this->T+(this->T-1)];
		output_id = this->T-1;
		found=true;
	}
	else{
		/*Find which interval we fall within: */
		for(int i=0;i<this->T-1;i++){
			if(time==this->timesteps[i]){
				/*We are right on one step time: */
				output = this->values[(row*this->N+column)*this->T+i];
				output_id = i;
				found=true;
				break; //we are done with the time interpolation.
			}
			else{
				if(this->timesteps[i]<time && time<this->timesteps[i+1]){
					/*ok, we have the interval [i:i+1]. Interpolate linearly for now: */
					IssmDouble deltat = this->timesteps[i+1]-this->timesteps[i];
					IssmDouble alpha  = (time-this->timesteps[i])/deltat;
					if(this->interpolation==true) output=(1.0-alpha)*this->values[(row*this->N+column)*this->T+i] + alpha*this->values[(row*this->N+column)*this->T+i+1];
					else output=this->values[(row*this->N+column)*this->T+i];
					found=true;
					output_id = i;
					break;
				}
				else continue; //keep looking on the next interval
			}
		}
	}
	if(!found)_error_("did not find time interval on which to interpolate values");
	if (pdouble != NULL) *pdouble=output;
	if (index !=NULL) *index=output_id;
}
/*}}}*/
void  TransientGriddedFieldParam::GetParameterValue(IssmDouble* pdouble,int row,int column,IssmDouble starttime,IssmDouble endtime){/*{{{*/
	/*compute average field between the given time period*/

	IssmDouble  output;
	IssmDouble  datastart, dataend;
	int         startid, endid;
	bool        found;
	IssmDouble yts=3600*24*365;
	_assert_(row>=0 && row<this->M);
	_assert_(column>=0 && column<this->N);
	_assert_(starttime<endtime);

	if(this->cycle) _error_("not implemented yet");

	/*Ok, we have the time and row, go through the timesteps, and figure out which interval we
	 *fall within. Then use trapzoidal rule to integrate over time and average: */
	if(endtime<=this->timesteps[0]){
		/*get values for the first time: */
		output=this->values[(row*this->N+column)*this->T];
		found=true;
	}
	else if(starttime>=this->timesteps[this->T-1]){
		/*get values for the last time: */
		output=this->values[(row*this->N+column)*this->T+(this->T-1)];
		found=true;
	}
	else{
		/*Find which interval we fall within: */
		this->GetParameterValue(&datastart,&startid,row,column,starttime);
		this->GetParameterValue(&dataend,&endid,row,column,endtime);
		/*include start, but the negative*/
		output = 0.5*(datastart+this->values[(row*this->N+column)*this->T+startid])*(this->timesteps[startid]-starttime);
		/*include end*/
		output += 0.5*(this->values[(row*this->N+column)*this->T+endid]+dataend)*(endtime-this->timesteps[endid]);
		/*integrate in betweem*/
		for (int i=startid;i<endid;i++) {
			output += 0.5*(this->values[(row*this->N+column)*this->T+i]+this->values[(row*this->N+column)*this->T+i+1])*(this->timesteps[i+1]-this->timesteps[i]);
		}
		output = output / (endtime-starttime);
		found = true;
	}

	if(!found)_error_("did not find time interval on which to interpolate values");
	*pdouble=output;
}
/*}}}*/
void  TransientGriddedFieldParam::GetParameterValue(IssmDouble** pIssmDoublearray,int* pM,int* pN,IssmDouble time){/*{{{*/

	IssmDouble* output = xNew<IssmDouble>(M*N);

	for (int i=0;i<N;i++) {
		for (int j=0;j<M;j++) {
			this->GetParameterValue(&output[j*N+i],j,i,time);
		}
	}
	/*Assign output pointers:*/
	if(pM) *pM=M;
	if(pN) *pN=N;
	*pIssmDoublearray=output;
}
/*}}}*/
void  TransientGriddedFieldParam::GetParameterValue(IssmDouble** pIssmDoublearray,int* pM,int* pN,IssmDouble starttime,IssmDouble endtime){/*{{{*/

	IssmDouble* output = xNew<IssmDouble>(M*N);

	for (int i=0;i<N;i++) {
		for (int j=0;j<M;j++) {
			this->GetParameterValue(&output[j*N+i],j,i,starttime,endtime);
		}
	}
	/*Assign output pointers:*/
	if(pM) *pM=M;
	if(pN) *pN=N;
	*pIssmDoublearray=output;
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


/*!\file StochasticForcingx
 * \brief: compute noise terms for the StochasticForcing fields
 */

#include "./StochasticForcingx.h"
#include "../../classes/Loads/Friction.h"
#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"
#include "../../shared/Random/random.h"

void StochasticForcingx(FemModel* femmodel){/*{{{*/

   /*Retrieve parameters*/
   bool randomflag;
   int M,N,numfields,numtcov,my_rank;
   int* fields                = NULL;
   int* dimensions            = NULL;
   IssmDouble* timecovariance = NULL;
   IssmDouble* covariance     = NULL;
   femmodel->parameters->FindParam(&randomflag,StochasticForcingRandomflagEnum);
   femmodel->parameters->FindParam(&numfields,StochasticForcingNumFieldsEnum);
   femmodel->parameters->FindParam(&numtcov,StochasticForcingNumTimesCovarianceEnum);
   femmodel->parameters->FindParam(&fields,&N,StochasticForcingFieldsEnum);    _assert_(N==numfields);
   femmodel->parameters->FindParam(&dimensions,&N,StochasticForcingDimensionsEnum);    _assert_(N==numfields);
   femmodel->parameters->FindParam(&timecovariance,&N,StochasticForcingTimeCovarianceEnum);    _assert_(N==numtcov);
   int dimtot=0;
   for(int i=0;i<numfields;i++) dimtot = dimtot+dimensions[i];
   femmodel->parameters->FindParam(&covariance,&M,&N,StochasticForcingCovarianceEnum); _assert_(M==numtcov); _assert_(N==dimtot*dimtot);

	/*Check if this is a timestep for new noiseterms computation*/
	bool isstepforstoch = false;
	IssmDouble time,dt,starttime,tstep_stoch;
   femmodel->parameters->FindParam(&time,TimeEnum);
   femmodel->parameters->FindParam(&dt,TimesteppingTimeStepEnum);
   femmodel->parameters->FindParam(&starttime,TimesteppingStartTimeEnum);
   femmodel->parameters->FindParam(&tstep_stoch,StochasticForcingTimestepEnum);

	/*Check if we use HydroarmaPw*/
	bool ispwHydro;
	femmodel->parameters->FindParam(&ispwHydro,HydrologyIsWaterPressureArmaEnum);

	#ifndef _HAVE_AD_
   if((fmod(time,tstep_stoch)<fmod((time-dt),tstep_stoch)) || (time<=starttime+dt) || tstep_stoch==dt) isstepforstoch = true;
   #else
   _error_("not implemented yet");
   #endif

   /*Compute noise terms*/
	IssmDouble* timestepcovariance = xNew<IssmDouble>(dimtot*dimtot);
	IssmDouble* noiseterms         = xNew<IssmDouble>(dimtot);
   if(isstepforstoch){
		/*Find covariance to be applied at current time step*/
		int itime;
		if(numtcov>1){
			for(int i=0;i<numtcov;i++){
				if(time>=timecovariance[i]) itime=i;
			}
		}
		else itime=0;
		for(int i=0;i<dimtot*dimtot;i++) timestepcovariance[i] = covariance[itime*dimtot*dimtot+i];
		my_rank=IssmComm::GetRank();
   	if(my_rank==0){
   	   int fixedseed;
			/*Determine whether random seed is fixed to time step (randomflag==false) or random seed truly random (randomflag==true)*/
   	   if(randomflag) fixedseed=-1;
   	   else fixedseed = reCast<int,IssmDouble>((time-starttime)/dt);
			/*multivariateNormal needs to be passed a NULL pointer to avoid memory leak issues*/
   	   IssmDouble* temparray = NULL;
   	   multivariateNormal(&temparray,dimtot,0.0,timestepcovariance,fixedseed);
   	   for(int i=0;i<dimtot;i++) noiseterms[i]=temparray[i];
			xDelete<IssmDouble>(temparray);
   	}
   	ISSM_MPI_Bcast(noiseterms,dimtot,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());
		femmodel->parameters->SetParam(noiseterms,dimtot,StochasticForcingNoisetermsEnum);
	}
	else{
		IssmDouble* temparray = NULL;
		femmodel->parameters->FindParam(&temparray,&N,StochasticForcingNoisetermsEnum); _assert_(N==dimtot);
		for(int i=0;i<dimtot;i++) noiseterms[i] = temparray[i];
		xDelete<IssmDouble>(temparray);
	}

	int i=0;
   for(int j=0;j<numfields;j++){
      int dimenum_type,noiseenum_type;
      IssmDouble* noisefield = xNew<IssmDouble>(dimensions[j]);
      for(int k=0;k<dimensions[j];k++){
         noisefield[k]=noiseterms[i+k];
      }

		int dimensionid;

		/*Deal with the ARMA models*/
		if(fields[j]==SMBarmaEnum || fields[j]==FrontalForcingsRignotarmaEnum || fields[j]==BasalforcingsDeepwaterMeltingRatearmaEnum || fields[j]==FrontalForcingsSubglacialDischargearmaEnum || (fields[j]==FrictionWaterPressureEnum && ispwHydro)){
			switch(fields[j]){
				case SMBarmaEnum:
					dimenum_type   = SmbBasinsIdEnum;
					noiseenum_type = SmbARMANoiseEnum;
					break;
				case FrontalForcingsRignotarmaEnum:
					dimenum_type   = FrontalForcingsBasinIdEnum;
					noiseenum_type = ThermalforcingARMANoiseEnum;
					break;
				case BasalforcingsDeepwaterMeltingRatearmaEnum:
					dimenum_type   = BasalforcingsLinearBasinIdEnum;
					noiseenum_type = BasalforcingsDeepwaterMeltingRateNoiseEnum;
					break;
				case FrontalForcingsSubglacialDischargearmaEnum:
					dimenum_type   = FrontalForcingsBasinIdEnum;
					noiseenum_type = SubglacialdischargeARMANoiseEnum;
					break;	
				case FrictionWaterPressureEnum:
					dimenum_type   = HydrologyBasinsIdEnum;
					noiseenum_type = FrictionWaterPressureNoiseEnum;
					break;	
			}
			for(Object* &object:femmodel->elements->objects){
            Element* element = xDynamicCast<Element*>(object);
            int numvertices  = element->GetNumberOfVertices();
            IssmDouble* noise_element = xNew<IssmDouble>(numvertices);
            element->GetInputValue(&dimensionid,dimenum_type);
            for(int i=0;i<numvertices;i++) noise_element[i] = noisefield[dimensionid];
            element->AddInput(noiseenum_type,noise_element,P0Enum);
            xDelete<IssmDouble>(noise_element);
			}
		}
		else{
			switch(fields[j]){
				case SMBarmaEnum:
				case FrontalForcingsRignotarmaEnum:
				case BasalforcingsDeepwaterMeltingRatearmaEnum:
				case FrontalForcingsSubglacialDischargearmaEnum:
					/*Already done above*/
					break;
				case BasalforcingsSpatialDeepwaterMeltingRateEnum:
               /*Delete BasalforcingsSpatialDeepwaterMeltingRateEnum at previous time step (required if it is transient)*/
               femmodel->inputs->DeleteInput(BasalforcingsSpatialDeepwaterMeltingRateEnum);
               for(Object* &object:femmodel->elements->objects){
                  Element* element = xDynamicCast<Element*>(object);
                  int numvertices  = element->GetNumberOfVertices();
                  IssmDouble baselinedeepwatermelt;
                  IssmDouble deepwatermelt_tot[10]; _assert_(numvertices<10);
                  Input* baselinedeepwatermelt_input  = NULL;
                  baselinedeepwatermelt_input = element->GetInput(BaselineBasalforcingsSpatialDeepwaterMeltingRateEnum); _assert_(baselinedeepwatermelt_input);
                  element->GetInputValue(&dimensionid,StochasticForcingDefaultIdEnum);
                  Gauss* gauss = element->NewGauss();
                  for(int i=0;i<numvertices;i++){
                     gauss->GaussVertex(i);
                     baselinedeepwatermelt_input->GetInputValue(&baselinedeepwatermelt,gauss);
                     deepwatermelt_tot[i] = baselinedeepwatermelt+noisefield[dimensionid];
                  }
                  element->AddInput(BasalforcingsSpatialDeepwaterMeltingRateEnum,&deepwatermelt_tot[0],P1DGEnum);
                  delete gauss;
               }
               break;
				case DefaultCalvingEnum:
					/*Delete CalvingCalvingrateEnum at previous time step (required if it is transient)*/
					femmodel->inputs->DeleteInput(CalvingCalvingrateEnum);
					for(Object* &object:femmodel->elements->objects){
						Element* element = xDynamicCast<Element*>(object);
						int numvertices  = element->GetNumberOfVertices();
						IssmDouble baselinecalvingrate;
						IssmDouble calvingrate_tot[10]; _assert_(numvertices<10);
						Input* baselinecalvingrate_input  = NULL;
						baselinecalvingrate_input = element->GetInput(BaselineCalvingCalvingrateEnum); _assert_(baselinecalvingrate_input);
						element->GetInputValue(&dimensionid,StochasticForcingDefaultIdEnum);
						Gauss* gauss = element->NewGauss();
						for(int i=0;i<numvertices;i++){
							gauss->GaussVertex(i);
							baselinecalvingrate_input->GetInputValue(&baselinecalvingrate,gauss);
							calvingrate_tot[i] = max(0.0,baselinecalvingrate+noisefield[dimensionid]);
						}
						element->AddInput(CalvingCalvingrateEnum,&calvingrate_tot[0],P1DGEnum);
						delete gauss;
					}
					break;
				case FloatingMeltRateEnum:
					/*Delete BasalforcingsFloatingiceMeltingRateEnum at previous time step (required if it is transient)*/
					femmodel->inputs->DeleteInput(BasalforcingsFloatingiceMeltingRateEnum);
					for(Object* &object:femmodel->elements->objects){
						Element* element = xDynamicCast<Element*>(object);
						int numvertices  = element->GetNumberOfVertices();
						IssmDouble baselinefloatingicemeltrate;
						IssmDouble floatingicemeltrate_tot[10]; _assert_(numvertices<10);
						Input* baselinefloatingicemeltrate_input  = NULL;
						baselinefloatingicemeltrate_input = element->GetInput(BaselineBasalforcingsFloatingiceMeltingRateEnum); _assert_(baselinefloatingicemeltrate_input);
						element->GetInputValue(&dimensionid,StochasticForcingDefaultIdEnum);
						Gauss* gauss = element->NewGauss();
						for(int i=0;i<numvertices;i++){
							gauss->GaussVertex(i);
							baselinefloatingicemeltrate_input->GetInputValue(&baselinefloatingicemeltrate,gauss);
							/*No check for positive melt rate because basal accretion is allowed*/
							floatingicemeltrate_tot[i] = baselinefloatingicemeltrate+noisefield[dimensionid];
						}
						element->AddInput(BasalforcingsFloatingiceMeltingRateEnum,&floatingicemeltrate_tot[0],P1DGEnum);
						delete gauss;
					}
					break;
				case SMBforcingEnum:
					/*Delete SmbMassBalanceEnum at previous time step (required if it is transient)*/
					femmodel->inputs->DeleteInput(SmbMassBalanceEnum);
					for(Object* &object:femmodel->elements->objects){
						Element* element = xDynamicCast<Element*>(object);
						int numvertices  = element->GetNumberOfVertices();
						IssmDouble baselinesmb;
						IssmDouble smb_tot[10]; _assert_(numvertices<10);
						Input* baselinesmb_input  = NULL;
						baselinesmb_input = element->GetInput(BaselineSmbMassBalanceEnum); _assert_(baselinesmb_input);
						element->GetInputValue(&dimensionid,StochasticForcingDefaultIdEnum);
						Gauss* gauss = element->NewGauss();
						for(int i=0;i<numvertices;i++){
							gauss->GaussVertex(i);
							baselinesmb_input->GetInputValue(&baselinesmb,gauss);
							smb_tot[i] = baselinesmb+noisefield[dimensionid];
						}
						element->AddInput(SmbMassBalanceEnum,&smb_tot[0],P1DGEnum);
						delete gauss;
					}
					break;
				case FrictionWaterPressureEnum:
					/*Specify that WaterPressure is stochastic*/ 
					femmodel->parameters->SetParam(true,StochasticForcingIsWaterPressureEnum);
					for(Object* &object:femmodel->elements->objects){
                  Element* element = xDynamicCast<Element*>(object);
                  int numvertices  = element->GetNumberOfVertices();
                  IssmDouble p_water_deterministic[10]; _assert_(numvertices<10);
                  IssmDouble p_water[10];
						element->GetInputValue(&dimensionid,StochasticForcingDefaultIdEnum);
						element->SubglacialWaterPressure(FrictionWaterPressureEnum);
                  element->GetInputListOnVertices(&p_water_deterministic[0],FrictionWaterPressureEnum);
                  for(int i=0;i<numvertices;i++) p_water[i] = p_water_deterministic[i] + noisefield[dimensionid];
                  element->AddInput(FrictionWaterPressureEnum,p_water,P1DGEnum);
					}
					break;
				default:
					_error_("Field "<<EnumToStringx(fields[j])<<" does not support stochasticity yet.");
			}
		}
		i=i+dimensions[j];
      xDelete<IssmDouble>(noisefield);
   }

	/*Cleanup*/
   xDelete<int>(fields);
   xDelete<int>(dimensions);
   xDelete<IssmDouble>(covariance);
   xDelete<IssmDouble>(timecovariance);
   xDelete<IssmDouble>(timestepcovariance);
   xDelete<IssmDouble>(noiseterms);
}/*}}}*/

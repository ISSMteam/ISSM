/*!\file FrontalForcingsx
 * \brief: compute ice frontal melting rate
 */

#include "./FrontalForcingsx.h"
#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"
#include "../../shared/Random/random.h"

void FrontalForcingsx(FemModel* femmodel){/*{{{*/

	/*Recover melt_parameterization*/
	int melt_parameterization;
	femmodel->parameters->FindParam(&melt_parameterization,FrontalForcingsParamEnum);

	/*Calculate melting rate*/
	switch(melt_parameterization){
		case FrontalForcingsDefaultEnum:
			break;
		case FrontalForcingsRignotarmaEnum:
			Thermalforcingarmax(femmodel);
			bool isdischargearma;
			femmodel->parameters->FindParam(&isdischargearma,FrontalForcingsIsDischargeARMAEnum);
			if(isdischargearma==true) Subglacialdischargearmax(femmodel);
			/*Do not break here, call IcefrontAreax(),RignotMeltParameterizationx()*/
		case FrontalForcingsRignotEnum:
			femmodel->IcefrontAreax();
			femmodel->RignotMeltParameterizationx();
			break;
		default:
			_error_("Frontal forcings "<<EnumToStringx(melt_parameterization)<<" not supported yet");
	}
}/*}}}*/
void Thermalforcingarmax(FemModel* femmodel){/*{{{*/

   /*Get time parameters*/
   IssmDouble time,dt,starttime,tstep_arma;
   femmodel->parameters->FindParam(&time,TimeEnum);
   femmodel->parameters->FindParam(&dt,TimesteppingTimeStepEnum);
   femmodel->parameters->FindParam(&starttime,TimesteppingStartTimeEnum);
   femmodel->parameters->FindParam(&tstep_arma,FrontalForcingsARMATimestepEnum);

   /*Determine if this is a time step for the ARMA model*/
   bool isstepforarma = false;

   #ifndef _HAVE_AD_
   if((fmod(time,tstep_arma)<fmod((time-dt),tstep_arma)) || (time<=starttime+dt) || tstep_arma==dt) isstepforarma = true;
   #else
   _error_("not implemented yet");
   #endif

   /*Load parameters*/
	bool isstochastic;
   bool istfstochastic = false;
	int M,N,arorder,maorder,numbasins,numparams,numbreaks,nummonthbreaks,my_rank;
   femmodel->parameters->FindParam(&numbasins,FrontalForcingsNumberofBasinsEnum);
   femmodel->parameters->FindParam(&numparams,FrontalForcingsNumberofParamsEnum);
   femmodel->parameters->FindParam(&numbreaks,FrontalForcingsNumberofBreaksEnum);
   femmodel->parameters->FindParam(&nummonthbreaks,FrontalForcingsNumberofMonthBreaksEnum);
   femmodel->parameters->FindParam(&arorder,FrontalForcingsARMAarOrderEnum);
   femmodel->parameters->FindParam(&maorder,FrontalForcingsARMAmaOrderEnum);
   IssmDouble* datebreaks        = NULL;
   IssmDouble* arlagcoefs        = NULL;
   IssmDouble* malagcoefs        = NULL;
   IssmDouble* monthlyeff        = NULL;
	IssmDouble* polyparams        = NULL;
	IssmDouble* monthdatebreaks   = NULL;
	IssmDouble* monthintercepts   = NULL;
	IssmDouble* monthtrends       = NULL;

   femmodel->parameters->FindParam(&datebreaks,&M,&N,FrontalForcingsARMAdatebreaksEnum);            _assert_(M==numbasins); _assert_(N==max(numbreaks,1));        
   femmodel->parameters->FindParam(&polyparams,&M,&N,FrontalForcingsARMApolyparamsEnum);            _assert_(M==numbasins); _assert_(N==(numbreaks+1)*numparams);        
   femmodel->parameters->FindParam(&arlagcoefs,&M,&N,FrontalForcingsARMAarlagcoefsEnum);            _assert_(M==numbasins); _assert_(N==arorder);
   femmodel->parameters->FindParam(&malagcoefs,&M,&N,FrontalForcingsARMAmalagcoefsEnum);            _assert_(M==numbasins); _assert_(N==maorder);
   femmodel->parameters->FindParam(&monthdatebreaks,&M,&N,FrontalForcingsARMAmonthdatebreaksEnum);  _assert_(M==numbasins); _assert_(N==max(nummonthbreaks,1));        
   femmodel->parameters->FindParam(&monthintercepts,&M,&N,FrontalForcingsARMAmonthinterceptsEnum);  _assert_(M==numbasins); _assert_(N==12*(nummonthbreaks+1)); 
   femmodel->parameters->FindParam(&monthtrends,&M,&N,FrontalForcingsARMAmonthtrendsEnum);          _assert_(M==numbasins); _assert_(N==12*(nummonthbreaks+1)); 

	femmodel->parameters->FindParam(&isstochastic,StochasticForcingIsStochasticForcingEnum);
	if(isstochastic){
		int  numstochasticfields;
		int* stochasticfields;
		femmodel->parameters->FindParam(&numstochasticfields,StochasticForcingNumFieldsEnum);
		femmodel->parameters->FindParam(&stochasticfields,&N,StochasticForcingFieldsEnum); _assert_(N==numstochasticfields);
		for(int i=0;i<numstochasticfields;i++){
			if(stochasticfields[i]==FrontalForcingsRignotarmaEnum) istfstochastic = true;
		}
		xDelete<int>(stochasticfields);
	}

   /*Loop over each element to compute Thermal Forcing at vertices*/
   for(Object* &object:femmodel->elements->objects){
      Element* element = xDynamicCast<Element*>(object);
		/*Compute ARMA*/
      element->ArmaProcess(isstepforarma,arorder,maorder,numparams,numbreaks,tstep_arma,polyparams,arlagcoefs,malagcoefs,datebreaks,istfstochastic,FrontalForcingsRignotarmaEnum);
		/*Compute monthly effects*/
		element->MonthlyPiecewiseLinearEffectBasin(nummonthbreaks,monthintercepts,monthtrends,monthdatebreaks,FrontalForcingsRignotarmaEnum);
	}

   /*Cleanup*/
   xDelete<IssmDouble>(arlagcoefs);
   xDelete<IssmDouble>(malagcoefs);
   xDelete<IssmDouble>(monthlyeff);
   xDelete<IssmDouble>(polyparams);
   xDelete<IssmDouble>(datebreaks);
   xDelete<IssmDouble>(monthdatebreaks);
   xDelete<IssmDouble>(monthintercepts);
   xDelete<IssmDouble>(monthtrends);
}/*}}}*/
void Subglacialdischargearmax(FemModel* femmodel){/*{{{*/

	/*Get time parameters*/
   IssmDouble time,dt,starttime,tstep_arma;
   femmodel->parameters->FindParam(&time,TimeEnum);
   femmodel->parameters->FindParam(&dt,TimesteppingTimeStepEnum);
   femmodel->parameters->FindParam(&starttime,TimesteppingStartTimeEnum);
   femmodel->parameters->FindParam(&tstep_arma,FrontalForcingsSdARMATimestepEnum);

   /*Determine if this is a time step for the ARMA model*/
   bool isstepforarma = false;

   #ifndef _HAVE_AD_
   if((fmod(time,tstep_arma)<fmod((time-dt),tstep_arma)) || (time<=starttime+dt) || tstep_arma==dt) isstepforarma = true;
   #else
   _error_("not implemented yet");
   #endif

   /*Load parameters*/
	bool isstochastic;
   bool isdischargestochastic = false;
	int M,N,arorder,maorder,numbasins,numparams,numbreaks,my_rank;
   femmodel->parameters->FindParam(&numbasins,FrontalForcingsNumberofBasinsEnum);
   femmodel->parameters->FindParam(&numparams,FrontalForcingsSdNumberofParamsEnum);
   femmodel->parameters->FindParam(&numbreaks,FrontalForcingsSdNumberofBreaksEnum);
   femmodel->parameters->FindParam(&arorder,FrontalForcingsSdarOrderEnum);
   femmodel->parameters->FindParam(&maorder,FrontalForcingsSdmaOrderEnum);
   IssmDouble* datebreaks        = NULL;
   IssmDouble* arlagcoefs        = NULL;
   IssmDouble* malagcoefs        = NULL;
   IssmDouble* monthlyfrac       = NULL;
	IssmDouble* polyparams        = NULL;

   femmodel->parameters->FindParam(&datebreaks,&M,&N,FrontalForcingsSddatebreaksEnum);            _assert_(M==numbasins); _assert_(N==max(numbreaks,1));        
   femmodel->parameters->FindParam(&polyparams,&M,&N,FrontalForcingsSdpolyparamsEnum);            _assert_(M==numbasins); _assert_(N==(numbreaks+1)*numparams);        
   femmodel->parameters->FindParam(&arlagcoefs,&M,&N,FrontalForcingsSdarlagcoefsEnum);            _assert_(M==numbasins); _assert_(N==arorder);
   femmodel->parameters->FindParam(&malagcoefs,&M,&N,FrontalForcingsSdmalagcoefsEnum);            _assert_(M==numbasins); _assert_(N==maorder);
   femmodel->parameters->FindParam(&monthlyfrac,&M,&N,FrontalForcingsSdMonthlyFracEnum);          _assert_(M==numbasins); _assert_(N==12); 

	femmodel->parameters->FindParam(&isstochastic,StochasticForcingIsStochasticForcingEnum);
	if(isstochastic){
		int  numstochasticfields;
		int* stochasticfields;
		femmodel->parameters->FindParam(&numstochasticfields,StochasticForcingNumFieldsEnum);
		femmodel->parameters->FindParam(&stochasticfields,&N,StochasticForcingFieldsEnum); _assert_(N==numstochasticfields);
		for(int i=0;i<numstochasticfields;i++){
			if(stochasticfields[i]==FrontalForcingsSubglacialDischargearmaEnum) isdischargestochastic = true;
		}
		xDelete<int>(stochasticfields);
	}

   /*Loop over each element to compute Subglacial Discharge at vertices*/
   for(Object* &object:femmodel->elements->objects){
      Element* element = xDynamicCast<Element*>(object);
		/*Compute ARMA*/
      element->ArmaProcess(isstepforarma,arorder,maorder,numparams,numbreaks,tstep_arma,polyparams,arlagcoefs,malagcoefs,datebreaks,isdischargestochastic,FrontalForcingsSubglacialDischargearmaEnum);
		/*Scale with monthly fractions*/
		element->MonthlyFactorBasin(monthlyfrac,FrontalForcingsSubglacialDischargearmaEnum);
	}

   /*Cleanup*/
   xDelete<IssmDouble>(arlagcoefs);
   xDelete<IssmDouble>(malagcoefs);
   xDelete<IssmDouble>(monthlyfrac);
   xDelete<IssmDouble>(polyparams);
   xDelete<IssmDouble>(datebreaks);
}/*}}}*/


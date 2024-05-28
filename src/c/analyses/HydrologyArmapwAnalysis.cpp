#include "./HydrologyArmapwAnalysis.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"

/*Model processing*/
void HydrologyArmapwAnalysis::CreateConstraints(Constraints* constraints,IoModel* iomodel){/*{{{*/

	return;

}/*}}}*/
void HydrologyArmapwAnalysis::CreateLoads(Loads* loads, IoModel* iomodel){/*{{{*/
	/*No loads*/
}/*}}}*/
void HydrologyArmapwAnalysis::CreateNodes(Nodes* nodes,IoModel* iomodel,bool isamr){/*{{{*/

	return;

}/*}}}*/
int  HydrologyArmapwAnalysis::DofsPerNode(int** doflist,int domaintype,int approximation){/*{{{*/
	return 0;
}/*}}}*/
void HydrologyArmapwAnalysis::UpdateElements(Elements* elements,Inputs* inputs,IoModel* iomodel,int analysis_counter,int analysis_type){/*{{{*/

	/*Fetch data needed: */
   int    hydrology_model,frictionlaw;
   iomodel->FindConstant(&hydrology_model,"md.hydrology.model");

   /*Now, do we really want armapw?*/
   if(hydrology_model!=HydrologyarmapwEnum) return;

   /*Add input to elements*/
   iomodel->FetchDataToInput(inputs,elements,"md.mask.ice_levelset",MaskIceLevelsetEnum);
   iomodel->FetchDataToInput(inputs,elements,"md.mask.ocean_levelset",MaskOceanLevelsetEnum);
   iomodel->FetchDataToInput(inputs,elements,"md.basalforcings.groundedice_melting_rate",BasalforcingsGroundediceMeltingRateEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.hydrology.basin_id",HydrologyBasinsIdEnum);
   iomodel->FetchDataToInput(inputs,elements,"md.initialization.watercolumn",WatercolumnEnum,0.);

}/*}}}*/
void HydrologyArmapwAnalysis::UpdateParameters(Parameters* parameters,IoModel* iomodel,int solution_enum,int analysis_enum){/*{{{*/

	/*retrieve some parameters: */
	int    hydrology_model;
	int    numoutputs;
	char** requestedoutputs = NULL;
	iomodel->FindConstant(&hydrology_model,"md.hydrology.model");

	/*Now, do we really want Armapw?*/
	if(hydrology_model!=HydrologyarmapwEnum) return;

	parameters->AddObject(new IntParam(HydrologyModelEnum,hydrology_model));
  /*Requested outputs*/
  iomodel->FindConstant(&requestedoutputs,&numoutputs,"md.hydrology.requested_outputs");
  parameters->AddObject(new IntParam(HydrologyNumRequestedOutputsEnum,numoutputs));
  if(numoutputs)parameters->AddObject(new StringArrayParam(HydrologyRequestedOutputsEnum,requestedoutputs,numoutputs));
  iomodel->DeleteData(&requestedoutputs,numoutputs,"md.hydrology.requested_outputs");

}/*}}}*/

/*Finite Element Analysis*/
void           HydrologyArmapwAnalysis::Core(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
void           HydrologyArmapwAnalysis::PreCore(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
ElementVector* HydrologyArmapwAnalysis::CreateDVector(Element* element){/*{{{*/
	/*Default, return NULL*/
	return NULL;
}/*}}}*/
ElementMatrix* HydrologyArmapwAnalysis::CreateJacobianMatrix(Element* element){/*{{{*/
	_error_("Not implemented");
}/*}}}*/
ElementMatrix* HydrologyArmapwAnalysis::CreateKMatrix(Element* element){/*{{{*/
	_error_("Not implemented");
}/*}}}*/
ElementVector* HydrologyArmapwAnalysis::CreatePVector(Element* element){/*{{{*/
	_error_("Not implemented");
}/*}}}*/
void           HydrologyArmapwAnalysis::GetSolutionFromInputs(Vector<IssmDouble>* solution,Element* element){/*{{{*/
	_error_("Not implemented");
}/*}}}*/
void           HydrologyArmapwAnalysis::GradientJ(Vector<IssmDouble>* gradient,Element*  element,int control_type,int control_interp,int control_index){/*{{{*/
	_error_("Not implemented");
}/*}}}*/
void           HydrologyArmapwAnalysis::InputUpdateFromSolution(IssmDouble* solution,Element* element){/*{{{*/
	_error_("Not implemented");
}/*}}}*/
void           HydrologyArmapwAnalysis::UpdateConstraints(FemModel* femmodel){/*{{{*/
	_error_("Not implemented");
}/*}}}*/

/*Additional methods*/
void HydrologyArmapwAnalysis::UpdateSubglacialWaterPressure(FemModel* femmodel){/*{{{*/

	/*Get time parameters*/
   IssmDouble time,dt,starttime,tstep_arma;
   femmodel->parameters->FindParam(&time,TimeEnum);
   femmodel->parameters->FindParam(&dt,TimesteppingTimeStepEnum);
   femmodel->parameters->FindParam(&starttime,TimesteppingStartTimeEnum);
   femmodel->parameters->FindParam(&tstep_arma,HydrologyarmaTimestepEnum);

	/*Determine if this is a time step for the ARMA model*/
   bool isstepforarma = false;

   #ifndef _HAVE_AD_
   if((fmod(time,tstep_arma)<fmod((time-dt),tstep_arma)) || (time<=starttime+dt) || tstep_arma==dt) isstepforarma = true;
   #else
   _error_("not implemented yet");
   #endif

   /*Load parameters*/
   bool isstochastic;
   bool ispwstochastic = false;
   int M,N,arorder,maorder,numbasins,numparams,numbreaks,my_rank;
   femmodel->parameters->FindParam(&numbasins,HydrologyNumBasinsEnum);
   femmodel->parameters->FindParam(&numparams,HydrologyarmaNumParamsEnum);
   femmodel->parameters->FindParam(&numbreaks,HydrologyarmaNumBreaksEnum);
   femmodel->parameters->FindParam(&arorder,HydrologyarmaarOrderEnum);
   femmodel->parameters->FindParam(&maorder,HydrologyarmamaOrderEnum);
   IssmDouble* datebreaks        = NULL;
   IssmDouble* arlagcoefs        = NULL;
   IssmDouble* malagcoefs        = NULL;
	IssmDouble* monthlyfactors    = NULL;
   IssmDouble* polyparams        = NULL;	

	femmodel->parameters->FindParam(&datebreaks,&M,&N,HydrologyarmadatebreaksEnum);            _assert_(M==numbasins); _assert_(N==max(numbreaks,1));
   femmodel->parameters->FindParam(&polyparams,&M,&N,HydrologyarmapolyparamsEnum);            _assert_(M==numbasins); _assert_(N==(numbreaks+1)*numparams);
   femmodel->parameters->FindParam(&arlagcoefs,&M,&N,HydrologyarmaarlagcoefsEnum);            _assert_(M==numbasins); _assert_(N==arorder);
   femmodel->parameters->FindParam(&malagcoefs,&M,&N,HydrologyarmamalagcoefsEnum);            _assert_(M==numbasins); _assert_(N==maorder);
	femmodel->parameters->FindParam(&monthlyfactors,&M,&N,HydrologyarmaMonthlyFactorsEnum);    _assert_(M==numbasins); _assert_(N==12);

	femmodel->parameters->FindParam(&isstochastic,StochasticForcingIsStochasticForcingEnum);
   if(isstochastic){
      int  numstochasticfields;
      int* stochasticfields;
      femmodel->parameters->FindParam(&numstochasticfields,StochasticForcingNumFieldsEnum);
      femmodel->parameters->FindParam(&stochasticfields,&N,StochasticForcingFieldsEnum); _assert_(N==numstochasticfields);
      for(int i=0;i<numstochasticfields;i++){
         if(stochasticfields[i]==FrictionWaterPressureEnum) ispwstochastic = true;
      }
      xDelete<int>(stochasticfields);
   }

	/*Check if seasonality is imposed*/
	bool isseasonality = false;
	for(int i=0;i<numbasins*12;i++){
		if(monthlyfactors[i]!=1) isseasonality = true;
	}

	/*Loop over each element to compute Subglacial Water Pressure at vertices*/
   for(Object* &object:femmodel->elements->objects){
      Element* element = xDynamicCast<Element*>(object);
      /*Compute ARMA perturbation values*/
      element->ArmaProcess(isstepforarma,arorder,maorder,numparams,numbreaks,tstep_arma,polyparams,arlagcoefs,malagcoefs,datebreaks,ispwstochastic,HydrologyarmapwEnum);
      /*Compute subglacial water pressure with the ARMA perturbation*/
		element->SubglacialWaterPressure(FrictionWaterPressureEnum);
		/*Scale with monthly factors*/
		if(isseasonality) element->MonthlyFactorBasin(monthlyfactors,HydrologyarmapwEnum);
   }

	/*Cleanup*/
   xDelete<IssmDouble>(arlagcoefs);
   xDelete<IssmDouble>(malagcoefs);
   xDelete<IssmDouble>(polyparams);
   xDelete<IssmDouble>(datebreaks);
   xDelete<IssmDouble>(monthlyfactors);
}/*}}}*/


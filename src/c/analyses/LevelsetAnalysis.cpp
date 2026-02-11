#ifdef HAVE_CONFIG_H
#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif
#include "./LevelsetAnalysis.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"
#include "../solutionsequences/solutionsequences.h"
#include <math.h>

void LevelsetAnalysis::CreateConstraints(Constraints* constraints,IoModel* iomodel){/*{{{*/

	/*intermediary: */
	int finiteelement;
	int         code,vector_layout;
	IssmDouble *spcdata = NULL;
	int         M,N;

	/*Get finite element type for this analysis*/
	iomodel->FindConstant(&finiteelement,"md.levelset.fe");

	/*First of, find the record for the enum, and get code  of data type: */
	iomodel->SetFilePointerToData(&code, &vector_layout,"md.levelset.spclevelset");
	if(code!=7)_error_("expecting a IssmDouble vector for constraints md.levelset.spclevelset");
	if(vector_layout!=1)_error_("expecting a nodal vector for constraints md.levelset.spclevelset");

	/*Fetch vector:*/
	iomodel->FetchData(&spcdata,&M,&N,"md.levelset.spclevelset");

	/*Call IoModelToConstraintsx*/
	if(N>1){
		/*If it is a time series, most likely we are forcing the ice front position and do not want to have a Dynamic Constraint*/
		_assert_(M==iomodel->numberofvertices+1);
		IoModelToConstraintsx(constraints,iomodel,spcdata,M,N,LevelsetAnalysisEnum,finiteelement);
	}
	else{
		/*This is not a time series, we probably have calving on, we need the levelset constraints to update as the levelset moves*/
		_assert_(N==1);
		_assert_(M==iomodel->numberofvertices);
		IoModelToDynamicConstraintsx(constraints,iomodel,spcdata,M,N,LevelsetAnalysisEnum,finiteelement);
	}

	/*Clean up*/
	xDelete<IssmDouble>(spcdata);
}
/*}}}*/
void LevelsetAnalysis::CreateLoads(Loads* loads, IoModel* iomodel){/*{{{*/
	return;
}/*}}}*/
void LevelsetAnalysis::CreateNodes(Nodes* nodes,IoModel* iomodel,bool isamr){/*{{{*/
	int finiteelement;
	iomodel->FindConstant(&finiteelement,"md.levelset.fe");
	if(iomodel->domaintype!=Domain2DhorizontalEnum) iomodel->FetchData(2,"md.mesh.vertexonbase","md.mesh.vertexonsurface");
	::CreateNodes(nodes,iomodel,LevelsetAnalysisEnum,finiteelement);
	iomodel->DeleteData(2,"md.mesh.vertexonbase","md.mesh.vertexonsurface");
}
/*}}}*/
int  LevelsetAnalysis::DofsPerNode(int** doflist,int domaintype,int approximation){/*{{{*/
	return 1;
}
/*}}}*/
void LevelsetAnalysis::UpdateElements(Elements* elements,Inputs* inputs,IoModel* iomodel,int analysis_counter,int analysis_type){/*{{{*/

	/*Finite element type*/
	int finiteelement;
	iomodel->FindConstant(&finiteelement,"md.levelset.fe");

	/*Update elements: */
	int counter=0;
	for(int i=0;i<iomodel->numberofelements;i++){
		if(iomodel->my_elements[i]){
			Element* element=(Element*)elements->GetObjectByOffset(counter);
			element->Update(inputs,i,iomodel,analysis_counter,analysis_type,finiteelement);
			counter++;
		}
	}

	iomodel->FetchDataToInput(inputs,elements,"md.mask.ice_levelset",MaskIceLevelsetEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.mask.ocean_levelset",MaskOceanLevelsetEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.initialization.vx",VxEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.initialization.vy",VyEnum);

	/*Get moving front parameters*/
	bool isstochastic;
   int  calvinglaw;
   iomodel->FindConstant(&calvinglaw,"md.calving.law");
   iomodel->FindConstant(&isstochastic,"md.stochasticforcing.isstochasticforcing");
   switch(calvinglaw){

		/*"Continuous" calving laws*/
      case DefaultCalvingEnum:
         iomodel->FetchDataToInput(inputs,elements,"md.calving.calvingrate",CalvingCalvingrateEnum);
         if(isstochastic){
				iomodel->FetchDataToInput(inputs,elements,"md.stochasticforcing.default_id",StochasticForcingDefaultIdEnum);
            iomodel->FetchDataToInput(inputs,elements,"md.calving.calvingrate",BaselineCalvingCalvingrateEnum);
         }
         break;	
		case CalvingLevermannEnum:
			iomodel->FetchDataToInput(inputs,elements,"md.calving.coeff",CalvinglevermannCoeffEnum);
			break;
		case CalvingVonmisesEnum:
			iomodel->FetchDataToInput(inputs,elements,"md.calving.stress_threshold_groundedice",CalvingStressThresholdGroundediceEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.calving.stress_threshold_floatingice",CalvingStressThresholdFloatingiceEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.geometry.bed",BedEnum);
			break;
		case CalvingVonmisesADEnum:
			iomodel->FetchDataToInput(inputs,elements,"md.calving.basin_id",CalvingBasinIdEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.geometry.bed",BedEnum);
			break;
		case CalvingDev2Enum:
			iomodel->FetchDataToInput(inputs,elements,"md.calving.stress_threshold_groundedice",CalvingStressThresholdGroundediceEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.calving.stress_threshold_floatingice",CalvingStressThresholdFloatingiceEnum);
			break;
		case CalvingTestEnum:
			break;
		case CalvingParameterizationEnum:
			iomodel->FetchDataToInput(inputs,elements,"md.geometry.bed",BedEnum);
			break;
		case CalvingCalvingMIPEnum:
			break;

		/*"Discrete" calving laws (need to specify rate as 0 so that we can still solve the level set equation)*/
		case CalvingMinthicknessEnum:
			iomodel->FetchDataToInput(inputs,elements,"md.geometry.bed",BedEnum);
			iomodel->ConstantToInput(inputs,elements,0.,CalvingratexEnum,P1Enum);
			iomodel->ConstantToInput(inputs,elements,0.,CalvingrateyEnum,P1Enum);
			break;
		case CalvingHabEnum:
			iomodel->FetchDataToInput(inputs,elements,"md.calving.flotation_fraction",CalvingHabFractionEnum);
			iomodel->ConstantToInput(inputs,elements,0.,CalvingratexEnum,P1Enum);
			iomodel->ConstantToInput(inputs,elements,0.,CalvingrateyEnum,P1Enum);
			break;
		case CalvingCrevasseDepthEnum:
			iomodel->FetchDataToInput(inputs,elements,"md.calving.water_height",WaterheightEnum);
			iomodel->ConstantToInput(inputs,elements,0.,CalvingratexEnum,P1Enum);
			iomodel->ConstantToInput(inputs,elements,0.,CalvingrateyEnum,P1Enum);
			break;
		case CalvingPollardEnum:
			break;

		default:
			_error_("Calving law "<<EnumToStringx(calvinglaw)<<" not supported yet");
	}

	/*Get frontal melt parameters*/
	int melt_parameterization;
	iomodel->FindConstant(&melt_parameterization,"md.frontalforcings.parameterization");
	switch(melt_parameterization){
		case FrontalForcingsDefaultEnum:
			iomodel->FetchDataToInput(inputs,elements,"md.frontalforcings.meltingrate",CalvingMeltingrateEnum);
			if ((calvinglaw == CalvingParameterizationEnum) || (calvinglaw == CalvingCalvingMIPEnum)) {
				iomodel->FetchDataToInput(inputs,elements,"md.frontalforcings.ablationrate",CalvingAblationrateEnum);
			}
			break;
		case FrontalForcingsRignotEnum:
         /*Retrieve thermal forcing only in the case of non-arma FrontalForcingsRignot*/
         iomodel->FetchDataToInput(inputs,elements,"md.frontalforcings.thermalforcing",ThermalForcingEnum);
         iomodel->FetchDataToInput(inputs,elements,"md.frontalforcings.basin_id",FrontalForcingsBasinIdEnum);
         iomodel->FetchDataToInput(inputs,elements,"md.frontalforcings.subglacial_discharge",FrontalForcingsSubglacialDischargeEnum);
			break;
		case FrontalForcingsRignotarmaEnum:
			bool isdischargearma;
			iomodel->FindConstant(&isdischargearma,"md.frontalforcings.isdischargearma");
         iomodel->FetchDataToInput(inputs,elements,"md.frontalforcings.basin_id",FrontalForcingsBasinIdEnum);
         if(isdischargearma==false) iomodel->FetchDataToInput(inputs,elements,"md.frontalforcings.subglacial_discharge",FrontalForcingsSubglacialDischargeEnum);
			break;	
		default:
			_error_("Frontal forcings"<<EnumToStringx(melt_parameterization)<<" not supported yet");
	}
}
/*}}}*/
void LevelsetAnalysis::UpdateParameters(Parameters* parameters,IoModel* iomodel,int solution_enum,int analysis_enum){/*{{{*/

	parameters->AddObject(iomodel->CopyConstantObject("md.levelset.stabilization",LevelsetStabilizationEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.levelset.reinit_frequency",LevelsetReinitFrequencyEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.levelset.kill_icebergs",LevelsetKillIcebergsEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.levelset.migration_max",MigrationMaxEnum));

	int  calvinglaw;
   IssmDouble *transparam = NULL;
   IssmDouble  yts;
   int         N,M;
   bool        interp,cycle;

	iomodel->FindConstant(&calvinglaw,"md.calving.law");
	switch(calvinglaw){
		case DefaultCalvingEnum:
		case CalvingLevermannEnum:
			break;
		case CalvingVonmisesEnum:
			parameters->AddObject(iomodel->CopyConstantObject("md.calving.min_thickness",CalvingMinthicknessEnum));
			break;
		case CalvingVonmisesADEnum:
			parameters->AddObject(iomodel->CopyConstantObject("md.calving.min_thickness",CalvingMinthicknessEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.calving.num_basins",CalvingNumberofBasinsEnum));

			iomodel->FetchData(&transparam,&M,&N,"md.calving.stress_threshold_groundedice");
         _assert_(M>=1 && N>=1);
         parameters->AddObject(new DoubleVecParam(CalvingADStressThresholdGroundediceEnum,transparam,M));
         xDelete<IssmDouble>(transparam);

         iomodel->FetchData(&transparam,&M,&N,"md.calving.stress_threshold_floatingice");
         _assert_(M>=1 && N>=1);
         parameters->AddObject(new DoubleVecParam(CalvingADStressThresholdFloatingiceEnum,transparam,M));
         xDelete<IssmDouble>(transparam);

			break;
		case CalvingMinthicknessEnum:
			parameters->AddObject(iomodel->CopyConstantObject("md.calving.min_thickness",CalvingMinthicknessEnum));
			break;
		case CalvingHabEnum:
			break;
		case CalvingCrevasseDepthEnum:
			parameters->AddObject(iomodel->CopyConstantObject("md.calving.crevasse_opening_stress",CalvingCrevasseDepthEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.calving.crevasse_threshold",CalvingCrevasseThresholdEnum));
			break;
		case CalvingDev2Enum:
			parameters->AddObject(iomodel->CopyConstantObject("md.calving.height_above_floatation",CalvingHeightAboveFloatationEnum));
			break;
		case CalvingTestEnum:
			iomodel->FindConstant(&interp,"md.timestepping.interp_forcing");
			iomodel->FindConstant(&cycle,"md.timestepping.cycle_forcing");
			iomodel->FetchData(&transparam,&N,&M,"md.calving.speedfactor");
			if(N==1){
				_assert_(M==1);
				parameters->AddObject(new DoubleParam(CalvingTestSpeedfactorEnum,transparam[0]));
         }
         else{
            _assert_(N==2);
            parameters->AddObject(new TransientParam(CalvingTestSpeedfactorEnum,&transparam[0],&transparam[M],interp,cycle,M));
         }
			xDelete<IssmDouble>(transparam);
			iomodel->FetchData(&transparam,&N,&M,"md.calving.independentrate");
			if(N==1){
				_assert_(M==1);
				parameters->AddObject(new DoubleParam(CalvingTestIndependentRateEnum,transparam[0]));
         }
         else{
            _assert_(N==2);
            parameters->AddObject(new TransientParam(CalvingTestIndependentRateEnum,&transparam[0],&transparam[M],interp,cycle,M));
         }
			xDelete<IssmDouble>(transparam);
			break;
		case CalvingParameterizationEnum:
			parameters->AddObject(iomodel->CopyConstantObject("md.calving.use_param",CalvingUseParamEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.calving.min_thickness",CalvingMinthicknessEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.calving.theta",CalvingThetaEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.calving.alpha",CalvingAlphaEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.calving.xoffset",CalvingXoffsetEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.calving.yoffset",CalvingYoffsetEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.calving.vel_lowerbound",CalvingVelLowerboundEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.calving.vel_threshold",CalvingVelThresholdEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.calving.vel_upperbound",CalvingVelUpperboundEnum));
			break;
		case CalvingPollardEnum:
			parameters->AddObject(iomodel->CopyConstantObject("md.calving.rc",CalvingRcEnum));
			break;
		case CalvingCalvingMIPEnum:
			parameters->AddObject(iomodel->CopyConstantObject("md.calving.experiment",CalvingUseParamEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.calving.min_thickness",CalvingMinthicknessEnum));
			break;
		default:
			_error_("Calving law "<<EnumToStringx(calvinglaw)<<" not supported yet");
	}

	/*Get frontal melt parameters*/
	int melt_parameterization;
	iomodel->FindConstant(&melt_parameterization,"md.frontalforcings.parameterization");
	switch(melt_parameterization){
		case FrontalForcingsDefaultEnum:
			break;
		case FrontalForcingsRignotarmaEnum:
			parameters->AddObject(iomodel->CopyConstantObject("md.frontalforcings.num_basins",FrontalForcingsNumberofBasinsEnum));
			/*Retrieve thermal forcing parameters*/
			parameters->AddObject(iomodel->CopyConstantObject("md.frontalforcings.num_params",FrontalForcingsNumberofParamsEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.frontalforcings.num_breaks",FrontalForcingsNumberofBreaksEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.frontalforcings.monthlyvals_numbreaks",FrontalForcingsNumberofMonthBreaksEnum));
         parameters->AddObject(iomodel->CopyConstantObject("md.frontalforcings.ar_order",FrontalForcingsARMAarOrderEnum));
         parameters->AddObject(iomodel->CopyConstantObject("md.frontalforcings.ma_order",FrontalForcingsARMAmaOrderEnum));
         parameters->AddObject(iomodel->CopyConstantObject("md.frontalforcings.arma_timestep",FrontalForcingsARMATimestepEnum));
         iomodel->FetchData(&transparam,&M,&N,"md.frontalforcings.polynomialparams");
         parameters->AddObject(new DoubleMatParam(FrontalForcingsARMApolyparamsEnum,transparam,M,N));
         xDelete<IssmDouble>(transparam);
         iomodel->FetchData(&transparam,&M,&N,"md.frontalforcings.datebreaks");
         parameters->AddObject(new DoubleMatParam(FrontalForcingsARMAdatebreaksEnum,transparam,M,N));
         xDelete<IssmDouble>(transparam);
         iomodel->FetchData(&transparam,&M,&N,"md.frontalforcings.arlag_coefs");
         parameters->AddObject(new DoubleMatParam(FrontalForcingsARMAarlagcoefsEnum,transparam,M,N));
         xDelete<IssmDouble>(transparam);
         iomodel->FetchData(&transparam,&M,&N,"md.frontalforcings.malag_coefs");
         parameters->AddObject(new DoubleMatParam(FrontalForcingsARMAmalagcoefsEnum,transparam,M,N));
         xDelete<IssmDouble>(transparam);
			iomodel->FetchData(&transparam,&M,&N,"md.frontalforcings.monthlyvals_datebreaks");
         parameters->AddObject(new DoubleMatParam(FrontalForcingsARMAmonthdatebreaksEnum,transparam,M,N));
         xDelete<IssmDouble>(transparam);
			iomodel->FetchData(&transparam,&M,&N,"md.frontalforcings.monthlyvals_intercepts");
         parameters->AddObject(new DoubleMatParam(FrontalForcingsARMAmonthinterceptsEnum,transparam,M,N));
         xDelete<IssmDouble>(transparam);
			iomodel->FetchData(&transparam,&M,&N,"md.frontalforcings.monthlyvals_trends");
         parameters->AddObject(new DoubleMatParam(FrontalForcingsARMAmonthtrendsEnum,transparam,M,N));
         xDelete<IssmDouble>(transparam);
			parameters->AddObject(iomodel->CopyConstantObject("md.frontalforcings.isdischargearma",FrontalForcingsIsDischargeARMAEnum));
			/*Retrieve subglacial discharge parameters */
			bool isdischargearma;
			parameters->FindParam(&isdischargearma,FrontalForcingsIsDischargeARMAEnum);
			if(isdischargearma==true){
				parameters->AddObject(iomodel->CopyConstantObject("md.frontalforcings.sd_num_params",FrontalForcingsSdNumberofParamsEnum));
				parameters->AddObject(iomodel->CopyConstantObject("md.frontalforcings.sd_num_breaks",FrontalForcingsSdNumberofBreaksEnum));
      	   parameters->AddObject(iomodel->CopyConstantObject("md.frontalforcings.sd_ar_order",FrontalForcingsSdarOrderEnum));
      	   parameters->AddObject(iomodel->CopyConstantObject("md.frontalforcings.sd_ma_order",FrontalForcingsSdmaOrderEnum));
      	   parameters->AddObject(iomodel->CopyConstantObject("md.frontalforcings.sd_arma_timestep",FrontalForcingsSdARMATimestepEnum));
      	   iomodel->FetchData(&transparam,&M,&N,"md.frontalforcings.sd_polynomialparams");
      	   parameters->AddObject(new DoubleMatParam(FrontalForcingsSdpolyparamsEnum,transparam,M,N));
      	   xDelete<IssmDouble>(transparam);
      	   iomodel->FetchData(&transparam,&M,&N,"md.frontalforcings.sd_datebreaks");
      	   parameters->AddObject(new DoubleMatParam(FrontalForcingsSddatebreaksEnum,transparam,M,N));
      	   xDelete<IssmDouble>(transparam);
      	   iomodel->FetchData(&transparam,&M,&N,"md.frontalforcings.sd_arlag_coefs");
      	   parameters->AddObject(new DoubleMatParam(FrontalForcingsSdarlagcoefsEnum,transparam,M,N));
      	   xDelete<IssmDouble>(transparam);
      	   iomodel->FetchData(&transparam,&M,&N,"md.frontalforcings.sd_malag_coefs");
      	   parameters->AddObject(new DoubleMatParam(FrontalForcingsSdmalagcoefsEnum,transparam,M,N));
      	   xDelete<IssmDouble>(transparam);
				iomodel->FetchData(&transparam,&M,&N,"md.frontalforcings.sd_monthlyfrac");
      	   parameters->AddObject(new DoubleMatParam(FrontalForcingsSdMonthlyFracEnum,transparam,M,N));
      	   xDelete<IssmDouble>(transparam);
			}
			break;
		case FrontalForcingsRignotEnum:
			parameters->AddObject(iomodel->CopyConstantObject("md.frontalforcings.num_basins",FrontalForcingsNumberofBasinsEnum));
			break;
		default:
			_error_("Frontal forcings "<<EnumToStringx(melt_parameterization)<<" not supported yet");
	}
}
/*}}}*/

/*Finite element Analysis*/
void           LevelsetAnalysis::Core(FemModel* femmodel){/*{{{*/

	/*parameters: */
	int  stabilization;
	femmodel->parameters->FindParam(&stabilization,LevelsetStabilizationEnum);

	/*activate formulation: */
	femmodel->SetCurrentConfiguration(LevelsetAnalysisEnum);

	if(stabilization==4){
		solutionsequence_fct(femmodel);
	}
	else{
		solutionsequence_linear(femmodel);
	}
}/*}}}*/
void           LevelsetAnalysis::PreCore(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
ElementVector* LevelsetAnalysis::CreateDVector(Element* element){/*{{{*/
	/*Default, return NULL*/
	return NULL;
}/*}}}*/
ElementMatrix* LevelsetAnalysis::CreateJacobianMatrix(Element* element){/*{{{*/
	/* Jacobian required for the Newton solver */
	_error_("not implemented yet");
}/*}}}*/
ElementMatrix* LevelsetAnalysis::CreateKMatrix(Element* element){/*{{{*/

	if(!element->IsOnBase()) return NULL;
	Element* basalelement = element->SpawnBasalElement();

	/*Intermediaries */
	int  stabilization,dim,domaintype;
	int i,j,k,row, col;
	IssmDouble kappa,factor;
	IssmDouble Jdet, dt, D_scalar;
	IssmDouble h,hx,hy,hz;
	IssmDouble vel,w[3];
	IssmDouble migrationmax;
	IssmDouble* xyz_list = NULL;

	/*Get problem dimension and whether there is moving front or not*/
	basalelement->FindParam(&domaintype,DomainTypeEnum);
	basalelement->FindParam(&stabilization,LevelsetStabilizationEnum);
	switch(domaintype){
		case Domain2DverticalEnum:   dim = 1; break;
		case Domain2DhorizontalEnum: dim = 2; break;
		case Domain3DEnum:           dim = 2; break;
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}
	/*Fetch number of nodes and dof for this finite element*/
	int numnodes    = basalelement->GetNumberOfNodes();

	/*Initialize Element vector and other vectors*/
	ElementMatrix* Ke       = basalelement->NewElementMatrix();
	IssmDouble*    basis    = xNew<IssmDouble>(numnodes);
	IssmDouble*    dbasis   = xNew<IssmDouble>(2*numnodes);

	/*Retrieve all inputs and parameters*/
	basalelement->GetVerticesCoordinates(&xyz_list);
	basalelement->FindParam(&dt,TimesteppingTimeStepEnum);
	basalelement->FindParam(&migrationmax,MigrationMaxEnum);

	h = basalelement->CharacteristicLength();

	Input* mf_vx_input        = NULL;
	Input* mf_vy_input        = NULL;

	/*Load velocities*/
	switch(domaintype){
		case Domain2DverticalEnum:
			mf_vx_input=basalelement->GetInput(MovingFrontalVxEnum); _assert_(mf_vx_input);
			break;
		case Domain2DhorizontalEnum:
			mf_vx_input=basalelement->GetInput(MovingFrontalVxEnum); _assert_(mf_vx_input);
			mf_vy_input=basalelement->GetInput(MovingFrontalVyEnum); _assert_(mf_vy_input);
			break;
		case Domain3DEnum:
			mf_vx_input=basalelement->GetInput(MovingFrontalVxEnum); _assert_(mf_vx_input);
			mf_vy_input=basalelement->GetInput(MovingFrontalVyEnum); _assert_(mf_vy_input);
			break;
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=basalelement->NewGauss(2);
	while(gauss->next()){

		basalelement->JacobianDeterminant(&Jdet,xyz_list,gauss);
		basalelement->NodalFunctions(basis,gauss);
		basalelement->NodalFunctionsDerivatives(dbasis,xyz_list,gauss);
		D_scalar=gauss->weight*Jdet;

		/* Transient */
		if(dt!=0.){
			for(i=0;i<numnodes;i++){
				for(j=0;j<numnodes;j++){
					Ke->values[i*numnodes+j] += D_scalar*basis[j]*basis[i];
				}
			}
			D_scalar=D_scalar*dt;
		}

		/* Levelset speed */
		mf_vx_input->GetInputValue(&w[0], gauss);
		mf_vy_input->GetInputValue(&w[1], gauss);

		/* Apply limiter to the migration rate */		
		vel = 0.;
		for(i=0;i<dim;i++) vel += w[i]*w[i];
		vel = sqrt(vel)+1e-14;
		/* !!NOTE: This is different from the previous version 25838 (and before). The current threshold restrict both advance and retreat velocity. */
		if (vel > migrationmax) {
			for(i=0;i<dim;i++) w[i] = w[i]/vel*migrationmax;
		}

		/*Compute D*/
		for(i=0;i<numnodes;i++){
			for(j=0;j<numnodes;j++){
				for(k=0;k<dim;k++){
					Ke->values[i*numnodes+j] += D_scalar*w[k]*dbasis[k*numnodes+j]*basis[i];
				}
			}
		}

		/* Stabilization */
		vel=0.;
		for(i=0;i<dim;i++) vel+=w[i]*w[i];
		vel=sqrt(vel)+1.e-14;
		switch(stabilization){
			case 0:
				/*Nothing to be done*/
				break;
			case 1:
				/* Artificial Diffusion */
				basalelement->ElementSizes(&hx,&hy,&hz);
				h=sqrt( pow(hx*w[0]/vel,2) + pow(hy*w[1]/vel,2) );
				kappa=h*vel/2.;
				for(i=0;i<numnodes;i++){
					for(j=0;j<numnodes;j++){
						for(k=0;k<dim;k++){
							Ke->values[i*numnodes+j] += D_scalar*kappa*dbasis[k*numnodes+j]*dbasis[k*numnodes+i];
						}
					}
				}
				break;
			case 2:
				  {
					/* Streamline Upwinding */
					mf_vx_input->GetInputAverage(&w[0]);
					mf_vy_input->GetInputAverage(&w[1]);
					vel=sqrt(w[0]*w[0]+w[1]*w[1])+1.e-8;
					IssmDouble tau=h/(2*vel);
					factor = dt*gauss->weight*Jdet*tau;
					for(int i=0;i<numnodes;i++){
						for(int j=0;j<numnodes;j++){
							Ke->values[i*numnodes+j]+=factor*(
										w[0]*dbasis[0*numnodes+i]+w[1]*dbasis[1*numnodes+i])*(w[0]*dbasis[0*numnodes+j]+w[1]*dbasis[1*numnodes+j]);
						}
					}
				  }
				break;
			case 5:{
				/*SUPG*/
				IssmDouble vx,vy;
				mf_vx_input->GetInputAverage(&vx);
				mf_vy_input->GetInputAverage(&vy);
				vel=sqrt(vx*vx+vy*vy)+1.e-8;
				IssmPDouble xi=1.;
				IssmDouble  tau=xi*h/(2*vel);

				/*Mass matrix - part 2*/
				factor = gauss->weight*Jdet*tau;
				for(int i=0;i<numnodes;i++){
					for(int j=0;j<numnodes;j++){
						Ke->values[i*numnodes+j]+=factor*basis[j]*(vx*dbasis[0*numnodes+i]+vy*dbasis[1*numnodes+i]);
					}
				}

				/*Advection matrix - part 2, A*/
				factor = dt*gauss->weight*Jdet*tau;
				for(int i=0;i<numnodes;i++){
					for(int j=0;j<numnodes;j++){
						Ke->values[i*numnodes+j]+=factor*(vx*dbasis[0*numnodes+j]+vy*dbasis[1*numnodes+j])*(vx*dbasis[0*numnodes+i]+vy*dbasis[1*numnodes+i]);
					}
				}

				break;
			}
			case 6:{
				/*SUPG*/
				IssmDouble vx,vy;
				mf_vx_input->GetInputAverage(&vx);
				mf_vy_input->GetInputAverage(&vy);
				vel=sqrt(vx*vx+vy*vy)+1.e-8;
				IssmDouble ECN, K;
				ECN = vel *dt /h;
				K = 1./tanh(ECN) - 1./ECN;
//				if (ECN<1e-6) K = ECN /3.0;

				/*According to Hilmar, xi=K is too large*/
				IssmDouble xi=0.1*K;

				IssmDouble  tau=xi*h/(2*vel);
				Input* levelset_input = NULL;

				IssmDouble kappa;
				IssmDouble p=4, q=4;
				IssmDouble phi[3];

			   levelset_input=basalelement->GetInput(MaskIceLevelsetEnum); _assert_(levelset_input);
				levelset_input->GetInputValue(&phi[0], gauss);

				IssmDouble dphidx=0., dphidy=0.;
				IssmDouble nphi;

				for(int i=0;i<numnodes;i++){
					dphidx += phi[i]*dbasis[0*numnodes+i];
					dphidy += phi[i]*dbasis[1*numnodes+i];
				}
				nphi = sqrt(dphidx*dphidx+dphidy*dphidy);

				if (nphi >= 1) {
					kappa = 1 - 1.0/nphi;
				}
				else {
					kappa = 0.5/M_PI *sin(2*M_PI*nphi)/nphi;
				}

				kappa = kappa * vel / h;

				/*Mass matrix - part 2*/
				factor = gauss->weight*Jdet*tau;
				for(int i=0;i<numnodes;i++){
					for(int j=0;j<numnodes;j++){
						Ke->values[i*numnodes+j]+=factor*basis[j]*(vx*dbasis[0*numnodes+i]+vy*dbasis[1*numnodes+i]);
					}
				}

				/*Advection matrix - part 2, A*/
				factor = dt*gauss->weight*Jdet*tau;
				for(int i=0;i<numnodes;i++){
					for(int j=0;j<numnodes;j++){
						Ke->values[i*numnodes+j]+=factor*(vx*dbasis[0*numnodes+j]+vy*dbasis[1*numnodes+j])*(vx*dbasis[0*numnodes+i]+vy*dbasis[1*numnodes+i]);
					}
				}
				/*Add the pertubation term \nabla\cdot(\kappa*\nabla\phi)*/
				factor = dt*gauss->weight*Jdet*kappa;
				for(int i=0;i<numnodes;i++){
					for(int j=0;j<numnodes;j++){
						for(int k=0;k<dim;k++){
								Ke->values[i*numnodes+j]+= factor*dbasis[k*numnodes+j]*dbasis[k*numnodes+i];
						}
					}
				}

				break;
			}
			default:
				_error_("unknown type of stabilization in LevelsetAnalysis.cpp");
		}
	}

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(basis);
	xDelete<IssmDouble>(dbasis);
	delete gauss;
	if(basalelement->IsSpawnedElement()){basalelement->DeleteMaterials(); delete basalelement;};
	return Ke;
}/*}}}*/
ElementVector* LevelsetAnalysis::CreatePVector(Element* element){/*{{{*/

	if(!element->IsOnBase()) return NULL;
	Element* basalelement = element->SpawnBasalElement();

	/*Intermediaries */
	int         domaintype,stabilization;
	IssmDouble  Jdet,dt;
	IssmDouble  lsf;
	IssmDouble* xyz_list = NULL;

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = basalelement->GetNumberOfNodes();
	basalelement->FindParam(&stabilization,LevelsetStabilizationEnum);

	/*Initialize Element vector*/
	ElementVector* pe = basalelement->NewElementVector();
	basalelement->FindParam(&dt,TimesteppingTimeStepEnum); _assert_(dt>0.);

	/*Initialize basis vector*/
	IssmDouble*    basis = xNew<IssmDouble>(numnodes);
	IssmDouble*    dbasis = NULL;
	if((stabilization==5) |(stabilization == 6)) dbasis= xNew<IssmDouble>(2*numnodes);

	/*Retrieve all inputs and parameters*/
	basalelement->GetVerticesCoordinates(&xyz_list);
	Input* levelset_input = basalelement->GetInput(MaskIceLevelsetEnum); _assert_(levelset_input);
	Input* mf_vx_input    = basalelement->GetInput(MovingFrontalVxEnum); _assert_(mf_vx_input);
	Input* mf_vy_input    = basalelement->GetInput(MovingFrontalVyEnum); _assert_(mf_vy_input);

	IssmDouble h=basalelement->CharacteristicLength();

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=basalelement->NewGauss(2);
	while(gauss->next()){
		basalelement->JacobianDeterminant(&Jdet,xyz_list,gauss);
		basalelement->NodalFunctions(basis,gauss);

		/* old function value */
		levelset_input->GetInputValue(&lsf,gauss);
		IssmDouble factor = Jdet*gauss->weight*lsf;
		for(int i=0;i<numnodes;i++) pe->values[i]+=factor*basis[i];

		if(stabilization==5){ /*SUPG*/
			IssmDouble vx,vy,vel;
			basalelement->NodalFunctionsDerivatives(dbasis,xyz_list,gauss);
			mf_vx_input->GetInputAverage(&vx);
			mf_vy_input->GetInputAverage(&vy);
			vel=sqrt(vx*vx+vy*vy)+1.e-8;
			IssmPDouble xi=1.;
			IssmDouble  tau=xi*h/(2*vel);

			/*Force vector - part 2*/
			factor = Jdet*gauss->weight*lsf;
			for(int i=0;i<numnodes;i++){
				pe->values[i]+=factor*(tau*vx*dbasis[0*numnodes+i]+tau*vy*dbasis[1*numnodes+i]);
			}
		}
		else if (stabilization ==6) {
			IssmDouble vx,vy,vel;
			basalelement->NodalFunctionsDerivatives(dbasis,xyz_list,gauss);
			mf_vx_input->GetInputAverage(&vx);
			mf_vy_input->GetInputAverage(&vy);
			vel=sqrt(vx*vx+vy*vy)+1.e-8;

			IssmDouble ECN, K;
			ECN = vel *dt /h;
			K = 1./tanh(ECN) - 1./ECN;
	//		if (ECN<1e-6) K = ECN /3.0;

			/*According to Hilmar, xi=K is too large*/
			IssmDouble xi=0.1*K;

			IssmDouble  tau=xi*h/(2*vel);

			/*Force vector - part 2*/
			factor = Jdet*gauss->weight*lsf*tau;
			for(int i=0;i<numnodes;i++){
				pe->values[i]+=factor*(vx*dbasis[0*numnodes+i]+vy*dbasis[1*numnodes+i]);
			}
		}
	}

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(basis);
	xDelete<IssmDouble>(dbasis);
	basalelement->FindParam(&domaintype,DomainTypeEnum);
	if(basalelement->IsSpawnedElement()){basalelement->DeleteMaterials(); delete basalelement;};
	delete gauss;

	return pe;
}/*}}}*/
IssmDouble     LevelsetAnalysis::GetDistanceToStraight(IssmDouble* q, IssmDouble* s0, IssmDouble* s1){/*{{{*/
	// returns distance d of point q to straight going through points s0, s1
	// d=|a x b|/|b|
	// with a=q-s0, b=s1-s0

	/* Intermediaries */
	const int dim=2;
	int i;
	IssmDouble a[dim], b[dim];
	IssmDouble norm_b;

	for(i=0;i<dim;i++){
		a[i]=q[i]-s0[i];
		b[i]=s1[i]-s0[i];
	}

	norm_b=0.;
	for(i=0;i<dim;i++)
	 norm_b+=b[i]*b[i];
	norm_b=sqrt(norm_b);
	_assert_(norm_b>0.);

	return fabs(a[0]*b[1]-a[1]*b[0])/norm_b;
}/*}}}*/
void           LevelsetAnalysis::GetSolutionFromInputs(Vector<IssmDouble>* solution,Element* element){/*{{{*/
	element->GetSolutionFromInputsOneDof(solution,MaskIceLevelsetEnum);
}/*}}}*/
void           LevelsetAnalysis::GradientJ(Vector<IssmDouble>* gradient,Element*  element,int control_type,int control_interp,int control_index){/*{{{*/
	_error_("Not implemented yet");
}/*}}}*/
void           LevelsetAnalysis::InputUpdateFromSolution(IssmDouble* solution,Element* element){/*{{{*/

	int domaintype;
	element->FindParam(&domaintype,DomainTypeEnum);
	switch(domaintype){
		case Domain2DhorizontalEnum:
			element->InputUpdateFromSolutionOneDof(solution,MaskIceLevelsetEnum);
			break;
		case Domain3DEnum:
			element->InputUpdateFromSolutionOneDofCollapsed(solution,MaskIceLevelsetEnum);
			break;
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}
}/*}}}*/
void           LevelsetAnalysis::PostProcess(FemModel* femmodel){/*{{{*/

	/*This function is only used by "discrete calving laws" for which we change
	 * the value of the levelset after the advection step (level set equation
	 * solve) based on the law*/

	/*Intermediaries*/
	int  calvinglaw;
	IssmDouble newlevelset[6];
	femmodel->parameters->FindParam(&calvinglaw,CalvingLawEnum);

	/*Apply minimum thickness criterion*/
	if(calvinglaw==CalvingMinthicknessEnum || calvinglaw==CalvingVonmisesEnum || calvinglaw==CalvingParameterizationEnum || calvinglaw==CalvingVonmisesADEnum || calvinglaw==CalvingCalvingMIPEnum){

		IssmDouble mig_max = femmodel->parameters->FindParam(MigrationMaxEnum);
		IssmDouble dt      = femmodel->parameters->FindParam(TimesteppingTimeStepEnum);

		/*Get current distance to terminus*/
		InputDuplicatex(femmodel,MaskIceLevelsetEnum,DistanceToCalvingfrontEnum);
		femmodel->DistanceToFieldValue(MaskIceLevelsetEnum,0,DistanceToCalvingfrontEnum);

		/*Intermediaries*/
		IssmDouble thickness,bed,sealevel,distance,levelset;
		IssmDouble min_thickness = femmodel->parameters->FindParam(CalvingMinthicknessEnum);

		/*Loop over all elements of this partition*/
		for(Object* & object : femmodel->elements->objects){
			Element* element  = xDynamicCast<Element*>(object);

			/*no need to postprocess an ice free element*/
			if(!element->IsIceInElement()) continue;

			int      numnodes     = element->GetNumberOfNodes(); _assert_(numnodes<7);
			Gauss*   gauss        = element->NewGauss();
			Input *H_input        = element->GetInput(ThicknessEnum);              _assert_(H_input);
			Input *b_input        = element->GetInput(BedEnum);                    _assert_(b_input);
			Input *sl_input       = element->GetInput(SealevelEnum);               _assert_(sl_input);
			Input *dis_input      = element->GetInput(DistanceToCalvingfrontEnum); _assert_(dis_input);
			Input *levelset_input = element->GetInput(MaskIceLevelsetEnum);        _assert_(levelset_input);

			/*Potentially constrain nodes of this element*/
			for(int in=0;in<numnodes;in++){
				gauss->GaussNode(element->GetElementType(),in);

				levelset_input->GetInputValue(&levelset,gauss);
				H_input->GetInputValue(&thickness,gauss);
				b_input->GetInputValue(&bed,gauss);
				sl_input->GetInputValue(&sealevel,gauss);
				dis_input->GetInputValue(&distance,gauss);

				if(thickness<min_thickness && bed<sealevel && fabs(distance)<mig_max*dt && levelset<0){
					newlevelset[in] = +400.; //Arbitrary > 0 number (i.e. deactivate this node)
				}
				else{
					newlevelset[in] = levelset;
				}
			}
			element->AddInput(MaskIceLevelsetEnum,&newlevelset[0],element->GetElementType());
			delete gauss;
		}
	}
}/*}}}*/
void           LevelsetAnalysis::UpdateConstraints(FemModel* femmodel){/*{{{*/

	/*Intermediaries*/
	int  calvinglaw;
	femmodel->parameters->FindParam(&calvinglaw,CalvingLawEnum);
	IssmDouble mig_max = femmodel->parameters->FindParam(MigrationMaxEnum);
	IssmDouble dt      = femmodel->parameters->FindParam(TimesteppingTimeStepEnum);

   /* Get current distance to terminus
	 * Only do this if necessary, PostProcess is already doing it for a few calving law
	 * Do not repeat the process is this function is particularly slow*/
	bool computedistance = true;
	if(
				calvinglaw==CalvingMinthicknessEnum ||
				calvinglaw==CalvingVonmisesEnum ||
				calvinglaw==CalvingParameterizationEnum ||
				calvinglaw==CalvingVonmisesADEnum ||
				calvinglaw==CalvingCalvingMIPEnum){
		int step;
		femmodel->parameters->FindParam(&step,StepEnum);
		if(step>1){
			computedistance = false;
		}
	}
	if(computedistance){
		InputDuplicatex(femmodel,MaskIceLevelsetEnum,DistanceToCalvingfrontEnum);
		femmodel->DistanceToFieldValue(MaskIceLevelsetEnum,0,DistanceToCalvingfrontEnum);
	}

   if(calvinglaw==CalvingHabEnum){

		/*Intermediaries*/
		IssmDouble  thickness,water_depth,distance,hab_fraction;

		/*Loop over all elements of this partition*/
		for(Object* & object : femmodel->elements->objects){
			Element* element  = xDynamicCast<Element*>(object);

			IssmDouble rho_ice   = element->FindParam(MaterialsRhoIceEnum);
			IssmDouble rho_water = element->FindParam(MaterialsRhoSeawaterEnum);

			int      numnodes           = element->GetNumberOfNodes();
			Gauss*   gauss              = element->NewGauss();
			Input*   H_input            = element->GetInput(ThicknessEnum); _assert_(H_input);
			Input*   bed_input          = element->GetInput(BedEnum); _assert_(bed_input);
			Input*   hab_fraction_input = element->GetInput(CalvingHabFractionEnum); _assert_(hab_fraction_input);
			Input*   dis_input           = element->GetInput(DistanceToCalvingfrontEnum); _assert_(dis_input);

			/*Potentially constrain nodes of this element*/
			for(int in=0;in<numnodes;in++){
				gauss->GaussNode(element->GetElementType(),in);
				Node* node=element->GetNode(in);
				if(!node->IsActive()) continue;

				H_input->GetInputValue(&thickness,gauss);
				bed_input->GetInputValue(&water_depth,gauss);
				dis_input->GetInputValue(&distance,gauss);
				hab_fraction_input->GetInputValue(&hab_fraction,gauss);

				if(thickness<((rho_water/rho_ice)*(1+hab_fraction)*-water_depth) && fabs(distance)<mig_max*dt){
					node->ApplyConstraint(0,+1.);
				}
				else {
					/* no ice, set no spc */
					node->DofInFSet(0);
				}
			}
			delete gauss;
		}
	}
   else if(calvinglaw==CalvingCrevasseDepthEnum){

		/*Intermediaries*/
		IssmDouble  levelset,crevassedepth,bed,surface_crevasse,thickness,surface;
		IssmDouble* constraint_nodes = NULL;

		/*Get the DistanceToCalvingfront*/
		InputDuplicatex(femmodel,MaskIceLevelsetEnum,DistanceToCalvingfrontEnum);
		femmodel->DistanceToFieldValue(MaskIceLevelsetEnum,0,DistanceToCalvingfrontEnum);

		/*Vector of size number of nodes*/
      int numnodes      = femmodel->nodes->NumberOfNodes();
      int localmasters  = femmodel->nodes->NumberOfNodesLocal();
      Vector<IssmDouble>* vec_constraint_nodes = vec_constraint_nodes=new Vector<IssmDouble>(localmasters,numnodes);

		IssmDouble crevasse_threshold = femmodel->parameters->FindParam(CalvingCrevasseThresholdEnum);

		for(Object* & object : femmodel->elements->objects){
			Element* element   = xDynamicCast<Element*>(object);
			int      numnodes  = element->GetNumberOfNodes();
			Gauss*   gauss     = element->NewGauss();

			Input*   crevassedepth_input    = element->GetInput(CrevasseDepthEnum); _assert_(crevassedepth_input);
			Input*   bed_input              = element->GetInput(BedEnum); _assert_(bed_input);
			Input*   surface_crevasse_input = element->GetInput(SurfaceCrevasseEnum); _assert_(surface_crevasse_input);
			Input*   thickness_input        = element->GetInput(ThicknessEnum); _assert_(thickness_input);
			Input*   surface_input          = element->GetInput(SurfaceEnum); _assert_(surface_input);

			/*First, look at ice front and figure out if any of the nodes will be calved*/
			if(element->IsIcefront()){
				for(int in=0;in<numnodes;in++){
					gauss->GaussNode(element->GetElementType(),in);
					Node* node=element->GetNode(in);
					if(!node->IsActive()) continue;

					crevassedepth_input->GetInputValue(&crevassedepth,gauss);
					bed_input->GetInputValue(&bed,gauss);
					surface_crevasse_input->GetInputValue(&surface_crevasse,gauss);
					thickness_input->GetInputValue(&thickness,gauss);
					surface_input->GetInputValue(&surface,gauss);

					if((surface_crevasse>surface || crevassedepth>(crevasse_threshold*thickness)-1e-10) && bed<0.){
						vec_constraint_nodes->SetValue(node->Pid(),1.0,INS_VAL);
					}
				}
			}
			delete gauss;
		}

		/*Assemble vector and serialize: */
		vec_constraint_nodes->Assemble();
      femmodel->GetLocalVectorWithClonesNodes(&constraint_nodes,vec_constraint_nodes);

		int nflipped=1;
		while(nflipped){
			int local_nflipped=0;
			for(Object* & object : femmodel->elements->objects){
				Element* element  = xDynamicCast<Element*>(object);
				int      numnodes = element->GetNumberOfNodes();

				Input *levelset_input         = element->GetInput(DistanceToCalvingfrontEnum); _assert_(levelset_input);
				Input *crevassedepth_input    = element->GetInput(CrevasseDepthEnum);          _assert_(crevassedepth_input);
				Input *bed_input              = element->GetInput(BedEnum);                    _assert_(bed_input);
				Input *surface_crevasse_input = element->GetInput(SurfaceCrevasseEnum);        _assert_(surface_crevasse_input);
				Input *thickness_input        = element->GetInput(ThicknessEnum);              _assert_(thickness_input);
				Input *surface_input          = element->GetInput(SurfaceEnum);                _assert_(surface_input);

				/*Is this element connected to a node that should be calved*/
				bool isconnected = false;
				for(int in=0;in<numnodes;in++){
					Node* node=element->GetNode(in);
					if(constraint_nodes[node->Lid()]>0.){
						isconnected = true;
						break;
					}
				}

				/*Check status if connected*/
				if(isconnected){
					Gauss* gauss = element->NewGauss();
					for(int in=0;in<numnodes;in++){
						gauss->GaussNode(element->GetElementType(),in);
						Node* node=element->GetNode(in);
						levelset_input->GetInputValue(&levelset,gauss);
						crevassedepth_input->GetInputValue(&crevassedepth,gauss);
						bed_input->GetInputValue(&bed,gauss);
						surface_crevasse_input->GetInputValue(&surface_crevasse,gauss);
						thickness_input->GetInputValue(&thickness,gauss);
						surface_input->GetInputValue(&surface,gauss);

                  /*FIXME: not sure about levelset<0. && fabs(levelset)>-mig_max*dt! SHould maybe be distance<mig_max*dt*/
                  if((surface_crevasse>surface || crevassedepth>(crevasse_threshold*thickness-1e-10)) && bed<0. && levelset<0. && levelset>-mig_max*dt && constraint_nodes[node->Lid()]==0.){
							local_nflipped++;
							vec_constraint_nodes->SetValue(node->Pid(),1.0,INS_VAL);
						}
					}
					delete gauss;
				}
			}

			/*Count how many new nodes were found*/
			ISSM_MPI_Allreduce(&local_nflipped,&nflipped,1,ISSM_MPI_INT,ISSM_MPI_SUM,IssmComm::GetComm());
			_printf0_("Found "<<nflipped<<" to flip\n");

			/*Assemble and serialize flag vector*/
			vec_constraint_nodes->Assemble();
			xDelete<IssmDouble>(constraint_nodes);
         femmodel->GetLocalVectorWithClonesNodes(&constraint_nodes,vec_constraint_nodes);
		}

		/*Free resources:*/
		delete vec_constraint_nodes;

		/*Contrain the nodes that will be calved*/
		for(Object* & object : femmodel->elements->objects){
			Element* element  = xDynamicCast<Element*>(object);
			int      numnodes = element->GetNumberOfNodes();
			Gauss*   gauss    = element->NewGauss();
			/*Potentially constrain nodes of this element*/
			for(int in=0;in<numnodes;in++){
				gauss->GaussNode(element->GetElementType(),in);
				Node* node=element->GetNode(in);
				if(!node->IsActive()) continue;

				if(constraint_nodes[node->Lid()]>0.){
					node->ApplyConstraint(0,+1.);
				}
				else {
					/* no ice, set no spc */
					node->DofInFSet(0);
				}
			}
			delete gauss;
		}
		xDelete<IssmDouble>(constraint_nodes);
	}

	/*Default, do nothing*/
	return;
}/*}}}*/

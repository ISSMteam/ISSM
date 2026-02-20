/*!\file Element.cpp
 * \brief: implementation of the Element object
 */
/*Headers:*/
/*{{{*/
#ifdef HAVE_CONFIG_H
#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "../classes.h"
#include "../../shared/shared.h"
#include "../../modules/SurfaceMassBalancex/SurfaceMassBalancex.h"
#include "../Inputs/BoolInput.h"
#include "../Inputs/TransientInput.h"
#include "../Inputs/ElementInput.h"
#include "../Inputs/PentaInput.h"
#include "../Inputs/DatasetInput.h"
#include "../Inputs/ControlInput.h"
#include "../Inputs/ArrayInput.h"
#include "../Inputs/IntArrayInput.h"
/*}}}*/
#define MAXVERTICES 6 /*Maximum number of vertices per element, currently Penta, to avoid dynamic mem allocation*/

#ifdef _HAVE_SEMIC_
/* SEMIC prototype {{{*/
extern "C" void run_semic_(IssmDouble *sf_in, IssmDouble *rf_in, IssmDouble *swd_in, IssmDouble *lwd_in, IssmDouble *wind_in, IssmDouble *sp_in, IssmDouble *rhoa_in,
			IssmDouble *qq_in, IssmDouble *tt_in, IssmDouble *tsurf_out, IssmDouble *smb_out, IssmDouble *saccu_out, IssmDouble *smelt_out);

extern "C" void run_semic_transient_(int *nx, int *ntime, int *nloop, 
			IssmDouble *sf_in, IssmDouble *rf_in, IssmDouble *swd_in, 
			IssmDouble *lwd_in, IssmDouble *wind_in, IssmDouble *sp_in, IssmDouble *rhoa_in,
			IssmDouble *qq_in, IssmDouble *tt_in, IssmDouble *tsurf_in, IssmDouble *qmr_in,
			IssmDouble *tstic,
			IssmDouble *hcrit, IssmDouble *rcrit,
			IssmDouble *mask, IssmDouble *hice, IssmDouble *hsnow,
			IssmDouble *albedo_in, IssmDouble *albedo_snow_in,
			int *alb_scheme, IssmDouble *alb_smax, IssmDouble *alb_smin, IssmDouble *albi, IssmDouble *albl,
			IssmDouble *Tamp, 
			IssmDouble *tmin, IssmDouble *tmax, IssmDouble *tmid, IssmDouble *mcrit, IssmDouble *wcrit, IssmDouble *tau_a, IssmDouble* tau_f, IssmDouble *afac, bool *verbose,
			IssmDouble *tsurf_out, IssmDouble *smb_out, IssmDouble *smbi_out, IssmDouble *smbs_out, IssmDouble *saccu_out, IssmDouble *smelt_out, IssmDouble *refr_out, IssmDouble *albedo_out, IssmDouble *albedo_snow_out, IssmDouble *hsnow_out, IssmDouble *hice_out, IssmDouble *qmr_out, IssmDouble *runoff_out, IssmDouble *subl_out);
#endif
// _HAVE_SEMIC_
/*}}}*/
/*Constructors/destructor/copy*/
Element::Element(){/*{{{*/
	this->id  = -1;
	this->sid = -1;
	this->lid = -1;
	this->inputs    = NULL;
	this->nodes      = NULL;
	this->vertices   = NULL;
	this->material   = NULL;
	this->parameters = NULL;
	this->element_type_list=NULL;
}/*}}}*/
Element::~Element(){/*{{{*/
	xDelete<int>(element_type_list);
}
/*}}}*/

/*Other*/
bool       Element::AnyFSet(){/*{{{*/

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = this->GetNumberOfNodes();

	for(int i=0;i<numnodes;i++){
		if(nodes[i]->FSize()) return true;
	}
	return false;
}/*}}}*/
void       Element::ArmaProcess(bool isstepforarma,int arorder,int maorder,int numparams,int numbreaks,IssmDouble tstep_arma,IssmDouble* polyparams,IssmDouble* arlagcoefs,IssmDouble* malagcoefs,IssmDouble* datebreaks,bool isfieldstochastic,int enum_type){/*{{{*/
   const int numvertices = this->GetNumberOfVertices();
	int         numperiods = numbreaks+1; 
   int         basinid,M,N,arenum_type,maenum_type,basinenum_type,noiseenum_type,outenum_type,indperiod;
   IssmDouble  time,dt,starttime,noiseterm;
   IssmDouble* arlagcoefs_basin     = xNew<IssmDouble>(arorder);
   IssmDouble* malagcoefs_basin     = xNew<IssmDouble>(maorder);
   IssmDouble* datebreaks_basin     = xNew<IssmDouble>(numbreaks);
   IssmDouble* polyparams_basin     = xNew<IssmDouble>(numperiods*numparams);
   IssmDouble* varlist              = xNew<IssmDouble>(numvertices);
   IssmDouble* sumpoly              = xNewZeroInit<IssmDouble>(arorder+1);
   IssmDouble* valuesautoregression = NULL;
   IssmDouble* valuesmovingaverage  = NULL;
   Input*      noiseterm_input      = NULL;

   /*Get field-specific enums*/
   switch(enum_type){
      case(SMBarmaEnum):
         arenum_type    = SmbValuesAutoregressionEnum;
         maenum_type    = SmbValuesMovingaverageEnum;
         basinenum_type = SmbBasinsIdEnum;
         noiseenum_type = SmbARMANoiseEnum;
         outenum_type   = SmbMassBalanceEnum;
         break;
      case(FrontalForcingsRignotarmaEnum):
         arenum_type    = ThermalforcingValuesAutoregressionEnum;
         maenum_type    = ThermalforcingValuesMovingaverageEnum;
         basinenum_type = FrontalForcingsBasinIdEnum;
         noiseenum_type = ThermalforcingARMANoiseEnum;
         outenum_type   = ThermalForcingEnum;
         break;
		case(BasalforcingsDeepwaterMeltingRatearmaEnum):
         arenum_type    = BasalforcingsDeepwaterMeltingRateValuesAutoregressionEnum;
         maenum_type    = BasalforcingsDeepwaterMeltingRateValuesMovingaverageEnum;
         basinenum_type = BasalforcingsLinearBasinIdEnum;
         noiseenum_type = BasalforcingsDeepwaterMeltingRateNoiseEnum;
         outenum_type   = BasalforcingsSpatialDeepwaterMeltingRateEnum;
         break;
		case(FrontalForcingsSubglacialDischargearmaEnum):
         arenum_type    = SubglacialdischargeValuesAutoregressionEnum;
         maenum_type    = SubglacialdischargeValuesMovingaverageEnum;
         basinenum_type = FrontalForcingsBasinIdEnum;
         noiseenum_type = SubglacialdischargeARMANoiseEnum;
         outenum_type   = FrontalForcingsSubglacialDischargeEnum;
         break;
		case(HydrologyarmapwEnum):
         arenum_type    = WaterPressureValuesAutoregressionEnum;
         maenum_type    = WaterPressureValuesMovingaverageEnum;
         basinenum_type = HydrologyBasinsIdEnum;
         noiseenum_type = FrictionWaterPressureNoiseEnum;
         outenum_type   = WaterPressureArmaPerturbationEnum;
         break;
	}

	/*Get time parameters*/
   this->parameters->FindParam(&time,TimeEnum);
   this->parameters->FindParam(&dt,TimesteppingTimeStepEnum);
   this->parameters->FindParam(&starttime,TimesteppingStartTimeEnum);

   /*Get basin coefficients*/
   this->GetInputValue(&basinid,basinenum_type);
	for(int i=0;i<arorder;i++) arlagcoefs_basin[i]   = arlagcoefs[basinid*arorder+i];
	for(int i=0;i<maorder;i++) malagcoefs_basin[i]   = malagcoefs[basinid*maorder+i];
	for(int i=0;i<numparams;i++){
		for(int j=0;j<numperiods;j++) polyparams_basin[i*numperiods+j] = polyparams[basinid*numparams*numperiods+i*numperiods+j];
	}
	if(numbreaks>0){
		for(int i=0;i<numbreaks;i++) datebreaks_basin[i] = datebreaks[basinid*numbreaks+i];
	}

	/*Compute terms from polynomial parameters from arorder steps back to present*/
	IssmDouble telapsed_break;
	IssmDouble tatstep;
	for(int s=0;s<arorder+1;s++){
		tatstep = time-s*tstep_arma;
		if(numbreaks>0){
			/*Find index of tatstep compared to the breakpoints*/
			indperiod = 0;
			for(int i=0;i<numbreaks;i++){
				if(tatstep>=datebreaks_basin[i]) indperiod = i+1;
			}
			/*Compute polynomial with parameters of indperiod*/
			if(indperiod==0) telapsed_break = tatstep-starttime;
			else             telapsed_break = tatstep-datebreaks_basin[indperiod-1];
			for(int j=0;j<numparams;j++)   sumpoly[s] = sumpoly[s]+polyparams_basin[indperiod+j*numperiods]*pow(telapsed_break,j);
		}
		else for(int j=0;j<numparams;j++) sumpoly[s] = sumpoly[s]+polyparams_basin[j*numperiods]*pow(tatstep-starttime,j);
	}

	/*Initialze autoregressive and moving-average values at first time step*/
	if(time<=starttime+dt){
		IssmDouble* initvaluesautoregression = xNewZeroInit<IssmDouble>(numvertices*arorder);
		IssmDouble* initvaluesmovingaverage  = xNewZeroInit<IssmDouble>(numvertices*maorder);
		for(int i=0;i<numvertices*arorder;i++) initvaluesautoregression[i]=polyparams_basin[0];
      this->inputs->SetArrayInput(arenum_type,this->lid,initvaluesautoregression,numvertices*arorder);
      this->inputs->SetArrayInput(maenum_type,this->lid,initvaluesmovingaverage,numvertices*maorder);
      xDelete<IssmDouble>(initvaluesautoregression);
      xDelete<IssmDouble>(initvaluesmovingaverage);
	}

   /*Get noise, autoregressive terms, moving-average terms*/
	if(isfieldstochastic){
      noiseterm_input = this->GetInput(noiseenum_type);
      Gauss* gauss = this->NewGauss();
      noiseterm_input->GetInputValue(&noiseterm,gauss);
      delete gauss;
   }
   else noiseterm = 0.0;
   this->inputs->GetArray(arenum_type,this->lid,&valuesautoregression,&M);
   this->inputs->GetArray(maenum_type,this->lid,&valuesmovingaverage,&M);

	/*If not ARMA model timestep: take the old values of variable*/
   if(isstepforarma==false){
      for(int i=0;i<numvertices;i++) varlist[i]=valuesautoregression[i];
   }
   /*If ARMA model timestep: compute variable values on vertices from ARMA*/
   else{
      for(int v=0;v<numvertices;v++){

         /*Compute autoregressive term*/
         IssmDouble autoregressionterm=0.;
         for(int ii=0;ii<arorder;ii++) autoregressionterm += arlagcoefs_basin[ii]*(valuesautoregression[v+ii*numvertices]-sumpoly[ii+1]);
			/*Compute moving-average term*/
         IssmDouble movingaverageterm=0.;
         for(int ii=0;ii<maorder;ii++) movingaverageterm  += malagcoefs_basin[ii]*valuesmovingaverage[v+ii*numvertices];

			/*Stochastic variable value*/
         varlist[v] = sumpoly[0]+autoregressionterm+movingaverageterm+noiseterm;

			/*Impose zero-bound*/
			if(outenum_type == ThermalForcingEnum || outenum_type == FrontalForcingsSubglacialDischargeEnum) varlist[v] = max(varlist[v],0.0);

		}

      /*Update autoregression and moving-average values*/
      IssmDouble* temparrayar = xNew<IssmDouble>(numvertices*arorder);
      IssmDouble* temparrayma = xNew<IssmDouble>(numvertices*maorder);
      /*Assign newest values and shift older values*/
      for(int i=0;i<numvertices;i++) temparrayar[i] = varlist[i];
      for(int i=0;i<numvertices;i++) temparrayma[i] = noiseterm;
      for(int i=0;i<(arorder-1)*numvertices;i++) temparrayar[i+numvertices] = valuesautoregression[i];
      for(int i=0;i<(maorder-1)*numvertices;i++) temparrayma[i+numvertices] = valuesmovingaverage[i];
		this->inputs->SetArrayInput(arenum_type,this->lid,temparrayar,numvertices*arorder);
      this->inputs->SetArrayInput(maenum_type,this->lid,temparrayma,numvertices*maorder);
      xDelete<IssmDouble>(temparrayar);
      xDelete<IssmDouble>(temparrayma);
   }

   /*Add input to element*/
   this->AddInput(outenum_type,varlist,P1Enum);

   /*Cleanup*/
   xDelete<IssmDouble>(arlagcoefs_basin);
   xDelete<IssmDouble>(malagcoefs_basin);
   xDelete<IssmDouble>(datebreaks_basin);
   xDelete<IssmDouble>(polyparams_basin);
   xDelete<IssmDouble>(sumpoly);
   xDelete<IssmDouble>(varlist);
   xDelete<IssmDouble>(valuesautoregression);
   xDelete<IssmDouble>(valuesmovingaverage);
}

/*}}}*/
void       Element::BasinLinearFloatingiceMeltingRate(IssmDouble* deepwaterel,IssmDouble* upperwatermelt,IssmDouble* upperwaterel,IssmDouble* perturbation){/*{{{*/

	const int NUM_VERTICES = this->GetNumberOfVertices();

	int basinid;
   IssmDouble alpha;
   IssmDouble base[MAXVERTICES];
   IssmDouble values[MAXVERTICES];
   IssmDouble deepwatermelt[MAXVERTICES];

	/*Get element-specific values*/
	this->GetInputValue(&basinid,BasalforcingsLinearBasinIdEnum);
	this->GetInputListOnVertices(&base[0],BaseEnum);
   this->GetInputListOnVertices(&deepwatermelt[0],BasalforcingsSpatialDeepwaterMeltingRateEnum);

	/*Compute melt rate at vertices according to basin-specific values of input arguments*/
   for(int i=0;i<NUM_VERTICES;i++){
		if(base[i]>=upperwaterel[basinid]){
         values[i]=upperwatermelt[basinid];
      }
      else if (base[i]<deepwaterel[basinid]){
         values[i]=deepwatermelt[i];
      }
      else{
         alpha = (base[i]-upperwaterel[basinid])/(deepwaterel[basinid]-upperwaterel[basinid]);
         values[i]=deepwatermelt[i]*alpha+(1.-alpha)*upperwatermelt[basinid];
      }
		if(xIsNan<IssmDouble>(values[i])) _error_("NaN found in melt");
		if(xIsInf<IssmDouble>(values[i])) _error_("Inf found in melt");
		if(fabs(values[i])>1.e+10) _error_("melt exceeds 1.e+10");
   }

   this->AddInput(BasalforcingsFloatingiceMeltingRateEnum,&values[0],P1Enum);
}/*}}}*/
void       Element::CalvingRateToVector(){/*{{{*/
	this->CalvingRateToVector(false);
}/*}}}*/
void       Element::CalvingRateToVector(bool isvelvector){/*{{{*/

	/*We are provided a calving rate, figure out the x/y components*/
	int         dim,domaintype;
	IssmDouble  vx,vy,vel,dphidx,dphidy,dphi,c;
	IssmDouble  calvingratex[MAXVERTICES];
	IssmDouble  calvingratey[MAXVERTICES];

	/*Get problem dimension and whether there is moving front or not*/
	this->FindParam(&domaintype,DomainTypeEnum);

	switch(domaintype){
		case Domain2DverticalEnum:   dim = 1; break;
		case Domain2DhorizontalEnum: dim = 2; break;
		case Domain3DEnum:           dim = 2; break;
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}
	if(dim==1) _error_("not implemented in 1D...");

	Input *calvingrate_input = this->GetInput(CalvingCalvingrateEnum);     _assert_(calvingrate_input);
	Input *vx_input          = this->GetInput(VxEnum);  _assert_(vx_input);
	Input *vy_input          = this->GetInput(VyEnum); _assert_(vy_input);
	Input *lsf_slopex_input  = this->GetInput(LevelsetfunctionSlopeXEnum); _assert_(lsf_slopex_input);
	Input *lsf_slopey_input  = this->GetInput(LevelsetfunctionSlopeYEnum); _assert_(lsf_slopey_input);

	/*Allocate arrays*/
	const int NUM_VERTICES = this->GetNumberOfVertices();

	/* Start looping on the number of vertices: */
	Gauss* gauss=this->NewGauss();
	for (int iv=0;iv<NUM_VERTICES;iv++){
		gauss->GaussVertex(iv);

      calvingrate_input->GetInputValue(&c,gauss);
		vx_input->GetInputValue(&vx,gauss);
		vy_input->GetInputValue(&vy,gauss);
		lsf_slopex_input->GetInputValue(&dphidx,gauss);
		lsf_slopey_input->GetInputValue(&dphidy,gauss);

		vel=sqrt(vx*vx + vy*vy) + 1e-14;
		dphi=sqrt(dphidx*dphidx+dphidy*dphidy)+ 1e-14;

		if(isvelvector){
			/*Velocity direction*/
			calvingratex[iv] = c*vx/vel;
         calvingratey[iv] = c*vy/vel;
		}
		else{
			/*Lelvelset gradient (perpendicular to ice front)*/
         calvingratex[iv] = c*dphidx/dphi;
         calvingratey[iv] = c*dphidy/dphi;
		}
	}

	/*Add to inputs*/
	this->AddInput(CalvingratexEnum,&calvingratex[0],P1DGEnum);
	this->AddInput(CalvingrateyEnum,&calvingratey[0],P1DGEnum);
	//this->AddInput(CalvingCalvingrateEnum,&calvingrate[0],P1DGEnum); /*Do not change calving rate, that's our input!*/

	/*Clean up and return*/
	delete gauss;
}/*}}}*/
void       Element::CalvingSetZeroRate(){/*{{{*/

	/*Set calving rate as 0, this is probably because we are  dealing with discrete calving laws*/
	IssmDouble  calvingratex[MAXVERTICES] = {0.};
	IssmDouble  calvingratey[MAXVERTICES] = {0.};
	IssmDouble  calvingrate[MAXVERTICES]  = {0.};
	this->AddInput(CalvingratexEnum,&calvingratex[0],P1DGEnum);
	this->AddInput(CalvingrateyEnum,&calvingratey[0],P1DGEnum);
	this->AddInput(CalvingCalvingrateEnum,&calvingrate[0],P1DGEnum);
}/*}}}*/
void       Element::ComputeLambdaS(){/*{{{*/

	/*Intermediaries*/
	IssmDouble vx,vy,vz,vmag;
	IssmDouble dvx[3],dvy[3],dvz[3],dvmag[3];
	IssmDouble eps[3][3],epseff,epsprime;
	IssmDouble lambdas[MAXVERTICES];
	int         dim;
	IssmDouble *xyz_list = NULL;

	/*Retrieve all inputs we will be needing: */
	this->GetVerticesCoordinates(&xyz_list);
	parameters->FindParam(&dim,DomainDimensionEnum);
	Input* vx_input=this->GetInput(VxEnum); _assert_(vx_input);
	Input* vy_input=this->GetInput(VyEnum); _assert_(vy_input);
	Input* vz_input=NULL;
	if(dim==3){vz_input=this->GetInput(VzEnum); _assert_(vz_input);}

	/*Allocate arrays*/
	const int NUM_VERTICES = this->GetNumberOfVertices();

	/* Start looping on the number of vertices: */
	Gauss* gauss=this->NewGauss();
	for (int iv=0;iv<NUM_VERTICES;iv++){
		gauss->GaussVertex(iv);

		/*Get velocity derivatives in all directions*/
		_assert_(dim>1);
		_assert_(vx_input);
		vx_input->GetInputValue(&vx,gauss);
		vx_input->GetInputDerivativeValue(&dvx[0],xyz_list,gauss);
		_assert_(vy_input);
		vy_input->GetInputValue(&vy,gauss);
		vy_input->GetInputDerivativeValue(&dvy[0],xyz_list,gauss);
		if(dim==3){
			_assert_(vz_input);
			vz_input->GetInputValue(&vz,gauss);
			vz_input->GetInputDerivativeValue(&dvz[0],xyz_list,gauss);
		}
		else{
			vz = 0.;
			dvz[0] = 0.; dvz[1] = 0.; dvz[2] = 0.;
			dvx[2]= 0.;
			dvy[2]= 0.;
		}
		/*Calculate velocity magnitude and its derivative*/
		vmag = sqrt(vx*vx+vy*vy+vz*vz);
		if(vmag<1e-12){
			vmag=1e-12;
			dvmag[0]=0;
			dvmag[1]=0;
			dvmag[2]=0;
		}
		else{
			dvmag[0]=1./(2*sqrt(vmag))*(2*vx*dvx[0]+2*vy*dvy[0]+2*vz*dvz[0]);
			dvmag[1]=1./(2*sqrt(vmag))*(2*vx*dvx[1]+2*vy*dvy[1]+2*vz*dvz[1]);
			dvmag[2]=1./(2*sqrt(vmag))*(2*vx*dvx[2]+2*vy*dvy[2]+2*vz*dvz[2]);
		}
		/*Build strain rate tensor*/
		eps[0][0] = dvx[0];             eps[0][1] = .5*(dvx[1]+dvy[0]); eps[0][2] = .5*(dvx[2]+dvz[0]);
		eps[1][0] = .5*(dvx[1]+dvy[0]); eps[1][1] = dvy[1];             eps[1][2] = .5*(dvy[2]+dvz[1]);
		eps[2][0] = .5*(dvx[2]+dvz[0]); eps[2][1] = .5*(dvy[2]+dvz[1]); eps[2][2] = dvz[2];

		/*effective strain rate*/
		epseff = 0.;
		/* eps_eff^2 = exx^2 + eyy^2 + exy^2 + exz^2 + eyz^2 + exx*eyy */
		epseff = sqrt(eps[0][0]*eps[0][0] + eps[1][1]*eps[1][1] + eps[0][1]*eps[0][1] +  eps[0][2]*eps[0][2] + eps[1][2]*eps[1][2] + eps[0][0]*eps[1][1]);

		EstarStrainrateQuantities(&epsprime,vx,vy,vz,vmag,&dvx[0],&dvy[0],&dvz[0],&dvmag[0]);
		lambdas[iv]=EstarLambdaS(epseff,epsprime);
	}

	/*Add Stress tensor components into inputs*/
	this->AddInput(LambdaSEnum,&lambdas[0],P1Enum);

	/*Clean up and return*/
	delete gauss;
	xDelete<IssmDouble>(xyz_list);
}/*}}}*/
void       Element::ComputeNewDamage(){/*{{{*/

	IssmDouble *xyz_list=NULL;
	IssmDouble  eps_xx,eps_xy,eps_yy,eps_xz,eps_yz,eps_zz,eps_eff;
	IssmDouble  epsmin=1.e-27;
	IssmDouble  eps_0,kappa,sigma_0,B,D,n,envelopeD;
	int         dim,counter=0;
	IssmDouble  k1,k2,threshold=1.e-12;

	/* Retrieve parameters */
	this->GetVerticesCoordinates(&xyz_list);
	this->ComputeStrainRate();
	parameters->FindParam(&dim,DomainDimensionEnum);
	parameters->FindParam(&kappa,DamageKappaEnum);
	parameters->FindParam(&sigma_0,DamageStressThresholdEnum);

	/* Retrieve inputs */
	Input* eps_xx_input=this->GetInput(StrainRatexxEnum); _assert_(eps_xx_input);
	Input* eps_yy_input=this->GetInput(StrainRateyyEnum); _assert_(eps_yy_input);
	Input* eps_xy_input=this->GetInput(StrainRatexyEnum); _assert_(eps_xy_input);
	Input* eps_xz_input=NULL;
	Input* eps_yz_input=NULL;
	Input* eps_zz_input=NULL;
	if(dim==3){
		eps_xz_input=this->GetInput(StrainRatexzEnum); _assert_(eps_xz_input);
		eps_yz_input=this->GetInput(StrainRateyzEnum); _assert_(eps_yz_input);
		eps_zz_input=this->GetInput(StrainRatezzEnum); _assert_(eps_zz_input);
	}

	/* Fetch number of nodes and allocate output*/
	int numnodes = this->GetNumberOfNodes();
	IssmDouble* newD = xNew<IssmDouble>(numnodes);

	/* Retrieve domain-dependent inputs */
	Input* n_input=this->GetInput(MaterialsRheologyNEnum); _assert_(n_input);
	Input* damage_input = NULL;
	Input* B_input = NULL;
	int domaintype;
	parameters->FindParam(&domaintype,DomainTypeEnum);
	if(domaintype==Domain2DhorizontalEnum){
		damage_input = this->GetInput(DamageDbarOldEnum);  _assert_(damage_input);
		B_input=this->GetInput(MaterialsRheologyBbarEnum); _assert_(B_input);
	}
	else{
		damage_input = this->GetInput(DamageDOldEnum);   _assert_(damage_input);
		B_input=this->GetInput(MaterialsRheologyBEnum); _assert_(B_input);
	}

	/* Start looping on the number of nodes: */
	Gauss* gauss=this->NewGauss();
	for (int i=0;i<numnodes;i++){
		gauss->GaussNode(this->GetElementType(),i);

		eps_xx_input->GetInputValue(&eps_xx,gauss);
		eps_yy_input->GetInputValue(&eps_yy,gauss);
		eps_xy_input->GetInputValue(&eps_xy,gauss);
		if(dim==3){
			eps_xz_input->GetInputValue(&eps_xz,gauss);
			eps_yz_input->GetInputValue(&eps_yz,gauss);
			eps_zz_input->GetInputValue(&eps_zz,gauss);
		}
		else{eps_xz=0; eps_yz=0; eps_zz=0;}

		/* eps_eff^2 = exx^2 + eyy^2 + exy^2 + exz^2 + eyz^2 + exx*eyy */
		eps_eff=sqrt(eps_xx*eps_xx+eps_yy*eps_yy+eps_xy*eps_xy+eps_xz*eps_xz+eps_yz*eps_yz+eps_xx*eps_yy+epsmin*epsmin);

		B_input->GetInputValue(&B,gauss);
		n_input->GetInputValue(&n,gauss);
		damage_input->GetInputValue(&D,gauss);

		/* Compute threshold strain rate from threshold stress */
		eps_0=pow(sigma_0/B,n);

		if(eps_eff>eps_0){
			/* Compute damage on envelope curve for existing level of effective strain rate */
			envelopeD=1.-pow(eps_0/eps_eff,1./n)*exp(-(eps_eff-eps_0)/(eps_0*(kappa-1.)));

			if(envelopeD>D){
				newD[i]=envelopeD;
			}
			else newD[i]=D;
		}
		else newD[i]=D;
	}

	/* Add new damage input to DamageEnum and NewDamageEnum */
	this->AddInput(NewDamageEnum,newD,P1DGEnum);
	if(domaintype==Domain2DhorizontalEnum){
		this->AddInput(DamageDbarEnum,newD,this->GetElementType());
	}
	else{
		this->AddInput(DamageDEnum,newD,this->GetElementType());
	}

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(newD);
	delete gauss;

}/*}}}*/
void       Element::ComputeStrainRate(){/*{{{*/

	int         dim;
	IssmDouble *xyz_list = NULL;
	IssmDouble  epsilon[6];

	/*Retrieve all inputs we will be needing: */
	this->GetVerticesCoordinates(&xyz_list);
	parameters->FindParam(&dim,DomainDimensionEnum);
	Input* vx_input=this->GetInput(VxEnum); _assert_(vx_input);
	Input* vy_input=this->GetInput(VyEnum); _assert_(vy_input);
	Input* vz_input=NULL;
	if(dim==3){vz_input=this->GetInput(VzEnum); _assert_(vz_input);}

	/*Allocate arrays*/
	const int NUM_VERTICES = this->GetNumberOfVertices();
	_assert_(NUM_VERTICES<=MAXVERTICES);

	IssmDouble eps_xx[MAXVERTICES];
	IssmDouble eps_yy[MAXVERTICES];
	IssmDouble eps_zz[MAXVERTICES];
	IssmDouble eps_xy[MAXVERTICES];
	IssmDouble eps_xz[MAXVERTICES];
	IssmDouble eps_yz[MAXVERTICES];
	IssmDouble eps_ef[MAXVERTICES];

	/* Start looping on the number of vertices: */
	Gauss* gauss=this->NewGauss();
	for (int iv=0;iv<NUM_VERTICES;iv++){
		gauss->GaussVertex(iv);

		/*Compute strain rate viscosity and pressure: */
		if(dim==2)
		 this->StrainRateSSA(&epsilon[0],xyz_list,gauss,vx_input,vy_input);
		else
		 this->StrainRateFS(&epsilon[0],xyz_list,gauss,vx_input,vy_input,vz_input);

		if(dim==2){
			/* epsilon=[exx,eyy,exy];*/
			eps_xx[iv]=epsilon[0];
			eps_yy[iv]=epsilon[1];
			eps_xy[iv]=epsilon[2];
			/* eps_eff^2 = 1/2 ( exx^2 + eyy^2 + 2*exy^2 )*/
			eps_ef[iv] = 1./sqrt(2.)*sqrt(epsilon[0]*epsilon[0] + epsilon[1]*epsilon[1] + 2.*epsilon[2]*epsilon[2]);
		}
		else{
			/*epsilon=[exx eyy ezz exy exz eyz]*/
			eps_xx[iv]=epsilon[0];
			eps_yy[iv]=epsilon[1];
			eps_zz[iv]=epsilon[2];
			eps_xy[iv]=epsilon[3];
			eps_xz[iv]=epsilon[4];
			eps_yz[iv]=epsilon[5];
			/* eps_eff^2 = exx^2 + eyy^2 + exy^2 + exz^2 + eyz^2 + exx*eyy */
			eps_ef[iv] = sqrt(epsilon[0]*epsilon[0] + epsilon[1]*epsilon[1] + epsilon[3]*epsilon[3] +  epsilon[4]*epsilon[4] + epsilon[5]*epsilon[5] + epsilon[0]*epsilon[1]);
		}
	}

	/*Add Stress tensor components into inputs*/
	this->AddInput(StrainRatexxEnum,&eps_xx[0],P1Enum);
	this->AddInput(StrainRatexyEnum,&eps_xy[0],P1Enum);
	this->AddInput(StrainRatexzEnum,&eps_xz[0],P1Enum);
	this->AddInput(StrainRateyyEnum,&eps_yy[0],P1Enum);
	this->AddInput(StrainRateyzEnum,&eps_yz[0],P1Enum);
	this->AddInput(StrainRatezzEnum,&eps_zz[0],P1Enum);
	this->AddInput(StrainRateeffectiveEnum,&eps_ef[0],P1Enum);

	/*Clean up and return*/
	delete gauss;
	xDelete<IssmDouble>(xyz_list);
}
/*}}}*/
void       Element::CoordinateSystemTransform(IssmDouble** ptransform,Node** nodes_list,int numnodes,int* cs_array){/*{{{*/

	/*Some checks in debugging mode*/
	_assert_(numnodes && nodes_list);

	/*Get total number of dofs*/
	int numdofs = 0;
	for(int i=0;i<numnodes;i++){
		switch(cs_array[i]){
			case PressureEnum: numdofs+=1; break;
			case XYEnum:       numdofs+=2; break;
			case XYZEnum:      numdofs+=3; break;
			default: _error_("Coordinate system " << EnumToStringx(cs_array[i]) << " not supported yet");
		}
	}

	/*Allocate and initialize transform matrix*/
	IssmDouble* transform=xNew<IssmDouble>(numdofs*numdofs);
	for(int i=0;i<numdofs*numdofs;i++) transform[i]=0.0;

	/*Create transform matrix for all nodes (x,y for 2d and x,y,z for 3d). It is a block matrix
	 *for 3 nodes:

	 *     | T1 0  0 |
	 * Q = | 0  T2 0 |
	 *     | 0  0  T3|
	 *
	 * Where T1 is the transform matrix for node 1. It is a simple copy of the coordinate system
	 * associated to this node*/
	IssmDouble  coord_system[3][3];
	int         counter=0;
	for(int i=0;i<numnodes;i++){
		nodes_list[i]->GetCoordinateSystem(&coord_system[0][0]);
		switch(cs_array[i]){
			case PressureEnum:
				/*DO NOT change anything*/
				transform[(numdofs)*(counter) + counter] = 1.;
				counter+=1;
				break;
			case XYEnum:
				  {
				/*We remove the z component, we need to renormalize x and y: x=[x1 x2 0] y=[-x2 x1 0]*/
				IssmDouble norm = sqrt( coord_system[0][0]*coord_system[0][0] + coord_system[1][0]*coord_system[1][0]); _assert_(norm>1.e-4);
				transform[(numdofs)*(counter+0) + counter+0] =   coord_system[0][0]/norm;
				transform[(numdofs)*(counter+0) + counter+1] = - coord_system[1][0]/norm;
				transform[(numdofs)*(counter+1) + counter+0] =   coord_system[1][0]/norm;
				transform[(numdofs)*(counter+1) + counter+1] =   coord_system[0][0]/norm;
				counter+=2;
				  }
				break;
			case XYZEnum:
				/*The 3 coordinates are changed (x,y,z)*/
				transform[(numdofs)*(counter+0) + counter+0] = coord_system[0][0];
				transform[(numdofs)*(counter+0) + counter+1] = coord_system[0][1];
				transform[(numdofs)*(counter+0) + counter+2] = coord_system[0][2];
				transform[(numdofs)*(counter+1) + counter+0] = coord_system[1][0];
				transform[(numdofs)*(counter+1) + counter+1] = coord_system[1][1];
				transform[(numdofs)*(counter+1) + counter+2] = coord_system[1][2];
				transform[(numdofs)*(counter+2) + counter+0] = coord_system[2][0];
				transform[(numdofs)*(counter+2) + counter+1] = coord_system[2][1];
				transform[(numdofs)*(counter+2) + counter+2] = coord_system[2][2];
				counter+=3;
				break;
			default:
				_error_("Coordinate system " << EnumToStringx(cs_array[i]) << " not supported yet");
		}
	}

	/*Assign output pointer*/
	*ptransform=transform;
}
/*}}}*/
void       Element::DeepEcho(void){/*{{{*/

	_printf_(EnumToStringx(this->ObjectEnum())<<" element:\n");
	_printf_("   id : "<<this->id <<"\n");
	_printf_("   sid: "<<this->sid<<"\n");
	_printf_("   lid: "<<this->lid<<"\n");
	if(vertices){
		const int NUM_VERTICES = this->GetNumberOfVertices();
		for(int i=0;i<NUM_VERTICES;i++) vertices[i]->Echo();
	}
	else _printf_("vertices = NULL\n");

	if(nodes){
		int numnodes = this->GetNumberOfNodes();
		for(int i=0;i<numnodes;i++) nodes[i]->DeepEcho();
	}
	else _printf_("nodes = NULL\n");

	if (material) material->DeepEcho();
	else _printf_("material = NULL\n");

	_printf_("   parameters\n");
	if (parameters) parameters->DeepEcho();
	else _printf_("parameters = NULL\n");

	_printf_("   inputs\n");
	if(inputs) inputs->DeepEcho();
	else _printf_("inputs=NULL\n");

	return;
}
/*}}}*/
void       Element::DeleteMaterials(void){/*{{{*/
	delete this->material;
}/*}}}*/
void       Element::Delta18oParameterization(void){/*{{{*/

	/*Are we on the base? If not, return*/
	if(!IsOnBase()) return;

	const int NUM_VERTICES = this->GetNumberOfVertices();
	const int NUM_VERTICES_MONTHS_PER_YEAR	= NUM_VERTICES * 12;

	int        i;
	int*        vertexlids=xNew<int>(NUM_VERTICES);
	IssmDouble* monthlytemperatures=xNew<IssmDouble>(NUM_VERTICES_MONTHS_PER_YEAR);
	IssmDouble* monthlyprec=xNew<IssmDouble>(NUM_VERTICES_MONTHS_PER_YEAR);
	IssmDouble* TemperaturesPresentday=xNew<IssmDouble>(NUM_VERTICES_MONTHS_PER_YEAR);
	IssmDouble* TemperaturesLgm=xNew<IssmDouble>(NUM_VERTICES_MONTHS_PER_YEAR);
	IssmDouble* PrecipitationsPresentday=xNew<IssmDouble>(NUM_VERTICES_MONTHS_PER_YEAR);
	IssmDouble* tmp=xNew<IssmDouble>(NUM_VERTICES);

	/*Recover parameters*/
	IssmDouble time,yts,finaltime;
	this->parameters->FindParam(&time,TimeEnum);
	this->parameters->FindParam(&yts,ConstantsYtsEnum);
	this->parameters->FindParam(&finaltime,TimesteppingFinalTimeEnum);
	this->GetVerticesLidList(vertexlids);
	IssmDouble time_yr=floor(time/yts)*yts;

	/*Recover present day temperature and precipitation*/
	DatasetInput* dinput1=this->GetDatasetInput(SmbTemperaturesPresentdayEnum);   _assert_(dinput1);
	DatasetInput* dinput2=this->GetDatasetInput(SmbTemperaturesLgmEnum);          _assert_(dinput2);
	DatasetInput* dinput3=this->GetDatasetInput(SmbPrecipitationsPresentdayEnum); _assert_(dinput3);

	/*loop over vertices: */
	Gauss* gauss=this->NewGauss();
	for(int month=0;month<12;month++){
		for(int iv=0;iv<NUM_VERTICES;iv++){
			gauss->GaussVertex(iv);
			dinput1->GetInputValue(&TemperaturesPresentday[iv*12+month],gauss,month);
			dinput2->GetInputValue(&TemperaturesLgm[iv*12+month],gauss,month);
			dinput3->GetInputValue(&PrecipitationsPresentday[iv*12+month],gauss,month);

			PrecipitationsPresentday[iv*12+month]=PrecipitationsPresentday[iv*12+month]*yts;
		}
	}

	/*Recover delta18o and Delta18oSurface at present day, lgm and at time t*/
	IssmDouble Delta18oPresent,Delta18oLgm,Delta18oTime;
	IssmDouble Delta18oSurfacePresent,Delta18oSurfaceLgm,Delta18oSurfaceTime;
	this->parameters->FindParam(&Delta18oPresent,SmbDelta18oEnum,finaltime);
	this->parameters->FindParam(&Delta18oLgm,SmbDelta18oEnum,(finaltime-(21000*yts)));
	this->parameters->FindParam(&Delta18oTime,SmbDelta18oEnum,time);
	this->parameters->FindParam(&Delta18oSurfacePresent,SmbDelta18oSurfaceEnum,finaltime);
	this->parameters->FindParam(&Delta18oSurfaceLgm,SmbDelta18oSurfaceEnum,(finaltime-(21000*yts)));
	this->parameters->FindParam(&Delta18oSurfaceTime,SmbDelta18oSurfaceEnum,time);

	/*Compute the temperature and precipitation*/
	for(int iv=0;iv<NUM_VERTICES;iv++){
		ComputeDelta18oTemperaturePrecipitation(Delta18oSurfacePresent, Delta18oSurfaceLgm, Delta18oSurfaceTime,
					Delta18oPresent, Delta18oLgm, Delta18oTime,
					&PrecipitationsPresentday[iv*12],
					&TemperaturesLgm[iv*12], &TemperaturesPresentday[iv*12],
					&monthlytemperatures[iv*12], &monthlyprec[iv*12]);
	}

	/*Update inputs*/
	for (int imonth=0;imonth<12;imonth++) {
		for(i=0;i<NUM_VERTICES;i++) tmp[i]=monthlytemperatures[i*12+imonth];
		switch(this->ObjectEnum()){
			case TriaEnum: this->inputs->SetTriaDatasetInput(SmbMonthlytemperaturesEnum,imonth,P1Enum,NUM_VERTICES,vertexlids,tmp); break;
			case PentaEnum: this->inputs->SetPentaDatasetInput(SmbMonthlytemperaturesEnum,imonth,P1Enum,NUM_VERTICES,vertexlids,tmp); break;
			default: _error_("Not implemented yet");
		}
		for(i=0;i<NUM_VERTICES;i++) tmp[i]=monthlyprec[i*12+imonth]/yts;
		switch(this->ObjectEnum()){
			case TriaEnum: this->inputs->SetTriaDatasetInput(SmbPrecipitationEnum,imonth,P1Enum,NUM_VERTICES,vertexlids,tmp); break;
			case PentaEnum: this->inputs->SetPentaDatasetInput(SmbPrecipitationEnum,imonth,P1Enum,NUM_VERTICES,vertexlids,tmp); break;
			default: _error_("Not implemented yet");
		}
	}

	switch(this->ObjectEnum()){
		case TriaEnum: break;
		case PentaEnum:
		case TetraEnum:
         this->DatasetInputExtrude(SmbMonthlytemperaturesEnum,-1);
         this->DatasetInputExtrude(SmbPrecipitationEnum,-1);
         break;
		default: _error_("Not implemented yet");
	}

	/*clean-up*/
	delete gauss;
	xDelete<IssmDouble>(monthlytemperatures);
	xDelete<IssmDouble>(monthlyprec);
	xDelete<IssmDouble>(TemperaturesPresentday);
	xDelete<IssmDouble>(TemperaturesLgm);
	xDelete<IssmDouble>(PrecipitationsPresentday);
	xDelete<IssmDouble>(tmp);
	xDelete<int>(vertexlids);

} /*}}}*/
void       Element::Delta18opdParameterization(void){/*{{{*/
	/*Are we on the base? If not, return*/
	if(!IsOnBase()) return;

	const int NUM_VERTICES 					= this->GetNumberOfVertices();
	const int NUM_VERTICES_MONTHS_PER_YEAR	= NUM_VERTICES * 12;

	int        	i,offset;
	int*        vertexlids=xNew<int>(NUM_VERTICES);
	IssmDouble* monthlytemperatures=xNew<IssmDouble>(NUM_VERTICES_MONTHS_PER_YEAR);
	IssmDouble* monthlyprec=xNew<IssmDouble>(NUM_VERTICES_MONTHS_PER_YEAR);
	IssmDouble* TemperaturesPresentday=xNew<IssmDouble>(NUM_VERTICES_MONTHS_PER_YEAR);
	IssmDouble* PrecipitationsPresentday=xNew<IssmDouble>(NUM_VERTICES_MONTHS_PER_YEAR);
	IssmDouble* TemperaturesReconstructed=xNew<IssmDouble>(NUM_VERTICES_MONTHS_PER_YEAR);
	IssmDouble* PrecipitationsReconstructed=xNew<IssmDouble>(NUM_VERTICES_MONTHS_PER_YEAR);
	IssmDouble* tmp=xNew<IssmDouble>(NUM_VERTICES);
	IssmDouble Delta18oTime;
	IssmDouble f;
	IssmDouble time,yts,time_yr,month,time_climt,time_climp,del_clim;
	this->parameters->FindParam(&time,TimeEnum);
	this->parameters->FindParam(&yts,ConstantsYtsEnum);
	this->parameters->FindParam(&f,SmbFEnum);
	this->GetVerticesLidList(vertexlids);
	time_yr=floor(time/yts)*yts;
	time_climt=ceil(time/yts + 1e-10)*yts;
	time_climp=ceil(time/yts + 1e-10)*yts;

	/*Get some pdd parameters*/
	bool isTemperatureScaled,isPrecipScaled;
	IssmDouble dpermil=this->FindParam(SmbDpermilEnum);
	this->parameters->FindParam(&isTemperatureScaled,SmbIstemperaturescaledEnum);
	this->parameters->FindParam(&isPrecipScaled,SmbIsprecipscaledEnum);

	/*Recover present day temperature and precipitation*/
	DatasetInput *dinput3 = NULL;
	DatasetInput *dinput4 = NULL;
	int            offset_t,offset_p,N;
	if(!isTemperatureScaled){
		IssmDouble* time_temp_scaled = NULL;
		parameters->FindParam(&time_temp_scaled,&N,SmbTemperaturesReconstructedYearsEnum);
		if(!binary_search(&offset_t,time_climt,time_temp_scaled,N)) _error_("time not sorted?");
		if(offset_t<0) offset_t=0;
		xDelete<IssmDouble>(time_temp_scaled);
		dinput3=this->GetDatasetInput(SmbTemperaturesReconstructedEnum); _assert_(dinput3);
	}
	if(!isPrecipScaled){
		IssmDouble* time_precip_scaled = NULL;
		parameters->FindParam(&time_precip_scaled,&N,SmbPrecipitationsReconstructedYearsEnum);
		if(!binary_search(&offset_p,time_climt,time_precip_scaled,N)) _error_("time not sorted?");
		if(offset_p<0) offset_p=0;
		xDelete<IssmDouble>(time_precip_scaled);
		dinput4=this->GetDatasetInput(SmbPrecipitationsReconstructedEnum); _assert_(dinput4);
	}

	/*Get present day temp and precip (monthly)*/
	DatasetInput *dinput1 = this->GetDatasetInput(SmbTemperaturesPresentdayEnum);   _assert_(dinput1);
	DatasetInput *dinput2 = this->GetDatasetInput(SmbPrecipitationsPresentdayEnum); _assert_(dinput2);

	/*loop over vertices: */
	Gauss* gauss=this->NewGauss();
	for(int month=0;month<12;month++) {
		for(int iv=0;iv<NUM_VERTICES;iv++) {
			gauss->GaussVertex(iv);
			dinput1->GetInputValue(&TemperaturesPresentday[iv*12+month],gauss,month);
			dinput2->GetInputValue(&PrecipitationsPresentday[iv*12+month],gauss,month);
			PrecipitationsPresentday[iv*12+month]=PrecipitationsPresentday[iv*12+month]*yts;

			if(!isTemperatureScaled){
				dinput3->GetInputValue(&TemperaturesReconstructed[iv*12+month],gauss,offset_t*12+month);
			}
			if(!isPrecipScaled){
				dinput4->GetInputValue(&PrecipitationsReconstructed[iv*12+month],gauss,offset_p*12+month);
				PrecipitationsReconstructed[iv*12+month]=PrecipitationsReconstructed[iv*12+month]*yts;
			}
		}
	}

	/*Recover interpolation parameters at time t*/
	this->parameters->FindParam(&Delta18oTime,SmbDelta18oEnum,time);

	/*Compute the temperature and precipitation*/
	for(int iv=0;iv<NUM_VERTICES;iv++){
		ComputeD18OTemperaturePrecipitationFromPD(Delta18oTime,dpermil,isTemperatureScaled,isPrecipScaled,
					f,&PrecipitationsPresentday[iv*12], &TemperaturesPresentday[iv*12],
					&PrecipitationsReconstructed[iv*12], &TemperaturesReconstructed[iv*12],
					&monthlytemperatures[iv*12], &monthlyprec[iv*12]);
	}

	/*Update inputs*/
	for (int imonth=0;imonth<12;imonth++) {
		for(i=0;i<NUM_VERTICES;i++) tmp[i]=monthlytemperatures[i*12+imonth];
		switch(this->ObjectEnum()){
			case TriaEnum:  this->inputs->SetTriaDatasetInput(SmbMonthlytemperaturesEnum,imonth,P1Enum,NUM_VERTICES,vertexlids,tmp); break;
			case PentaEnum: this->inputs->SetPentaDatasetInput(SmbMonthlytemperaturesEnum,imonth,P1Enum,NUM_VERTICES,vertexlids,tmp); break;
			default: _error_("Not implemented yet");
		}
		for(i=0;i<NUM_VERTICES;i++) tmp[i]=monthlyprec[i*12+imonth]/yts;
		switch(this->ObjectEnum()){
			case TriaEnum:  this->inputs->SetTriaDatasetInput(SmbPrecipitationEnum,imonth,P1Enum,NUM_VERTICES,vertexlids,tmp); break;
			case PentaEnum: this->inputs->SetPentaDatasetInput(SmbPrecipitationEnum,imonth,P1Enum,NUM_VERTICES,vertexlids,tmp); break;
			default: _error_("Not implemented yet");
		}
	}

	switch(this->ObjectEnum()){
		case TriaEnum: break;
		case PentaEnum:
		case TetraEnum:
         this->DatasetInputExtrude(SmbMonthlytemperaturesEnum,-1);
         this->DatasetInputExtrude(SmbPrecipitationEnum,-1);
         break;
		default: _error_("Not implemented yet");
	}

	/*clean-up*/
	delete gauss;
	xDelete<IssmDouble>(monthlytemperatures);
	xDelete<IssmDouble>(monthlyprec);
	xDelete<IssmDouble>(TemperaturesPresentday);
	xDelete<IssmDouble>(PrecipitationsPresentday);
	xDelete<IssmDouble>(TemperaturesReconstructed);
	xDelete<IssmDouble>(PrecipitationsReconstructed);
	xDelete<IssmDouble>(tmp);
	xDelete<int>(vertexlids);

} /*}}}*/
void       Element::SmbGradCompParameterization(void){/*{{{*/

	/*Are we on the base? If not, return*/
	if(!IsOnBase()) return;
	const int NUM_VERTICES = this->GetNumberOfVertices();

	IssmDouble accuref, runoffref;  //reference values at given altitude
	IssmDouble accugrad, runoffgrad; //gradients from reference altitude

	IssmDouble smb[MAXVERTICES];
	IssmDouble surf[MAXVERTICES];
	IssmDouble accu[MAXVERTICES];
	IssmDouble runoff[MAXVERTICES];

	/*Get material parameters :*/
	IssmDouble rho_water = this->FindParam(MaterialsRhoFreshwaterEnum);
	IssmDouble rho_ice   = this->FindParam(MaterialsRhoIceEnum);

	/*Recover parameters*/
	IssmDouble time       = this->FindParam(TimeEnum);
	IssmDouble accualti   = this->FindParam(SmbAccualtiEnum);
	IssmDouble runoffalti = this->FindParam(SmbRunoffaltiEnum);

	/*Recover reference values at current time*/
	parameters->FindParam(&accugrad,SmbAccugradEnum,time);
	parameters->FindParam(&runoffgrad,SmbRunoffgradEnum,time);
	parameters->FindParam(&accuref,SmbAccurefEnum,time);
	parameters->FindParam(&runoffref,SmbRunoffrefEnum,time);

	/*Recover surface elevation*/
	GetInputListOnVertices(&surf[0],SurfaceEnum);

	/*Compute the temperature and precipitation*/
	for(int iv=0;iv<NUM_VERTICES;iv++){
		accu[iv]   = max(0.,(accuref+(surf[iv]-accualti)*accugrad));
		runoff[iv] = max(0.,(runoffref+(surf[iv]-runoffalti)*runoffgrad));
		smb[iv]    = (accu[iv]-runoff[iv])*rho_ice/rho_water;
	}

	switch(this->ObjectEnum()){
	case TriaEnum:
		this->AddInput(SmbMassBalanceSubstepEnum,&smb[0],P1Enum);
		this->AddInput(SmbRunoffSubstepEnum,&runoff[0],P1Enum);
		break;
	case PentaEnum:
		this->AddInput(SmbMassBalanceSubstepEnum,&smb[0],P1Enum);
		this->AddInput(SmbRunoffSubstepEnum,&runoff[0],P1Enum);
		this->InputExtrude(SmbMassBalanceSubstepEnum,-1);
		this->InputExtrude(SmbRunoffSubstepEnum,-1);
		break;
	default: _error_("Not implemented yet");
	}
}
/*}}}*/
IssmDouble Element::Divergence(void){/*{{{*/
	/*Compute element divergence*/

	/*Intermediaries*/
	int        dim;
	IssmDouble Jdet;
	IssmDouble divergence=0.;
	IssmDouble dvx[3],dvy[3],dvz[3];
	IssmDouble *xyz_list = NULL;

	/*Get inputs and parameters*/
	this->FindParam(&dim,DomainDimensionEnum);
	Input* vx_input = this->GetInput(VxEnum); _assert_(vx_input);
	Input* vy_input = this->GetInput(VyEnum); _assert_(vy_input);
	Input* vz_input = NULL;
	if(dim==3){
		vz_input = this->GetInput(VzEnum); _assert_(vz_input);
	}
	this->GetVerticesCoordinates(&xyz_list);

	Gauss* gauss=this->NewGauss(5);
	while(gauss->next()){

		this->JacobianDeterminant(&Jdet,xyz_list,gauss);

		/*Get strain rate assuming that epsilon has been allocated*/
		vx_input->GetInputDerivativeValue(&dvx[0],xyz_list,gauss);
		vy_input->GetInputDerivativeValue(&dvy[0],xyz_list,gauss);
		if(dim==2){
			divergence += (dvx[0]+dvy[1])*gauss->weight*Jdet;
		}
		else{
			vz_input->GetInputDerivativeValue(&dvz[0],xyz_list,gauss);
			divergence += (dvx[0]+dvy[1]+dvz[2])*gauss->weight*Jdet;
		}

	}

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	delete gauss;
	return divergence;
}/*}}}*/
void       Element::dViscositydBFS(IssmDouble* pdmudB,int dim,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input,Input* vz_input){/*{{{*/

	/*Intermediaries*/
	int materialstype;
	IssmDouble dmudB;
	IssmDouble epsilon3d[6];/* epsilon=[exx,eyy,exy,exy,exz,eyz];    */
	IssmDouble epsilon2d[3];/* epsilon=[exx,eyy,exy];    */
	IssmDouble eps_eff;
	IssmDouble eps0=1.e-27;

	if(dim==3){
		/* eps_eff^2 = exx^2 + eyy^2 + exy^2 + exz^2 + eyz^2 + exx*eyy */
		this->StrainRateFS(&epsilon3d[0],xyz_list,gauss,vx_input,vy_input,vz_input);
		eps_eff = sqrt(epsilon3d[0]*epsilon3d[0] + epsilon3d[1]*epsilon3d[1] + epsilon3d[3]*epsilon3d[3] +  epsilon3d[4]*epsilon3d[4] + epsilon3d[5]*epsilon3d[5] + epsilon3d[0]*epsilon3d[1]+eps0*eps0);
	}
	else{
		/* eps_eff^2 = 1/2 ( exx^2 + eyy^2 + 2*exy^2 )*/
		this->StrainRateSSA(&epsilon2d[0],xyz_list,gauss,vx_input,vy_input);
		eps_eff = 1./sqrt(2.)*sqrt(epsilon2d[0]*epsilon2d[0] + epsilon2d[1]*epsilon2d[1] + 2.*epsilon2d[2]*epsilon2d[2]);
	}
	/*Get viscosity*/
	materialstype=this->material->ObjectEnum();
	switch(materialstype){
		case MaticeEnum:
			material->GetViscosity_B(&dmudB,eps_eff,gauss);
			break;
		case MatestarEnum:
			material->ViscosityBFS(&dmudB,dim,xyz_list,gauss,vx_input,vy_input,vz_input,eps_eff);
			break;
		default: _error_("not supported");
	}

	/*Assign output pointer*/
	*pdmudB=dmudB;

}
/*}}}*/
void       Element::dViscositydBHO(IssmDouble* pdmudB,int dim,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input){/*{{{*/

	/*Intermediaries*/
	int materialstype;
	IssmDouble dmudB;
	IssmDouble epsilon3d[5];/* epsilon=[exx,eyy,exy,exy,exz,eyz];    */
	IssmDouble epsilon2d[2];/* epsilon=[exx,eyy,exy];    */
	IssmDouble eps_eff;
	IssmDouble eps0=1.e-27;

	if(dim==3){
		/* eps_eff^2 = exx^2 + eyy^2 + exy^2 + exz^2 + eyz^2 + exx*eyy */
		this->StrainRateHO(&epsilon3d[0],xyz_list,gauss,vx_input,vy_input);
		eps_eff = sqrt(epsilon3d[0]*epsilon3d[0] + epsilon3d[1]*epsilon3d[1] + epsilon3d[2]*epsilon3d[2] + epsilon3d[3]*epsilon3d[3] +  epsilon3d[4]*epsilon3d[4] + epsilon3d[0]*epsilon3d[1]+eps0*eps0);
	}
	else{
		/* eps_eff^2 = 1/2 ( exx^2 + eyy^2 + 2*exy^2 )*/
		this->StrainRateHO2dvertical(&epsilon2d[0],xyz_list,gauss,vx_input,vy_input);
		eps_eff = 1./sqrt(2.)*sqrt(epsilon2d[0]*epsilon2d[0] + 2.*epsilon2d[1]*epsilon2d[1] + eps0*eps0);
	}
	/*Get viscosity*/
	materialstype=this->material->ObjectEnum();
	switch(materialstype){
		case MaticeEnum:
			material->GetViscosity_B(&dmudB,eps_eff,gauss);
			break;
		case MatestarEnum:
			material->ViscosityBHO(&dmudB,dim,xyz_list,gauss,vx_input,vy_input,eps_eff);
			break;
		default: _error_("not supported");
	}

	/*Assign output pointer*/
	*pdmudB=dmudB;

}
/*}}}*/
void       Element::dViscositydBMOLHO(IssmDouble* pdmudB,int dim,IssmDouble* xyz_list,Gauss* gauss,Input* vxbase_input,Input* vybase_input, Input* vxshear_input ,Input* vyshear_input,Input* thickness_input,Input* n_input, IssmDouble zeta){/*{{{*/

	/*Intermediaries*/
	int materialstype;
	IssmDouble dmudB;
	IssmDouble epsilon[5];/* epsilon=[exx,eyy,exy,exy,exz,eyz];    */
	IssmDouble eps_eff;
	IssmDouble eps0=1.e-27;

	/* eps_eff^2 = exx^2 + eyy^2 + exy^2 + exz^2 + eyz^2 + exx*eyy */
   this->StrainRateMOLHO(&epsilon[0],xyz_list,gauss,
				vxbase_input,vybase_input,vxshear_input,vyshear_input,thickness_input,n_input,zeta);
   eps_eff=sqrt(epsilon[0]*epsilon[0] + epsilon[1]*epsilon[1] + epsilon[2]*epsilon[2]
				+ epsilon[3]*epsilon[3] + epsilon[4]*epsilon[4] + epsilon[0]*epsilon[1] + eps0*eps0);

	/*Get viscosity*/
	materialstype=this->material->ObjectEnum();
	switch(materialstype){
		case MaticeEnum:
			material->GetViscosity_B(&dmudB,eps_eff,gauss);
			break;
		default: _error_("not supported");
	}

	/*Assign output pointer*/
	*pdmudB=dmudB;

}
/*}}}*/
void       Element::dViscositydBSSA(IssmDouble* pdmudB,int dim,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input){/*{{{*/

	/*Intermediaries*/
	int materialstype;
	IssmDouble dmudB;
	IssmDouble epsilon2d[3];/* epsilon=[exx,eyy,exy];    */
	IssmDouble epsilon1d;   /* epsilon=[exx];    */
	IssmDouble eps_eff;

	if(dim==2){
		/* eps_eff^2 = exx^2 + eyy^2 + exy^2 + exx*eyy*/
		this->StrainRateSSA(&epsilon2d[0],xyz_list,gauss,vx_input,vy_input);
		eps_eff = sqrt(epsilon2d[0]*epsilon2d[0] + epsilon2d[1]*epsilon2d[1] + epsilon2d[2]*epsilon2d[2] + epsilon2d[0]*epsilon2d[1]);
	}
	else{
		/* eps_eff^2 = 1/2 exx^2*/
		this->StrainRateSSA1d(&epsilon1d,xyz_list,gauss,vx_input);
		eps_eff = sqrt(epsilon1d*epsilon1d/2.);
	}

	/*Get viscosity*/
	materialstype=this->material->ObjectEnum();
	switch(materialstype){
		case MaticeEnum:
			material->GetViscosity_B(&dmudB,eps_eff,gauss);
			break;
		case MatestarEnum:
			material->ViscosityBSSA(&dmudB,dim,xyz_list,gauss,vx_input,vy_input,eps_eff);
			break;
		default: _error_("not supported");
	}

	/*Assign output pointer*/
	*pdmudB=dmudB;

}
/*}}}*/
void       Element::dViscositydDSSA(IssmDouble* pdmudB,int dim,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input){/*{{{*/

	/*Intermediaries*/
	IssmDouble dmudB;
	IssmDouble epsilon2d[3];/* epsilon=[exx,eyy,exy];    */
	IssmDouble epsilon1d;   /* epsilon=[exx];    */
	IssmDouble eps_eff;

	if(dim==2){
		/* eps_eff^2 = exx^2 + eyy^2 + exy^2 + exx*eyy*/
		this->StrainRateSSA(&epsilon2d[0],xyz_list,gauss,vx_input,vy_input);
		eps_eff = sqrt(epsilon2d[0]*epsilon2d[0] + epsilon2d[1]*epsilon2d[1] + epsilon2d[2]*epsilon2d[2] + epsilon2d[0]*epsilon2d[1]);
	}
	else{
		/* eps_eff^2 = 1/2 exx^2*/
		this->StrainRateSSA1d(&epsilon1d,xyz_list,gauss,vx_input);
		eps_eff = sqrt(epsilon1d*epsilon1d/2.);
	}

	/*Get viscosity*/
	material->GetViscosity_D(&dmudB,eps_eff,gauss);

	/*Assign output pointer*/
	*pdmudB=dmudB;

}
/*}}}*/
void       Element::Echo(void){/*{{{*/
	_printf_(EnumToStringx(this->ObjectEnum())<<" element:\n");
	_printf_("   id : "<<this->id <<"\n");
	_printf_("   sid: "<<this->sid<<"\n");
	_printf_("   lid: "<<this->lid<<"\n");
	if(vertices){
		const int NUM_VERTICES = this->GetNumberOfVertices();
		for(int i=0;i<NUM_VERTICES;i++) vertices[i]->Echo();
	}
	else _printf_("vertices = NULL\n");

	if(nodes){
		int numnodes = this->GetNumberOfNodes();
		for(int i=0;i<numnodes;i++) {
			_printf_("nodes[" << i << "] = " << nodes[i]);
			nodes[i]->Echo();
		}
	}
	else _printf_("nodes = NULL\n");

	if (material) material->Echo();
	else _printf_("material = NULL\n");

	_printf_("   parameters\n");
	if (parameters) parameters->Echo();
	else _printf_("parameters = NULL\n");

	_printf_("   inputs\n");
	if (inputs) inputs->Echo();
	else _printf_("inputs=NULL\n");
}
/*}}}*/
void       Element::FindParam(bool* pvalue,int paramenum){/*{{{*/
	this->parameters->FindParam(pvalue,paramenum);
}/*}}}*/
void       Element::FindParam(int* pvalue,int paramenum){/*{{{*/
	this->parameters->FindParam(pvalue,paramenum);
}/*}}}*/
void       Element::FindParam(IssmDouble* pvalue,int paramenum){/*{{{*/
	this->parameters->FindParam(pvalue,paramenum);
}/*}}}*/
IssmDouble Element::FindParam(int paramenum){/*{{{*/
	return this->parameters->FindParam(paramenum);
}/*}}}*/
void       Element::FindParam(int** pvalues,int* psize,int paramenum){/*{{{*/
	this->parameters->FindParam(pvalues,psize,paramenum);
}/*}}}*/
IssmDouble Element::FloatingArea(IssmDouble* mask, bool scaled){/*{{{*/

	/*Retrieve values of the mask defining the element: */
	for(int i=0;i<this->GetNumberOfVertices();i++){
		if(mask[this->vertices[i]->Sid()]<=0.){
			return 0.;
		}
	}

	/*Return: */
	return this->FloatingArea(scaled);
}
/*}}}*/
void       Element::FrictionAlpha2CreateInput(void){/*{{{*/

	/*Return if element is inactive*/
	if(this->IsAllFloating() || !this->IsIceInElement()) return;

	/*Intermediaries*/
	int      domaintype, dim;
	Element* basalelement = NULL;

	/*Get basal element*/
	this->FindParam(&domaintype,DomainTypeEnum);
	switch(domaintype){
		case Domain2DhorizontalEnum:
			basalelement = this;
			dim = 2;
			break;
		case Domain3DEnum: case Domain2DverticalEnum:
			if(!this->IsOnBase()) return;
			basalelement = this->SpawnBasalElement(true);
			dim = 2;
			break;
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}

	int numvertices = basalelement->GetNumberOfVertices();

	/*build friction object, used later on: */
	Friction* friction=new Friction(basalelement,dim);
	IssmDouble alpha2_list[MAXVERTICES];
	IssmDouble alpha2;

	Gauss* gauss=basalelement->NewGauss();
	for(int i=0;i<numvertices;i++){
		gauss->GaussVertex(i);

		friction->GetAlpha2(&alpha2,gauss);
		alpha2_list[i] = alpha2;
	}

	/*Clean up and return*/
	delete gauss;
	delete friction;
	if(basalelement->IsSpawnedElement()){basalelement->DeleteMaterials(); delete basalelement;};
	this->AddBasalInput(FrictionAlpha2Enum,&alpha2_list[0],P1Enum);
}
/*}}}*/
void       Element::GetDofList(int** pdoflist,int approximation_enum,int setenum,bool hideclones){/*{{{*/

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = this->GetNumberOfNodes();

	/*First, figure out size of doflist and create it: */
	int numberofdofs=0;
	for(int i=0;i<numnodes;i++) numberofdofs+=nodes[i]->GetNumberOfDofs(approximation_enum,GsetEnum);

	/*Allocate output*/
	int* doflist=xNew<int>(numberofdofs);

	/*Populate: */
	int count=0;
	for(int i=0;i<numnodes;i++){
		nodes[i]->GetDofList(&doflist[count],approximation_enum,setenum,hideclones);
		count+=nodes[i]->GetNumberOfDofs(approximation_enum,GsetEnum);
	}

	/*Assign output pointers:*/
	*pdoflist=doflist;
}
/*}}}*/
void       Element::GetDofListLocal(int** pdoflist,int approximation_enum,int setenum){/*{{{*/

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = this->GetNumberOfNodes();

	/*First, figure out size of doflist and create it: */
	int numberofdofs=0;
	for(int i=0;i<numnodes;i++) numberofdofs+=nodes[i]->GetNumberOfDofs(NoneApproximationEnum,GsetEnum);

	/*Allocate output*/
	int* doflist=xNew<int>(numberofdofs);

	/*Populate: */
	int count=0;
	for(int i=0;i<numnodes;i++){
		_assert_(count<numberofdofs);
		nodes[i]->GetDofListLocal(doflist+count,approximation_enum,setenum);
		count+=nodes[i]->GetNumberOfDofs(approximation_enum,GsetEnum);
	}

	/*Assign output pointers:*/
	*pdoflist=doflist;
}
/*}}}*/
void       Element::GetDofListPressure(int** pdoflist,int setenum){/*{{{*/

	/*Fetch number of nodes and dof for this finite element*/
	int vnumnodes = this->NumberofNodesVelocity();
	int pnumnodes = this->NumberofNodesPressure();

	/*First, figure out size of doflist and create it: */
	int numberofdofs=0;
	for(int i=vnumnodes;i<vnumnodes+pnumnodes;i++) numberofdofs+=nodes[i]->GetNumberOfDofs(FSApproximationEnum,setenum);

	/*Allocate output*/
	int* doflist=xNew<int>(numberofdofs);

	/*Populate: */
	int count=0;
	for(int i=vnumnodes;i<vnumnodes+pnumnodes;i++){
		nodes[i]->GetDofList(doflist+count,FSApproximationEnum,setenum);
		count+=nodes[i]->GetNumberOfDofs(FSApproximationEnum,setenum);
	}

	/*Assign output pointers:*/
	*pdoflist=doflist;
}
/*}}}*/
void       Element::GetDofListVelocity(int** pdoflist,int setenum){/*{{{*/

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = this->NumberofNodesVelocity();

	/*First, figure out size of doflist and create it: */
	int numberofdofs=0;
	for(int i=0;i<numnodes;i++) numberofdofs+=nodes[i]->GetNumberOfDofs(FSvelocityEnum,setenum);

	/*Allocate output*/
	int* doflist=xNew<int>(numberofdofs);

	/*Populate: */
	int count=0;
	for(int i=0;i<numnodes;i++){
		nodes[i]->GetDofList(doflist+count,FSvelocityEnum,setenum);
		count+=nodes[i]->GetNumberOfDofs(FSvelocityEnum,setenum);
	}

	/*Assign output pointers:*/
	*pdoflist=doflist;
}
/*}}}*/
void       Element::GetDofListLocalPressure(int** pdoflist,int setenum){/*{{{*/

	/*Fetch number of nodes and dof for this finite element*/
	int vnumnodes = this->NumberofNodesVelocity();
	int pnumnodes = this->NumberofNodesPressure();

	/*First, figure out size of doflist and create it: */
	int numberofdofs=0;
	for(int i=vnumnodes;i<vnumnodes+pnumnodes;i++) numberofdofs+=nodes[i]->GetNumberOfDofs(FSApproximationEnum,setenum);

	/*Allocate output*/
	int* doflist=xNew<int>(numberofdofs);

	/*Populate: */
	int count=0;
	for(int i=vnumnodes;i<vnumnodes+pnumnodes;i++){
		nodes[i]->GetDofListLocal(doflist+count,FSApproximationEnum,setenum);
		count+=nodes[i]->GetNumberOfDofs(FSApproximationEnum,setenum);
	}

	/*Assign output pointers:*/
	*pdoflist=doflist;
}
/*}}}*/
void       Element::GetDofListLocalVelocity(int** pdoflist,int setenum){/*{{{*/

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = this->NumberofNodesVelocity();

	/*First, figure out size of doflist and create it: */
	int numberofdofs=0;
	for(int i=0;i<numnodes;i++) numberofdofs+=nodes[i]->GetNumberOfDofs(FSvelocityEnum,setenum);

	/*Allocate output*/
	int* doflist=xNew<int>(numberofdofs);

	/*Populate: */
	int count=0;
	for(int i=0;i<numnodes;i++){
		nodes[i]->GetDofListLocal(doflist+count,FSvelocityEnum,setenum);
		count+=nodes[i]->GetNumberOfDofs(FSvelocityEnum,setenum);
	}

	/*Assign output pointers:*/
	*pdoflist=doflist;
}
/*}}}*/
void       Element::GetInputListOnNodes(IssmDouble* pvalue,int enumtype,IssmDouble defaultvalue){/*{{{*/
	Input *input    = this->GetInput(enumtype);
	this->GetInputListOnNodes(pvalue,input,defaultvalue);
}
/*}}}*/
void       Element::GetInputListOnNodes(IssmDouble* pvalue,int enumtype){/*{{{*/

	Input *input    = this->GetInput(enumtype);
	if(!input) _error_("Input " << EnumToStringx(enumtype) << " not found in element");
	this->GetInputListOnNodes(pvalue,input,0.);

}
/*}}}*/
void       Element::GetInputListOnNodesVelocity(IssmDouble* pvalue,int enumtype){/*{{{*/

	_assert_(pvalue);

	int     numnodes = this->NumberofNodesVelocity();
	Input *input    = this->GetInput(enumtype);
	if(!input) _error_("Input " << EnumToStringx(enumtype) << " not found in element");

	/* Start looping on the number of vertices: */
	Gauss* gauss=this->NewGauss();
	for(int iv=0;iv<numnodes;iv++){
		gauss->GaussNode(this->VelocityInterpolation(),iv);
		input->GetInputValue(&pvalue[iv],gauss);
	}
	delete gauss;
}
/*}}}*/
void       Element::GetInputListOnVertices(IssmDouble* pvalue,int enumtype){/*{{{*/

	/*Recover input*/
	Input* input2=this->GetInput(enumtype);
	if(!input2) _error_("input "<<EnumToStringx(enumtype)<<" not found in element");
	this->GetInputListOnVertices(pvalue,input2,0.);
}
/*}}}*/
void       Element::GetInputListOnVerticesAtTime(IssmDouble* pvalue, int enumtype, IssmDouble time){/*{{{*/

	/*Recover input*/
	Input* input=this->GetInput(enumtype,time);
	if (!input) _error_("Input " << EnumToStringx(enumtype) << " not found in element");
	this->GetInputListOnVertices(pvalue,input,0.);
}
/*}}}*/
void       Element::GetInputListOnVertices(IssmDouble* pvalue,int enumtype,IssmDouble defaultvalue){/*{{{*/
	Input* input=this->GetInput(enumtype);
	this->GetInputListOnVertices(pvalue,input,defaultvalue);
}
/*}}}*/
void       Element::GetInputLocalMinMaxOnNodes(IssmDouble* min,IssmDouble* max,IssmDouble* ug){/*{{{*/

	/*Get number of nodes for this element*/
	int numnodes = this->GetNumberOfNodes();

	/*Some checks to avoid segmentation faults*/
	_assert_(ug);
	_assert_(numnodes>0);
	_assert_(nodes);

	/*Get element minimum/maximum*/
	IssmDouble input_min = ug[nodes[0]->GetDof(0,GsetEnum)];
	IssmDouble input_max = input_min;
	for(int i=1;i<numnodes;i++){
		if(ug[nodes[i]->GetDof(0,GsetEnum)] < input_min) input_min = ug[nodes[i]->GetDof(0,GsetEnum)];
		if(ug[nodes[i]->GetDof(0,GsetEnum)] > input_max) input_max = ug[nodes[i]->GetDof(0,GsetEnum)];
	}

	/*Second loop to reassign min and max with local extrema*/
	for(int i=0;i<numnodes;i++){
		if(min[nodes[i]->GetDof(0,GsetEnum)]>input_min) min[nodes[i]->GetDof(0,GsetEnum)] = input_min;
		if(max[nodes[i]->GetDof(0,GsetEnum)]<input_max) max[nodes[i]->GetDof(0,GsetEnum)] = input_max;
	}
}
/*}}}*/
void       Element::GetInputValue(bool* pvalue,int inputenum){/*{{{*/

	this->inputs->GetInputValue(pvalue,inputenum,this->lid);

}/*}}}*/
void       Element::GetInputValue(int* pvalue,int inputenum){/*{{{*/

	this->inputs->GetInputValue(pvalue,inputenum,this->lid);

}/*}}}*/
void       Element::GetInputValue(IssmDouble* pvalue,int inputenum){/*{{{*/

	this->inputs->GetInputValue(pvalue,inputenum,this->lid);

}/*}}}*/
void       Element::GetInputValue(IssmDouble* pvalue,Gauss* gauss,int inputenum){/*{{{*/

	Input* input=this->GetInput(inputenum);
	if(!input) _error_("Input " << EnumToStringx(inputenum) << " not found in element");
	input->GetInputValue(pvalue,gauss);

}/*}}}*/
Node*      Element::GetNode(int nodeindex){/*{{{*/
	_assert_(nodeindex>=0);
	_assert_(nodeindex<this->GetNumberOfNodes(this->element_type));
	return this->nodes[nodeindex];
}/*}}}*/
int        Element::GetNodeIndex(Node* node){/*{{{*/

	_assert_(this->nodes);
	int numnodes = this->GetNumberOfNodes(this->element_type);

	for(int i=0;i<numnodes;i++){
		if(node==nodes[i]) return i;
	}
	_error_("Node provided not found among element nodes");

}
/*}}}*/
void       Element::GetNodesLidList(int* lidlist){/*{{{*/

	_assert_(lidlist);
	_assert_(nodes);
	int numnodes = this->GetNumberOfNodes();
	for(int i=0;i<numnodes;i++){
		lidlist[i]=nodes[i]->Lid();
	}
}
/*}}}*/
void       Element::GetNodesSidList(int* sidlist){/*{{{*/

	_assert_(sidlist);
	_assert_(nodes);
	int numnodes = this->GetNumberOfNodes();
	for(int i=0;i<numnodes;i++){
		sidlist[i]=nodes[i]->Sid();
	}
}
/*}}}*/
void       Element::GetPhi(IssmDouble* phi, IssmDouble*  epsilon, IssmDouble viscosity){/*{{{*/
	/*Compute deformational heating from epsilon and viscosity */

	IssmDouble epsilon_matrix[3][3];
	IssmDouble epsilon_eff;
	IssmDouble epsilon_sqr[3][3];

	/* Build epsilon matrix */
	epsilon_matrix[0][0]=epsilon[0];
	epsilon_matrix[1][0]=epsilon[3];
	epsilon_matrix[2][0]=epsilon[4];
	epsilon_matrix[0][1]=epsilon[3];
	epsilon_matrix[1][1]=epsilon[1];
	epsilon_matrix[2][1]=epsilon[5];
	epsilon_matrix[0][2]=epsilon[4];
	epsilon_matrix[1][2]=epsilon[5];
	epsilon_matrix[2][2]=epsilon[2];

	/* Effective value of epsilon_matrix */
	epsilon_sqr[0][0]=epsilon_matrix[0][0]*epsilon_matrix[0][0];
	epsilon_sqr[1][0]=epsilon_matrix[1][0]*epsilon_matrix[1][0];
	epsilon_sqr[2][0]=epsilon_matrix[2][0]*epsilon_matrix[2][0];
	epsilon_sqr[0][1]=epsilon_matrix[0][1]*epsilon_matrix[0][1];
	epsilon_sqr[1][1]=epsilon_matrix[1][1]*epsilon_matrix[1][1];
	epsilon_sqr[2][1]=epsilon_matrix[2][1]*epsilon_matrix[2][1];
	epsilon_sqr[0][2]=epsilon_matrix[0][2]*epsilon_matrix[0][2];
	epsilon_sqr[1][2]=epsilon_matrix[1][2]*epsilon_matrix[1][2];
	epsilon_sqr[2][2]=epsilon_matrix[2][2]*epsilon_matrix[2][2];
	epsilon_eff=1/sqrt(2.)*sqrt(epsilon_sqr[0][0]+epsilon_sqr[0][1]+ epsilon_sqr[0][2]+ epsilon_sqr[1][0]+ epsilon_sqr[1][1]+ epsilon_sqr[1][2]+ epsilon_sqr[2][0]+ epsilon_sqr[2][1]+ epsilon_sqr[2][2]);

	/*Phi = Tr(sigma * eps)
	 *    = Tr(sigma'* eps)
	 *    = 2 * eps_eff * sigma'_eff
	 *    = 4 * mu * eps_eff ^2*/
	*phi=4.*epsilon_eff*epsilon_eff*viscosity;
}
/*}}}*/
void       Element::GetSolutionFromInputsOneDof(Vector<IssmDouble>* solution, int enum_type){/*{{{*/

	int        *doflist = NULL;
	IssmDouble  value;

	/*Fetch number of nodes for this finite element*/
	int numnodes = this->GetNumberOfNodes();

	/*Fetch dof list and allocate solution vector*/
	GetDofList(&doflist,NoneApproximationEnum,GsetEnum, true);
	IssmDouble* values = xNew<IssmDouble>(numnodes);

	/*Get inputs*/
	Input* enum_input=this->GetInput(enum_type); _assert_(enum_input);

	/*Ok, we have the values, fill in the array: */
	Gauss* gauss=this->NewGauss();
	for(int i=0;i<numnodes;i++){
		gauss->GaussNode(this->element_type,i);

		enum_input->GetInputValue(&value,gauss);
		values[i]=value;
	}

	solution->SetValues(numnodes,doflist,values,INS_VAL);

	/*Free resources:*/
	xDelete<int>(doflist);
	xDelete<IssmDouble>(values);
	delete gauss;
}
/*}}}*/
void       Element::GetVectorFromInputs(Vector<IssmDouble>* vector,int input_enum,int type){/*{{{*/

	switch(type){
		case ElementSIdEnum:{
         IssmDouble  value;
			Input* input=this->GetInput(input_enum); _assert_(input);
			input->GetInputAverage(&value);
			vector->SetValue(this->sid,value,INS_VAL);
                          }
			break;
		case VertexPIdEnum:{
         int        doflist[MAXVERTICES];
         IssmDouble values[MAXVERTICES];
			const int  NUM_VERTICES = this->GetNumberOfVertices();
			/*Fill in values*/
			this->GetVerticesPidList(&doflist[0]);
			/*Take care of Clones*/
			for(int i=0;i<NUM_VERTICES;i++){
				if(vertices[i]->clone) doflist[i] = -1;
			}
			this->GetInputListOnVertices(&values[0],input_enum);
			vector->SetValues(NUM_VERTICES,doflist,values,INS_VAL);
                         }
			break;
		case VertexSIdEnum:{
         int        doflist[MAXVERTICES];
         IssmDouble values[MAXVERTICES];
			const int  NUM_VERTICES = this->GetNumberOfVertices();
			/*Fill in values*/
			this->GetVerticesSidList(doflist);
			this->GetInputListOnVertices(values,input_enum);
			vector->SetValues(NUM_VERTICES,doflist,values,INS_VAL);
                         }
			break;
		case NodesEnum:{
         int         numnodes = this->GetNumberOfNodes();
         int        *doflist  = xNew<int>(numnodes);
         IssmDouble *values   = xNew<IssmDouble>(numnodes);
			/*Fill in values*/
			this->GetInputListOnNodes(values,input_enum);
			this->GetDofList(&doflist,NoneApproximationEnum,GsetEnum);
			vector->SetValues(numnodes,doflist,values,INS_VAL);
         xDelete<int>(doflist);
         xDelete<IssmDouble>(values);
                     }
			break;
		case NodeSIdEnum:{
         int         numnodes = this->GetNumberOfNodes();
         int        *doflist  = xNew<int>(numnodes);
         IssmDouble *values   = xNew<IssmDouble>(numnodes);
			/*Fill in values*/
			this->GetNodesSidList(doflist);
			this->GetInputListOnNodes(values,input_enum);
			vector->SetValues(numnodes,doflist,values,INS_VAL);
         xDelete<int>(doflist);
         xDelete<IssmDouble>(values);
                       }
			break;
		default:
			_error_("type " << type << " (" << EnumToStringx(type) << ") not implemented yet");
	}

}/*}}}*/
void       Element::GetVectorFromInputs(Vector<IssmDouble>* vector,int input_enum,int type, IssmDouble time){/*{{{*/

	switch(type){
		case VertexPIdEnum:{
			int        doflist[MAXVERTICES];
			IssmDouble values[MAXVERTICES];
			const int  NUM_VERTICES = this->GetNumberOfVertices();
			/*Fill in values*/
			this->GetVerticesPidList(doflist);
			this->GetInputListOnVerticesAtTime(values,input_enum,time);
			vector->SetValues(NUM_VERTICES,doflist,&values[0],INS_VAL);
								 }
			break;
		case VertexSIdEnum:{
			int        doflist[MAXVERTICES];
			IssmDouble values[MAXVERTICES];
			const int  NUM_VERTICES = this->GetNumberOfVertices();
			/*Fill in values*/
			this->GetVerticesSidList(doflist);
			this->GetInputListOnVerticesAtTime(values,input_enum,time);
			vector->SetValues(NUM_VERTICES,doflist,&values[0],INS_VAL);
								 }
			break;
		default:
			_error_("type " << type << " (" << EnumToStringx(type) << ") not implemented yet");
	}
}
/*}}}*/
void       Element::GetVerticesLidList(int* lidlist){/*{{{*/

	const int NUM_VERTICES = this->GetNumberOfVertices();
	for(int i=0;i<NUM_VERTICES;i++) lidlist[i]=vertices[i]->Lid();

}
/*}}}*/
void       Element::GetVerticesPidList(int* pidlist){/*{{{*/

	const int NUM_VERTICES = this->GetNumberOfVertices();
	for(int i=0;i<NUM_VERTICES;i++) pidlist[i]=vertices[i]->Pid();

}
/*}}}*/
void       Element::GetVerticesConnectivityList(int* connectivity){/*{{{*/

	const int NUM_VERTICES = this->GetNumberOfVertices();
	for(int i=0;i<NUM_VERTICES;i++) connectivity[i]=this->vertices[i]->Connectivity();
}
/*}}}*/
void       Element::GetVerticesCoordinates(IssmDouble** pxyz_list){/*{{{*/

	const int NUM_VERTICES = this->GetNumberOfVertices();

	IssmDouble* xyz_list    = xNew<IssmDouble>(NUM_VERTICES*3);
	::GetVerticesCoordinates(xyz_list,this->vertices,NUM_VERTICES);

	*pxyz_list = xyz_list;

}/*}}}*/
void       Element::GetVerticesSidList(int* sidlist){/*{{{*/

	const int NUM_VERTICES = this->GetNumberOfVertices();
	for(int i=0;i<NUM_VERTICES;i++) sidlist[i]=this->vertices[i]->Sid();
}
/*}}}*/
IssmDouble Element::GetXcoord(IssmDouble* xyz_list,Gauss* gauss){/*{{{*/

	/*output*/
	IssmDouble x;

	/*Create list of x*/
	const int NUM_VERTICES = this->GetNumberOfVertices();
	IssmDouble x_list[MAXVERTICES];
	for(int i=0;i<NUM_VERTICES;i++) x_list[i]=xyz_list[i*3+0];

	/*Get value at gauss point*/
	ValueP1OnGauss(&x,&x_list[0],gauss);
	return x;
}/*}}}*/
IssmDouble Element::GetYcoord(IssmDouble* xyz_list,Gauss* gauss){/*{{{*/

	/*output*/
	IssmDouble y;

	/*Create list of y*/
	const int NUM_VERTICES = this->GetNumberOfVertices();
	IssmDouble y_list[MAXVERTICES];
	for(int i=0;i<NUM_VERTICES;i++) y_list[i]=xyz_list[i*3+1];

	/*Get value at gauss point*/
	ValueP1OnGauss(&y,&y_list[0],gauss);
	return y;
}/*}}}*/
IssmDouble Element::GetZcoord(IssmDouble* xyz_list,Gauss* gauss){/*{{{*/

	/*output*/
	IssmDouble z;

	/*Create list of z*/
	const int NUM_VERTICES = this->GetNumberOfVertices();
	IssmDouble z_list[MAXVERTICES];
	for(int i=0;i<NUM_VERTICES;i++) z_list[i]=xyz_list[i*3+2];

	/*Get value at gauss point*/
	ValueP1OnGauss(&z,&z_list[0],gauss);
	return z;
}/*}}}*/
void       Element::GradientIndexing(int* indexing,int control_index){/*{{{*/

	/*Get number of controls*/
	int num_controls;
	bool isautodiff;
	int* N=NULL;
	int* M=NULL;
	int* interp=NULL;
	parameters->FindParam(&num_controls,InversionNumControlParametersEnum);
	parameters->FindParam(&N,NULL,ControlInputSizeNEnum);
	parameters->FindParam(&M,NULL,ControlInputSizeMEnum);
	parameters->FindParam(&interp,NULL,ControlInputInterpolationEnum);

	/*Get number of vertices*/
	const int NUM_VERTICES = this->GetNumberOfVertices();

	/*Get starting index of gradient for this control*/
	int start = 0;
	for(int n=0;n<control_index;n++) start+=N[n]*M[n];

	/*Create index*/
	if(interp[control_index]==P1Enum){
		for(int n=0;n<N[control_index];n++){
			for(int i=0;i<NUM_VERTICES;i++){
				indexing[i+n*NUM_VERTICES]= start + n*M[control_index] + this->vertices[i]->Sid();
			}
		}
	}
	else if(interp[control_index]==P0Enum){
		for(int n=0;n<N[control_index];n++){
				indexing[n]= start + n*M[control_index] + this->Sid();
		}
	}
	else{
		_error_("interpolation not supported");
	}

	xDelete<int>(M);
	xDelete<int>(N);
	xDelete<int>(interp);
}
/*}}}*/
IssmDouble Element::GroundedArea(IssmDouble* mask, bool scaled){/*{{{*/

	/*Retrieve values of the mask defining the element: */
	for(int i=0;i<this->GetNumberOfVertices();i++){
		if(mask[this->vertices[i]->Sid()]<=0.){
			return 0.;
		}
	}

	/*Return: */
	return this->GroundedArea(scaled);
}
/*}}}*/
IssmDouble Element::GroundinglineMassFlux(IssmDouble* mask, bool scaled){/*{{{*/

	/*Retrieve values of the mask defining the element: */
	for(int i=0;i<this->GetNumberOfVertices();i++){
		if(mask[this->vertices[i]->Sid()]<=0.){
			return 0.;
		}
	}

	/*Return: */
	return this->GroundinglineMassFlux(scaled);
}
/*}}}*/
bool       Element::HasNodeOnBase(){/*{{{*/
	Input* input=this->GetInput(MeshVertexonbaseEnum); _assert_(input);
	return (input->GetInputMax()>0.);
}/*}}}*/
bool       Element::HasNodeOnSurface(){/*{{{*/
	Input* input=this->GetInput(MeshVertexonsurfaceEnum); _assert_(input);
	return (input->GetInputMax()>0.);
}/*}}}*/
IssmDouble Element::IceMass(bool scaled){/*{{{*/

	IssmDouble rho_ice;

	if(!IsIceInElement())return 0.; //do not contribute to the volume of the ice!

	/*recover ice density: */
	rho_ice=FindParam(MaterialsRhoIceEnum);

	return rho_ice*this->IceVolume(scaled);
}
/*}}}*/
IssmDouble Element::IceMass(IssmDouble* mask, bool scaled){/*{{{*/

	/*Retrieve values of the mask defining the element: */
	for(int i=0;i<this->GetNumberOfVertices();i++){
		if(mask[this->vertices[i]->Sid()]<=0.){
			return 0.;
		}
	}

	/*Return: */
	return this->IceMass(scaled);
}
/*}}}*/
IssmDouble Element::IceVolume(IssmDouble* mask, bool scaled){/*{{{*/

	/*Retrieve values of the mask defining the element: */
	for(int i=0;i<this->GetNumberOfVertices();i++){
		if(mask[this->vertices[i]->Sid()]<=0.){
			return 0.;
		}
	}

	/*Return: */
	return this->IceVolume(scaled);
}
/*}}}*/
IssmDouble Element::IceVolumeAboveFloatation(IssmDouble* mask, bool scaled){/*{{{*/

	/*Retrieve values of the mask defining the element: */
	for(int i=0;i<this->GetNumberOfVertices();i++){
		if(mask[this->vertices[i]->Sid()]<=0.){
			return 0.;
		}
	}

	/*Return: */
	return this->IceVolumeAboveFloatation(scaled);
}
/*}}}*/
int        Element::Id(){/*{{{*/

	return this->id;

}
/*}}}*/
void       Element::InputCreate(IssmDouble* vector,Inputs* inputs,IoModel* iomodel,int M,int N,int vector_type,int vector_enum,int code){/*{{{*/

	/*Branch on type of vector: nodal or elementary: */
	if(vector_type==1){ //nodal vector

		const int NUM_VERTICES = this->GetNumberOfVertices();

		int        vertexids[MAXVERTICES];
		int        vertexlids[MAXVERTICES];
		IssmDouble values[MAXVERTICES];

		/*Recover vertices ids needed to initialize inputs*/
		_assert_(iomodel->elements);
		for(int i=0;i<NUM_VERTICES;i++){
			vertexids[i] =reCast<int>(iomodel->elements[NUM_VERTICES*this->Sid()+i]); //ids for vertices are in the elements array from Matlab
			vertexlids[i]=iomodel->my_vertices_lids[vertexids[i]-1];
		}

		/*Are we in transient or static? */
		if(M==1){
			if(N!=1) _error_("Size of Input "<<EnumToStringx(vector_enum)<<" "<<M<<"x"<<N<<" not supported");
			_assert_(N==1);
			this->SetElementInput(inputs,vector_enum,vector[0]);
		}

		else if(M==iomodel->numberofvertices){
			if(N!=1) _error_("Size of Input "<<EnumToStringx(vector_enum)<<" "<<M<<"x"<<N<<" not supported");
			for(int i=0;i<NUM_VERTICES;i++) values[i]=vector[vertexids[i]-1];
			this->SetElementInput(inputs,NUM_VERTICES,vertexlids,values,vector_enum);
		}
		else if(M==iomodel->numberofvertices+1){
			/*create transient input: */
			IssmDouble* times = xNew<IssmDouble>(N);
			for(int t=0;t<N;t++) times[t] = vector[(M-1)*N+t];
			inputs->SetTransientInput(vector_enum,times,N);
			TransientInput* transientinput = inputs->GetTransientInput(vector_enum);
			for(int t=0;t<N;t++){
				for(int i=0;i<NUM_VERTICES;i++) values[i]=vector[N*(vertexids[i]-1)+t];
				switch(this->ObjectEnum()){
					case TriaEnum:  transientinput->AddTriaTimeInput( t,NUM_VERTICES,vertexlids,&values[0],P1Enum); break;
					case PentaEnum: transientinput->AddPentaTimeInput(t,NUM_VERTICES,vertexlids,&values[0],P1Enum); break;
					default: _error_("Not implemented yet");
				}
			}
			xDelete<IssmDouble>(times);
		}
		else if(M==iomodel->numberofelements){

			/*This is a Patch!*/
			IssmDouble* evalues = xNew<IssmDouble>(N);
			for(int j=0;j<N;j++) evalues[j]=vector[this->Sid()*N+j];

			if (N==this->GetNumberOfNodes(P1Enum)){
				this->SetElementInput(inputs,NUM_VERTICES,vertexlids,evalues,vector_enum);
			}
			else if(N==this->GetNumberOfNodes(P0Enum)){
				this->SetElementInput(inputs,vector_enum,evalues[0]);
			}
			else if(N==this->GetNumberOfNodes(P1xP2Enum)){ _assert_(this->ObjectEnum()==PentaEnum);
				inputs->SetPentaInput(vector_enum,P1xP2Enum,this->lid,N,evalues);
			}
			else if(N==this->GetNumberOfNodes(P1xP3Enum)){ _assert_(this->ObjectEnum()==PentaEnum);
				inputs->SetPentaInput(vector_enum,P1xP3Enum,this->lid,N,evalues);
			}
			else{
				_error_("Size of Input "<<EnumToStringx(vector_enum)<<" "<<M<<"x"<<N<<" not supported");
			}
			xDelete<IssmDouble>(evalues);

		}
		else{
			_error_("Size of Input "<<EnumToStringx(vector_enum)<<" "<<M<<"x"<<N<<" not supported");
		}
	}
	else if(vector_type==2){ //element vector

		/*Are we in transient or static? */
		if(M==1){
			if(N!=1) _error_("Size of Input "<<EnumToStringx(vector_enum)<<" "<<M<<"x"<<N<<" not supported");
			this->SetElementInput(inputs,vector_enum,vector[0]);
		}
		else if(M==2){
			/*create transient input: */
			IssmDouble* times = xNew<IssmDouble>(N);
			for(int t=0;t<N;t++) times[t] = vector[(M-1)*N+t];

			inputs->SetTransientInput(vector_enum,times,N);
			TransientInput* transientinput = inputs->GetTransientInput(vector_enum);

			for(int t=0;t<N;t++){
				IssmDouble value=vector[t]; //values are on the first line, times are on the second line
				switch(this->ObjectEnum()){
					case TriaEnum:  transientinput->AddTriaTimeInput( t,1,&(this->lid),&value,P0Enum); break;
					case PentaEnum: transientinput->AddPentaTimeInput(t,1,&(this->lid),&value,P0Enum); break;
					default: _error_("Not implemented yet");
				}
			}
			xDelete<IssmDouble>(times);
		}
		else if(M==iomodel->numberofelements){
			if(N!=1) _error_("Size of Input "<<EnumToStringx(vector_enum)<<" "<<M<<"x"<<N<<" not supported");
			if (code==5){ //boolean
				this->SetBoolInput(inputs,vector_enum,reCast<bool>(vector[this->Sid()]));
			}
			else if (code==6){ //integer
				this->SetIntInput(inputs,vector_enum,reCast<int>(vector[this->Sid()]));
			}
			else if (code==7){ //IssmDouble
				this->SetElementInput(inputs,vector_enum,vector[this->Sid()]);
			}
			else _error_("could not recognize nature of vector from code " << code);
		}
		else if(M==iomodel->numberofelements+1){
			/*create transient input: */
			IssmDouble* times = xNew<IssmDouble>(N);
			for(int t=0;t<N;t++) times[t] = vector[(M-1)*N+t];
			inputs->SetTransientInput(vector_enum,times,N);
			TransientInput* transientinput = inputs->GetTransientInput(vector_enum);
			for(int t=0;t<N;t++){
				IssmDouble value=vector[N*this->Sid()+t];
				switch(this->ObjectEnum()){
					case TriaEnum:  transientinput->AddTriaTimeInput( t,1,&(this->lid),&value,P0Enum); break;
					case PentaEnum: transientinput->AddPentaTimeInput(t,1,&(this->lid),&value,P0Enum); break;
					default: _error_("Not implemented yet");
				}
			}
			xDelete<IssmDouble>(times);
		}

		else{
			_error_("Size of Input "<<EnumToStringx(vector_enum)<<" "<<M<<"x"<<N<<" not supported");
		}
	}
	else if(vector_type==3){ //Double array matrix

		/*For right now we are static */
		if(M==iomodel->numberofelements){
			IssmDouble* layers = xNewZeroInit<IssmDouble>(N);
			for(int t=0;t<N;t++) layers[t] = vector[N*this->Sid()+t];
			inputs->SetArrayInput(vector_enum,this->lid,layers,N);
			xDelete<IssmDouble>(layers);
		}
		else{
			_error_("Size of Input "<<EnumToStringx(vector_enum)<<" "<<M<<"x"<<N<<" not supported");
		}
	}
	else{
		_error_("Cannot add input for vector type " << vector_type << " (not supported)");
	}
}
/*}}}*/
void       Element::InputCreateP1FromConstant(Inputs* inputs,IoModel* iomodel,IssmDouble value_in,int vector_enum){/*{{{*/

	const int NUM_VERTICES = this->GetNumberOfVertices();

	int        vertexids[MAXVERTICES];
	int        vertexlids[MAXVERTICES];
	IssmDouble values[MAXVERTICES];

	/*Recover vertices ids needed to initialize inputs*/
	_assert_(iomodel->elements);
	for(int i=0;i<NUM_VERTICES;i++){
		vertexids[i] =reCast<int>(iomodel->elements[NUM_VERTICES*this->Sid()+i]); //ids for vertices are in the elements array from Matlab
		vertexlids[i]=iomodel->my_vertices_lids[vertexids[i]-1];
	}

	for(int i=0;i<NUM_VERTICES;i++) values[i]=value_in;
	this->SetElementInput(inputs,NUM_VERTICES,vertexlids,values,vector_enum);

}
/*}}}*/
void       Element::InputCreateP0FromConstant(Inputs* inputs,IoModel* iomodel,IssmDouble value_in,int vector_enum){/*{{{*/

	this->SetElementInput(inputs,vector_enum,value_in);

}
/*}}}*/
void       Element::ControlInputCreate(IssmDouble* vector,IssmDouble* min_vector,IssmDouble* max_vector,Inputs* inputs,IoModel* iomodel,int M,int N,IssmDouble scale,int input_enum,int id){/*{{{*/

	/*Intermediaries*/
	const int numvertices = this->GetNumberOfVertices();
	IssmDouble  value,value_min,value_max;

	int         vertexids[MAXVERTICES];
	int         vertexlids[MAXVERTICES];
	IssmDouble  values[MAXVERTICES];
	IssmDouble  values_min[MAXVERTICES];
	IssmDouble  values_max[MAXVERTICES];

	/*Some sanity checks*/
	_assert_(vector);
	_assert_(min_vector);
	_assert_(max_vector);

	/*Recover vertices ids needed to initialize inputs*/
	_assert_(iomodel->elements);
	for(int i=0;i<numvertices;i++){
		vertexids[i]=reCast<int>(iomodel->elements[numvertices*this->Sid()+i]); //ids for vertices are in the elements array from Matlab
		vertexlids[i]=iomodel->my_vertices_lids[vertexids[i]-1];
	}

	/*Are we in transient or static? */
	if(M==iomodel->numberofvertices && N==1){
		for(int i=0;i<numvertices;i++){
			values[i]     = vector[vertexids[i]-1];
			values_min[i] = scale*min_vector[vertexids[i]-1];
			values_max[i] = scale*max_vector[vertexids[i]-1];
		}
		this->AddControlInput(input_enum,inputs,iomodel,&values[0],&values_min[0],&values_max[0],P1Enum,id);
	}
	else if(M==iomodel->numberofelements && N==1){
		values[0]     = vector[this->Sid()];
		values_min[0] = scale*min_vector[this->Sid()];
		values_max[0] = scale*max_vector[this->Sid()];
		this->AddControlInput(input_enum,inputs,iomodel,&values[0],&values_min[0],&values_max[0],P0Enum,id);
	}
	else if(M==iomodel->numberofelements+1 && N>1){

		/*create transient input: */
		IssmDouble* times = xNew<IssmDouble>(N);
		for(int t=0;t<N;t++) times[t] = vector[(M-1)*N+t];
		inputs->SetTransientControlInput(input_enum,id,times,N);

		/*Get input*/
		ControlInput* control_input = inputs->GetControlInput(input_enum); _assert_(control_input);
		TransientInput* transientinput_val = control_input->GetTransientInput("value");
		TransientInput* transientinput_min = control_input->GetTransientInput("lowerbound");
		TransientInput* transientinput_max = control_input->GetTransientInput("upperbound");
		for(int t=0;t<N;t++){
			value     = vector[N*this->Sid()+t];
			value_min = scale*min_vector[N*this->Sid()+t];
			value_max = scale*max_vector[N*this->Sid()+t];
			switch(this->ObjectEnum()){
				case TriaEnum:
					transientinput_val->AddTriaTimeInput( t,1,&(this->lid),&value,P0Enum);
					transientinput_min->AddTriaTimeInput( t,1,&(this->lid),&value_min,P0Enum);
					transientinput_max->AddTriaTimeInput( t,1,&(this->lid),&value_max,P0Enum);
					break;
				case PentaEnum:
					//transientinput->AddPentaTimeInput(t,1,&(this->lid),&value,P0Enum);
					_error_("to be implemented");
					break;
				default: _error_("Not implemented yet");
			}
		}
		xDelete<IssmDouble>(times);
	}
	else if(M==iomodel->numberofvertices+1 && N>1){

		/*create transient input: */
		IssmDouble* times = xNew<IssmDouble>(N);
		for(int t=0;t<N;t++) times[t] = vector[(M-1)*N+t];
		inputs->SetTransientControlInput(input_enum,id,times,N);

		/*Get input*/
		ControlInput* control_input = inputs->GetControlInput(input_enum); _assert_(control_input);
		TransientInput* transientinput_val = control_input->GetTransientInput("value");
		TransientInput* transientinput_min = control_input->GetTransientInput("lowerbound");
		TransientInput* transientinput_max = control_input->GetTransientInput("upperbound");
		for(int t=0;t<N;t++){
			for(int i=0;i<numvertices;i++){
				values[i]     = vector[N*(vertexids[i]-1)+t];
				values_min[i] = scale*min_vector[N*(vertexids[i]-1)+t];
				values_max[i] = scale*max_vector[N*(vertexids[i]-1)+t];
			}
			switch(this->ObjectEnum()){
				case TriaEnum:
					transientinput_val->AddTriaTimeInput( t,numvertices,vertexlids,&values[0],P1Enum);
					transientinput_min->AddTriaTimeInput( t,numvertices,vertexlids,&values_min[0],P1Enum);
					transientinput_max->AddTriaTimeInput( t,numvertices,vertexlids,&values_max[0],P1Enum);
					break;
				case PentaEnum:
					//transientinput->AddPentaTimeInput(t,1,&(this->lid),&value,P0Enum);
					_error_("to be implemented");
					break;
				default: _error_("Not implemented yet");
			}
		}
		xDelete<IssmDouble>(times);
	}
	else _error_("not currently supported type of M and N attempted");
}/*}}}*/
void       Element::DatasetInputAdd(int enum_type,IssmDouble* vector,Inputs* inputs,IoModel* iomodel,int M,int N,int vector_type,int input_enum,int input_id){/*{{{*/
	/*enum_type: the name of the DatasetInput (eg Outputdefinition1)
	 * vector: information being stored (eg observations)
	 * vector_type: is if by element or by vertex
	 * input_enum: is the name of the vector being stored
	 */

	/*Branch on type of vector: nodal or elementary: */
	if(vector_type==1){ //nodal vector

		const int NUM_VERTICES = this->GetNumberOfVertices();

		int        vertexids[MAXVERTICES];
		int        vertexlids[MAXVERTICES];
		IssmDouble values[MAXVERTICES];

		/*Recover vertices ids needed to initialize inputs*/
		_assert_(iomodel->elements);
		for(int i=0;i<NUM_VERTICES;i++){
			vertexids[i] =reCast<int>(iomodel->elements[NUM_VERTICES*this->Sid()+i]); //ids for vertices are in the elements array from Matlab
			vertexlids[i]=iomodel->my_vertices_lids[vertexids[i]-1];
		}

		/*Are we in transient or static? */
		if(M==1){
			values[0]=vector[0];
			//this->AddInput(vector_enum,values,P0Enum);
			_error_("not implemented yet");
		}
		else if(M==iomodel->numberofvertices){
			for(int i=0;i<NUM_VERTICES;i++) values[i]=vector[vertexids[i]-1];
			switch(this->ObjectEnum()){
				case TriaEnum:  inputs->SetTriaDatasetInput(enum_type,input_id,P1Enum,NUM_VERTICES,vertexlids,values); break;
				case PentaEnum: inputs->SetPentaDatasetInput(enum_type,input_id,P1Enum,NUM_VERTICES,vertexlids,values); break;
				default: _error_("Not implemented yet for "<<this->ObjectEnum());
			}
		}
		else if(M==iomodel->numberofvertices+1){
			/*create transient input: */
			IssmDouble* times = xNew<IssmDouble>(N);
			for(int t=0;t<N;t++) times[t] = vector[(M-1)*N+t];
			TransientInput* transientinput = inputs->SetDatasetTransientInput(enum_type,input_id,times,N);
			for(int t=0;t<N;t++){
				for(int i=0;i<NUM_VERTICES;i++) values[i]=vector[N*(vertexids[i]-1)+t];
				switch(this->ObjectEnum()){
					case TriaEnum:  transientinput->AddTriaTimeInput( t,NUM_VERTICES,vertexlids,values,P1Enum); break;
					case PentaEnum: transientinput->AddPentaTimeInput(t,NUM_VERTICES,vertexlids,values,P1Enum); break;
					default: _error_("Not implemented yet");
				}
			}
			xDelete<IssmDouble>(times);
		}
		else{
			_error_("not implemented yet (M="<<M<<")");
		}
	}
	else if(vector_type==2){ //element vector
		_error_("not supported");

		IssmDouble value;

		/*Are we in transient or static? */
		if(M==iomodel->numberofelements){
			_error_("not implemented");
		}
		else if(M==iomodel->numberofelements+1){
			_error_("not supported");
		}
		else _error_("element vector is either numberofelements or numberofelements+1 long. Field provided (" << EnumToStringx(input_enum) << ") is " << M << " long");
	}
	else if(vector_type==3){ //element vector
		_error_("not supported");

		///*For right now we are static */
		//if(M==iomodel->numberofelements){
		//	/*create transient input: */
		//	IssmDouble* layers = xNewZeroInit<IssmDouble>(N);;
		//	for(t=0;t<N;t++) layers[t] = vector[N*this->Sid()+t];
		//	DoubleArrayInput* arrayinput=new DoubleArrayInput(input_enum,layers,N);
		//	xDelete<IssmDouble>(layers);
		//}
		//else _error_("element vector is either numberofelements or numberofelements+1 long. Field provided (" << EnumToStringx(input_enum) << ") is " << M << " long");
	}
	else{
		_error_("Cannot add input for vector type " << vector_type << " (not supported)");
	}
}
/*}}}*/
void       Element::InputUpdateFromConstant(int constant, int name){/*{{{*/

	/*Check that name is an element input*/
	if(!IsInputEnum(name)) _error_("Enum "<<EnumToStringx(name)<<" is not in IsInput");

	/*update input*/
	this->SetIntInput(this->inputs,name,constant);
}
/*}}}*/
void       Element::InputUpdateFromConstant(IssmDouble constant, int name){/*{{{*/

	/*Check that name is an element input*/
	if(!IsInputEnum(name)) return;

	/*update input*/
	this->SetElementInput(name,constant);
}
/*}}}*/
void       Element::InputUpdateFromConstant(IssmDouble constant, int name, int type){/*{{{*/

	/*Check that name is an element input*/
	if(!IsInputEnum(name)) return;

	/*update input*/
	this->SetElementInput(name,constant,type);
}
/*}}}*/
void       Element::InputUpdateFromConstant(bool constant, int name){/*{{{*/

	/*Check that name is an element input*/
	if(!IsInputEnum(name)) return;

	/*update input*/
	this->SetBoolInput(this->inputs,name,constant);
}
/*}}}*/
bool       Element::IsOnSurface(){/*{{{*/
	return this->isonsurface;
}/*}}}*/
bool       Element::IsOnBase(){/*{{{*/
	return this->isonbase;
}/*}}}*/
bool       Element::IsAllFloating(){/*{{{*/

	bool shelf;
	int  migration_style;
	parameters->FindParam(&migration_style,GroundinglineMigrationEnum);

	Input* input = this->GetInput(MaskOceanLevelsetEnum); _assert_(input);

	if(migration_style==SubelementMigrationEnum){ //Floating if all nodes are floating
		if(input->GetInputMax() <= 0.) shelf=true;
		else shelf=false;
	}
	else if(migration_style==ContactEnum){
		if(input->GetInputMin() < 0.) shelf=true;
		else shelf=false;
	}
	else if(migration_style==NoneEnum || migration_style==AggressiveMigrationEnum || migration_style==SoftMigrationEnum || migration_style==GroundingOnlyEnum){ //Floating if all nodes are floating
		if(input->GetInputMin() > 0.) shelf=false;
		else shelf=true;
	}
	else _error_("migration_style not implemented yet");

	return shelf;
}/*}}}*/
bool       Element::IsAllGrounded(){/*{{{*/

	Input* input=this->GetInput(MaskOceanLevelsetEnum); _assert_(input);
	if(input->GetInputMin() < 0.){
		return false;
	}
	else{
		return true;
	}
}/*}}}*/
bool       Element::IsGrounded(){/*{{{*/
	/*At least ONE node is grounded (partially grounded returns true)*/

	Input* input=this->GetInput(MaskOceanLevelsetEnum); _assert_(input);
	if(input->GetInputMax() > 0.){
		return true;
	}
	else{
		return false;
	}
}/*}}}*/
bool       Element::IsIceInElement(){/*{{{*/
	Input* input=this->GetInput(MaskIceLevelsetEnum); _assert_(input);
	return (input->GetInputMin()<0.);
}
/*}}}*/
bool       Element::IsIceOnlyInElement(){/*{{{*/
	Input* input=this->GetInput(MaskIceLevelsetEnum); _assert_(input);
	return (input->GetInputMax()<=0.);
}
/*}}}*/
bool       Element::IsLandInElement(){/*{{{*/
	Input* input=this->GetInput(MaskOceanLevelsetEnum); _assert_(input);
	return (input->GetInputMax()>0.);
}
/*}}}*/
bool       Element::IsOceanInElement(){/*{{{*/
	Input* input=this->GetInput(MaskOceanLevelsetEnum); _assert_(input);
	return (input->GetInputMin()<0.);
}
/*}}}*/
bool       Element::IsOceanOnlyInElement(){/*{{{*/
	Input* input=this->GetInput(MaskOceanLevelsetEnum); _assert_(input);
	return (input->GetInputMax()<=0.);
}
/*}}}*/
bool       Element::IsAllMinThicknessInElement(){/*{{{*/

	IssmDouble minthickness = this->FindParam(MasstransportMinThicknessEnum);

	Input* input=this->GetInput(ThicknessEnum); _assert_(input);
	if(input->GetInputMax()<=(minthickness+0.00000001)){
		return true;
	}
	else{
                return false;
        }
}
/*}}}*/
void       Element::Ismip6FloatingiceMeltingRate(){/*{{{*/

	if(!this->IsIceInElement() || !this->IsAllFloating() || !this->IsOnBase()) return;

	int         basinid,num_basins,M,N;
	IssmDouble  tf,gamma0,base,delta_t_basin,mean_tf_basin,absval,meltanomaly;
	bool        islocal;
	IssmDouble* delta_t = NULL;
	IssmDouble* mean_tf = NULL;
	IssmDouble* depths  = NULL;

	/*Allocate some arrays*/
	const int numvertices = this->GetNumberOfVertices();
	IssmDouble basalmeltrate[MAXVERTICES];

	/*Get variables*/
	IssmDouble rhoi = this->FindParam(MaterialsRhoIceEnum);
	IssmDouble rhow = this->FindParam(MaterialsRhoSeawaterEnum);
	IssmDouble lf   = this->FindParam(MaterialsLatentheatEnum);
	IssmDouble cp   = this->FindParam(MaterialsMixedLayerCapacityEnum);

	/*Hard code sea water density to be consistent with ISMIP6 documentation*/
	rhow = 1028.;

	/* Get parameters and inputs */
	this->GetInputValue(&basinid,BasalforcingsIsmip6BasinIdEnum);
	this->parameters->FindParam(&num_basins,BasalforcingsIsmip6NumBasinsEnum);
	this->parameters->FindParam(&gamma0,BasalforcingsIsmip6Gamma0Enum);
	this->parameters->FindParam(&delta_t,&M,BasalforcingsIsmip6DeltaTEnum);    _assert_(M==num_basins);
	this->parameters->FindParam(&islocal,BasalforcingsIsmip6IsLocalEnum);
	if(!islocal) {
		this->parameters->FindParam(&mean_tf,&N,BasalforcingsIsmip6AverageTfEnum); _assert_(N==num_basins);
	}
	Input* tf_input = this->GetInput(BasalforcingsIsmip6TfShelfEnum);              _assert_(tf_input);
	Input* meltanomaly_input = this->GetInput(BasalforcingsIsmip6MeltAnomalyEnum); _assert_(meltanomaly_input);
	delta_t_basin = delta_t[basinid];
	if(!islocal) mean_tf_basin = mean_tf[basinid];

	/*Compute melt rate for Local and Nonlocal parameterizations*/
	Gauss* gauss=this->NewGauss();
	for(int i=0;i<numvertices;i++){
		gauss->GaussVertex(i);
		tf_input->GetInputValue(&tf,gauss);
		meltanomaly_input->GetInputValue(&meltanomaly,gauss);
		if(!islocal) {
			absval = mean_tf_basin+delta_t_basin;
			if (absval<0) {absval = -1.*absval;}
			basalmeltrate[i] = gamma0*pow((rhow*cp)/(rhoi*lf),2)*(tf+delta_t_basin)*absval + meltanomaly;
		}
		else{
			basalmeltrate[i] = gamma0*pow((rhow*cp)/(rhoi*lf),2)*pow(max(tf+delta_t_basin,0.),2) + meltanomaly;
		}
	}

	/*Return basal melt rate*/
	this->AddInput(BasalforcingsFloatingiceMeltingRateEnum,basalmeltrate,P1DGEnum);

	/*Cleanup and return*/
	delete gauss;
	xDelete<IssmDouble>(delta_t);
	xDelete<IssmDouble>(mean_tf);
	xDelete<IssmDouble>(depths);

}/*}}}*/
void       Element::LapseRateBasinSMB(int numelevbins, IssmDouble* lapserates, IssmDouble* elevbins,IssmDouble* refelevation){/*{{{*/

	/*Variable declaration*/
   const int numvertices = this->GetNumberOfVertices();
   bool isadjustsmb = false;
	int basinid,bb1,bb2,mindex;
	IssmDouble ela,refelevation_b,time,dt,fracyear,yts;
   IssmDouble monthsteps[12]  = {0.,1./12,2./12,3./12,4./12,5./12,6./12,7./12,8./12,9./12,10./12,11./12};
   IssmDouble* surfacelist  = xNew<IssmDouble>(numvertices);
   IssmDouble* smblist      = xNew<IssmDouble>(numvertices);
   /* numelevbins values of lapse rates at current month */
	IssmDouble* lapserates_b = xNew<IssmDouble>(numelevbins);
   /* (numelevbins-1) limits between elevation bins at current month (be cautious with indexing) */
	IssmDouble* elevbins_b   = xNew<IssmDouble>(numelevbins-1);

	/*Find month of current time step*/
	this->parameters->FindParam(&yts,ConstantsYtsEnum);
   this->parameters->FindParam(&time,TimeEnum);
   this->parameters->FindParam(&dt,TimesteppingTimeStepEnum); 
   fracyear     = time/yts-floor(time/yts);
   for(int i=1;i<12;i++){
		if(fracyear>=monthsteps[i-1]) mindex = i-1;
	}
   if(fracyear>=monthsteps[11]) mindex = 11;

   /*Retrieve SMB values non-adjusted for SMB lapse rate*/
   this->GetInputListOnVertices(smblist,SmbMassBalanceEnum);
	/*Get surface elevation on vertices*/
	this->GetInputListOnVertices(surfacelist,SurfaceEnum);
   /*Get basin-specific parameters of the element*/
   this->GetInputValue(&basinid,SmbBasinsIdEnum);
   refelevation_b = refelevation[basinid];
	/*Retrieve bins and laps rates for this basin at this month*/
	for(int ii=0;ii<(numelevbins-1);ii++) elevbins_b[ii] = elevbins[basinid*(numelevbins-1)*12+mindex*(numelevbins-1)+ii];
	for(int ii=0;ii<numelevbins;ii++){
		lapserates_b[ii] = lapserates[basinid*numelevbins*12+mindex*numelevbins+ii];
		if(lapserates_b[ii]!=0) isadjustsmb=true;
	}

	/*Adjust SMB if any lapse rate value is non-zero*/
	if(isadjustsmb){

		_assert_(dt<yts);
	   for(int v=0;v<numvertices;v++){
	      /*Find elevation bin of Reference elevation and of Vertex*/
			bb1 = 0;
			bb2 = 0;
			for(int ii=0;ii<(numelevbins-1);ii++){
				if(surfacelist[v]>=elevbins_b[ii]) bb1 = ii+1;
				if(refelevation_b>=elevbins_b[ii]) bb2 = ii+1;
			}

			/*Vertex and Reference elevation in same elevation bin*/
			if(bb1==bb2){
				smblist[v] = smblist[v]+(surfacelist[v]-refelevation_b)*lapserates_b[bb2];
			}
			/*Vertex in lower elevation bin than Reference elevation*/
			if(bb1<bb2){
				smblist[v] = smblist[v]+(elevbins_b[bb2-1]-refelevation_b)*lapserates_b[bb2];
				for(int ii=bb2-1;ii>bb1;ii--) smblist[v] = smblist[v]+(elevbins_b[ii-1]-elevbins_b[ii])*lapserates_b[ii];
				smblist[v] = smblist[v]+(surfacelist[v]-elevbins_b[bb1])*lapserates_b[bb1];
			}
			/*Vertex in higher elevation bin than Reference elevation*/
			if(bb1>bb2){
				smblist[v] = smblist[v]+(elevbins_b[bb2]-refelevation_b)*lapserates_b[bb2];
				for(int ii=bb2+1;ii<bb1;ii++) smblist[v] = smblist[v]+(elevbins_b[ii]-elevbins_b[ii-1])*lapserates_b[ii];
				smblist[v] = smblist[v]+(surfacelist[v]-elevbins_b[bb1-1])*lapserates_b[bb1];
			}
		}
	}

   /*Add input to element*/
   this->AddInput(SmbMassBalanceEnum,smblist,P1Enum);

   /*Cleanup*/
   xDelete<IssmDouble>(lapserates_b);
   xDelete<IssmDouble>(elevbins_b);
   xDelete<IssmDouble>(surfacelist);
   xDelete<IssmDouble>(smblist);
}/*}}}*/
void       Element::LinearFloatingiceMeltingRate(){/*{{{*/

	const int NUM_VERTICES = this->GetNumberOfVertices();

	IssmDouble deepwaterel,upperwaterel,deepwatermelt,upperwatermelt;
	IssmDouble base[MAXVERTICES];
	IssmDouble perturbation[MAXVERTICES];
	IssmDouble values[MAXVERTICES];
	IssmDouble time;

	parameters->FindParam(&time,TimeEnum);
	parameters->FindParam(&deepwaterel,BasalforcingsDeepwaterElevationEnum,time);
	parameters->FindParam(&deepwatermelt,BasalforcingsDeepwaterMeltingRateEnum,time);
	parameters->FindParam(&upperwaterel,BasalforcingsUpperwaterElevationEnum,time);
	parameters->FindParam(&upperwatermelt,BasalforcingsUpperwaterMeltingRateEnum,time);
	_assert_(upperwaterel>deepwaterel);

	this->GetInputListOnVertices(&base[0],BaseEnum);
	this->GetInputListOnVertices(&perturbation[0],BasalforcingsPerturbationMeltingRateEnum);
	for(int i=0;i<NUM_VERTICES;i++){
		if(base[i]>=upperwaterel){
			values[i]=upperwatermelt;
		}
		else if (base[i]<deepwaterel){
			values[i]=deepwatermelt;
		}
		else{
			IssmDouble alpha = (base[i]-upperwaterel)/(deepwaterel-upperwaterel);
			values[i]=deepwatermelt*alpha+(1.-alpha)*upperwatermelt;
		}

		/*Add perturbation*/
		values[i] += perturbation[i];
	}

	this->AddInput(BasalforcingsFloatingiceMeltingRateEnum,values,P1Enum);

}/*}}}*/
void       Element::SpatialLinearFloatingiceMeltingRate(){/*{{{*/

	const int NUM_VERTICES = this->GetNumberOfVertices();

	IssmDouble alpha;
	IssmDouble deepwatermelt[MAXVERTICES];
	IssmDouble deepwaterel[MAXVERTICES];
	IssmDouble upperwatermelt[MAXVERTICES];
	IssmDouble upperwaterel[MAXVERTICES];
	IssmDouble perturbation[MAXVERTICES];
	IssmDouble base[MAXVERTICES];
	IssmDouble values[MAXVERTICES];

	this->GetInputListOnVertices(&base[0],BaseEnum);
	this->GetInputListOnVertices(&deepwatermelt[0],BasalforcingsSpatialDeepwaterMeltingRateEnum);
	this->GetInputListOnVertices(&deepwaterel[0],BasalforcingsSpatialDeepwaterElevationEnum);
	this->GetInputListOnVertices(&upperwatermelt[0],BasalforcingsSpatialUpperwaterMeltingRateEnum);
	this->GetInputListOnVertices(&upperwaterel[0],BasalforcingsSpatialUpperwaterElevationEnum);
   this->GetInputListOnVertices(&perturbation[0],BasalforcingsPerturbationMeltingRateEnum);

	for(int i=0;i<NUM_VERTICES;i++){
		if(base[i]>=upperwaterel[i]){
			values[i]=upperwatermelt[i];
		}
		else if (base[i]<deepwaterel[i]){
			values[i]=deepwatermelt[i];
		}
		else{
			alpha = (base[i]-upperwaterel[i])/(deepwaterel[i]-upperwaterel[i]);
			values[i]=deepwatermelt[i]*alpha+(1.-alpha)*upperwatermelt[i];
		}
		values[i] += perturbation[i];
	}

	this->AddInput(BasalforcingsFloatingiceMeltingRateEnum,&values[0],P1Enum);
}/*}}}*/
void       Element::MantlePlumeGeothermalFlux(){/*{{{*/

	const int NUM_VERTICES = this->GetNumberOfVertices();
	IssmDouble  mantleconductivity,nusselt,dtbg,plumeradius,topplumedepth,bottomplumedepth,plumex,plumey;
	IssmDouble  crustthickness,uppercrustthickness,uppercrustheat,lowercrustheat;
	IssmDouble  crustheat,plumeheat,dt,middleplumedepth,a,e,eprime,A0,lambda,Alambda,dAlambda;
	IssmDouble  x,y,z,c;
	IssmDouble  values[MAXVERTICES];
	IssmDouble *xyz_list = NULL;

	parameters->FindParam(&mantleconductivity,BasalforcingsMantleconductivityEnum);
	parameters->FindParam(&nusselt,BasalforcingsNusseltEnum);
	parameters->FindParam(&dtbg,BasalforcingsDtbgEnum);
	parameters->FindParam(&plumeradius,BasalforcingsPlumeradiusEnum);
	parameters->FindParam(&topplumedepth,BasalforcingsTopplumedepthEnum);
	parameters->FindParam(&bottomplumedepth,BasalforcingsBottomplumedepthEnum);
	parameters->FindParam(&plumex,BasalforcingsPlumexEnum);
	parameters->FindParam(&plumey,BasalforcingsPlumeyEnum);
	parameters->FindParam(&crustthickness,BasalforcingsCrustthicknessEnum);
	parameters->FindParam(&uppercrustthickness,BasalforcingsUppercrustthicknessEnum);
	parameters->FindParam(&uppercrustheat,BasalforcingsUppercrustheatEnum);
	parameters->FindParam(&lowercrustheat,BasalforcingsLowercrustheatEnum);

	this->GetVerticesCoordinates(&xyz_list);
	c=plumeradius;
	a=(bottomplumedepth-topplumedepth)/2.;
	e=pow(a*a-c*c,1./2.)/a;
	A0=(1-pow(e,2.))/pow(e,3.)*(1./2.*log((1+e)/(1-e))-e);
	for(int i=0;i<NUM_VERTICES;i++){
		y=xyz_list[i*3+0]-plumex;
		z=xyz_list[i*3+1]-plumey;
		x=-(a+topplumedepth+crustthickness);
		lambda=(-pow(a,2)-pow(c,2)+pow(x,2)+pow(y,2)+pow(z,2)+sqrt(pow(-pow(a,2)-pow(c,2)+pow(x,2)+pow(y,2)+pow(z,2),2)+4*(pow(c,2)*pow(x,2)+pow(a,2)*(-pow(c,2)+pow(y,2)+pow(z,2)))))/2;
		dAlambda=(-8*a*pow(c,2)*x*(-2*pow(a,2)+2*pow(c,2)+sqrt(2)*sqrt((a-c)*(a+c))*sqrt(pow(a,2)-pow(c,2)+pow(x,2)+pow(y,2)+pow(z,2)+sqrt(pow(-pow(a,2)-pow(c,2)+pow(x,2)+pow(y,2)+pow(z,2),2)+4*(pow(c,2)*pow(x,2)+pow(a,2)*(-pow(c,2)+pow(y,2)+pow(z,2))))))*(pow(a,4)*(pow(y,2)+pow(z,2))+pow(c,4)*(pow(y,2)+pow(z,2))+pow(pow(x,2)+pow(y,2)+pow(z,2),2)*(pow(x,2)+pow(y,2)+pow(z,2)+sqrt(pow(-pow(a,2)-pow(c,2)+pow(x,2)+pow(y,2)+pow(z,2),2)+4*(pow(c,2)*pow(x,2)+pow(a,2)*(-pow(c,2)+pow(y,2)+pow(z,2)))))+pow(c,2)*(pow(x,4)-pow(x,2)*(pow(y,2)+pow(z,2))-(pow(y,2)+pow(z,2))*(2*pow(y,2)+2*pow(z,2)+sqrt(pow(-pow(a,2)-pow(c,2)+pow(x,2)+pow(y,2)+pow(z,2),2)+4*(pow(c,2)*pow(x,2)+pow(a,2)*(-pow(c,2)+pow(y,2)+pow(z,2))))))+pow(a,2)*(-pow(x,4)+pow(x,2)*(pow(y,2)+pow(z,2))+(pow(y,2)+pow(z,2))*(-2*pow(c,2)+2*(pow(y,2)+pow(z,2))+sqrt(pow(-pow(a,2)-pow(c,2)+pow(x,2)+pow(y,2)+pow(z,2),2)+4*(pow(c,2)*pow(x,2)+pow(a,2)*(-pow(c,2)+pow(y,2)+pow(z,2))))))))/(sqrt((a-c)*(a+c))*sqrt(pow(-pow(a,2)-pow(c,2)+pow(x,2)+pow(y,2)+pow(z,2),2)+4*(pow(c,2)*pow(x,2)+pow(a,2)*(-pow(c,2)+pow(y,2)+pow(z,2))))*pow(pow(a,2)-pow(c,2)+pow(x,2)+pow(y,2)+pow(z,2)+sqrt(pow(-pow(a,2)-pow(c,2)+pow(x,2)+pow(y,2)+pow(z,2),2)+4*(pow(c,2)*pow(x,2)+pow(a,2)*(-pow(c,2)+pow(y,2)+pow(z,2)))),3.5)*pow(-(sqrt(2)*sqrt((a-c)*(a+c)))+sqrt(pow(a,2)-pow(c,2)+pow(x,2)+pow(y,2)+pow(z,2)+sqrt(pow(-pow(a,2)-pow(c,2)+pow(x,2)+pow(y,2)+pow(z,2),2)+4*(pow(c,2)*pow(x,2)+pow(a,2)*(-pow(c,2)+pow(y,2)+pow(z,2))))),2)*(sqrt(2)*sqrt((a-c)*(a+c))+sqrt(pow(a,2)-pow(c,2)+pow(x,2)+pow(y,2)+pow(z,2)+sqrt(pow(-pow(a,2)-pow(c,2)+pow(x,2)+pow(y,2)+pow(z,2),2)+4*(pow(c,2)*pow(x,2)+pow(a,2)*(-pow(c,2)+pow(y,2)+pow(z,2)))))));
		eprime=pow((a*a-plumeradius*plumeradius)/(a*a+lambda),1./2.);
		Alambda=(1.-e*e)/(e*e*e)*(1./2.*log((1.+eprime)/(1.-eprime))-eprime);
		dt=dtbg-(nusselt-1.)/(1.+A0*(nusselt-1.))*(Alambda*dtbg+x*dtbg*dAlambda);
		plumeheat=mantleconductivity*dt;
		crustheat=uppercrustheat*uppercrustthickness+lowercrustheat*(crustthickness-uppercrustthickness);
		values[i]=crustheat+plumeheat;
	}

	this->AddInput(BasalforcingsGeothermalfluxEnum,&values[0],P1Enum);
	xDelete<IssmDouble>(xyz_list);

}/*}}}*/
void       Element::MarshallElement2(MarshallHandle* marshallhandle,int numanalyses){/*{{{*/

	_assert_(this);
	if(marshallhandle->OperationNumber()==MARSHALLING_LOAD){
		this->nodes = NULL;
	}

	int object_enum = ElementEnum;
	marshallhandle->call(object_enum);

	marshallhandle->call(this->id);
	marshallhandle->call(this->sid);
	marshallhandle->call(this->lid);
	marshallhandle->call(this->element_type);
	marshallhandle->call(this->element_type_list,numanalyses);
}
/*}}}*/
void       Element::MigrateGroundingLine(IssmDouble* phi_ungrounding){/*{{{*/

	const int  NUM_VERTICES = this->GetNumberOfVertices();
	int        migration_style;
	IssmDouble bed_hydro,yts;
	IssmDouble melting[MAXVERTICES];
	IssmDouble phi[MAXVERTICES];
	IssmDouble h[MAXVERTICES];
	IssmDouble s[MAXVERTICES];
	IssmDouble b[MAXVERTICES];
	IssmDouble r[MAXVERTICES];
	IssmDouble sl[MAXVERTICES];

	/*Recover info at the vertices: */
	parameters->FindParam(&migration_style,GroundinglineMigrationEnum);
	parameters->FindParam(&yts,ConstantsYtsEnum);
	GetInputListOnVertices(&h[0],ThicknessEnum);
	GetInputListOnVertices(&s[0],SurfaceEnum);
	GetInputListOnVertices(&b[0],BaseEnum);
	GetInputListOnVertices(&r[0],BedEnum);
	GetInputListOnVertices(&sl[0],SealevelEnum);
	GetInputListOnVertices(&phi[0],MaskOceanLevelsetEnum);
	IssmDouble rho_water   = FindParam(MaterialsRhoSeawaterEnum);
	IssmDouble rho_ice     = FindParam(MaterialsRhoIceEnum);
	IssmDouble density     = rho_ice/rho_water;

	/*go through vertices, and update inputs, considering them to be TriaVertex type: */
	for(int i=0;i<NUM_VERTICES;i++){
		/* Contact FS*/
		if(migration_style == ContactEnum){
			phi[i]=phi_ungrounding[vertices[i]->Pid()];
			if(phi[i]>=0.) b[i]=r[i];
		}
		else if(migration_style == GroundingOnlyEnum && b[i]<r[i]){
			/*Ice shelf: if bed below bathymetry, impose it at the bathymetry and update surface, elso do nothing */
			b[i]=r[i];
		}
		else if(phi[i]<=0.){
			if(b[i]<=r[i]){
				b[i]        = r[i];
				s[i]        = b[i]+h[i];
			}
		}
		/*Ice sheet: if hydrostatic bed above bathymetry, ice sheet starts to unground, elso do nothing */
		/*Change only if AggressiveMigration or if the ice sheet is in contact with the ocean*/
		else{ // phi>0
			bed_hydro=-density*h[i]+sl[i];
			if (bed_hydro>r[i]){
				/*Unground only if the element is connected to the ice shelf*/
				if(migration_style==AggressiveMigrationEnum || migration_style==SubelementMigrationEnum){
					s[i] = (1.-density)*h[i]+sl[i];
					b[i] = -density*h[i]+sl[i];
				}
				else if(migration_style==SoftMigrationEnum && phi_ungrounding[vertices[i]->Pid()]<0.){
					s[i] = (1.-density)*h[i]+sl[i];
					b[i] = -density*h[i]+sl[i];
				}
				else{
					if(migration_style!=SoftMigrationEnum && migration_style!=ContactEnum && migration_style!=GroundingOnlyEnum) _error_("Error: migration should be Aggressive, Soft, Subelement, Contact or GroundingOnly");
				}
			}
		}
	}

	/*Recalculate phi*/
	for(int i=0;i<NUM_VERTICES;i++){
		if(migration_style==SoftMigrationEnum){
			bed_hydro=-density*h[i]+sl[i];
			if(phi[i]<0. || bed_hydro<=r[i] || phi_ungrounding[vertices[i]->Pid()]<0.){
				phi[i]=h[i]+(r[i]-sl[i])/density;
			}
		}
		else if(migration_style!=ContactEnum){
			phi[i]=h[i]+(r[i]-sl[i])/density;
		}
		else{
			/*do nothing*/
		}
	}
	this->AddInput(MaskOceanLevelsetEnum,&phi[0],P1Enum);

	/*Update inputs*/
	this->AddInput(SurfaceEnum,&s[0],P1Enum);
	this->AddInput(BaseEnum,&b[0],P1Enum);

}/*}}}*/
void       Element::MismipFloatingiceMeltingRate(){/*{{{*/

	IssmDouble thresholdthickness,upperdepthmelt;
	IssmDouble base[MAXVERTICES];
	IssmDouble bed[MAXVERTICES];
	IssmDouble meltratefactor[MAXVERTICES];
	IssmDouble values[MAXVERTICES];

	const int NUM_VERTICES = this->GetNumberOfVertices();

	parameters->FindParam(&thresholdthickness,BasalforcingsThresholdThicknessEnum);
	parameters->FindParam(&upperdepthmelt,BasalforcingsUpperdepthMeltEnum);
	this->GetInputListOnVertices(base,BaseEnum);
	this->GetInputListOnVertices(bed,BedEnum);
	this->GetInputListOnVertices(meltratefactor,BasalforcingsMeltrateFactorEnum);
	for(int i=0;i<NUM_VERTICES;i++){
		if(base[i]>upperdepthmelt){
			values[i]=0;
		}
		else{
			values[i]=meltratefactor[i]*tanh((base[i]-bed[i])/thresholdthickness)*(upperdepthmelt-base[i]);
		}
	}

	this->AddInput(BasalforcingsFloatingiceMeltingRateEnum,values,P1Enum);

}/*}}}*/
void       Element::MonthlyFactorBasin(IssmDouble* monthlyfac,int enum_type){/*{{{*/

	/*Variable declaration*/
	bool ratevariable;
   const int numvertices = this->GetNumberOfVertices();
	int basinid,mindex,mindexnext,basinenum_type,varenum_type,indperiod;
   IssmDouble time,dt,fracyear,fracyearnext,fracmonth,fracmonthnext,yts; 
   IssmDouble monthsteps[12]  = {0.,1./12,2./12,3./12,4./12,5./12,6./12,7./12,8./12,9./12,10./12,11./12};
   IssmDouble* monthlyfac_b   = xNew<IssmDouble>(12);
   IssmDouble* monthlyrate_b  = xNew<IssmDouble>(12);
	IssmDouble* fracdtinmonth  = xNew<IssmDouble>(12);
	IssmDouble* rateinmonth    = xNew<IssmDouble>(numvertices*12);
	IssmDouble* varlistinput   = xNew<IssmDouble>(numvertices);
	IssmDouble* varlist        = xNewZeroInit<IssmDouble>(numvertices);

	/*Get field-specific enums*/
   switch(enum_type){
      case(FrontalForcingsSubglacialDischargearmaEnum):
         basinenum_type = FrontalForcingsBasinIdEnum;
         varenum_type   = FrontalForcingsSubglacialDischargeEnum;
         ratevariable   = true;
			break;
		case(HydrologyarmapwEnum):
         basinenum_type = HydrologyBasinsIdEnum;
         varenum_type   = FrictionWaterPressureEnum;
         ratevariable   = false;
			break;
	}

	/*Evaluate the month index now and at (now-timestepjump)*/
	this->parameters->FindParam(&yts,ConstantsYtsEnum);
	this->parameters->FindParam(&time,TimeEnum);
   this->parameters->FindParam(&dt,TimesteppingTimeStepEnum); _assert_(dt<yts);
	fracyear     = time/yts-floor(time/yts);
	fracyearnext = (time+dt)/yts-floor((time+dt)/yts);
	for(int i=1;i<12;i++){
		if(fracyear>=monthsteps[i-1])     mindex     = i-1;
		if(fracyearnext>=monthsteps[i-1]) mindexnext = i-1;
	}
	if(fracyear>=monthsteps[11])         mindex     = 11;
	if(fracyearnext>=monthsteps[11])     mindexnext = 11;

	/*Calculate fraction of the time step spent in each month*/
	for(int i=0;i<12;i++){
		if(mindex<i && mindexnext>i)                            fracdtinmonth[i] = 1.0/dt*yts/12.0;
		else if(mindex<i && mindexnext<i && mindexnext<mindex)  fracdtinmonth[i] = 1.0/dt*yts/12.0;
		else if(mindex>i && mindexnext<mindex && mindexnext>i)  fracdtinmonth[i] = 1.0/dt*yts/12.0;
		else if(mindex>i && mindexnext<mindex && mindexnext==i) fracdtinmonth[i] = 1.0/dt*yts*(fracyearnext-monthsteps[i]);
		else if(mindex==i && mindexnext==i)                     fracdtinmonth[i] = 1.0/dt*yts*(fracyearnext-fracyear); 
		else if(mindex==i && mindexnext!=mindex)                fracdtinmonth[i] = 1.0/dt*yts*(1.0/12-(fracyear-monthsteps[i]));
		else if(mindexnext==i && mindex!=mindexnext)            fracdtinmonth[i] = 1.0/dt*yts*(fracyearnext-monthsteps[i]);
		else	                                                  fracdtinmonth[i] = 0.0;
	}

	/*Get basin-specific parameters of the element*/
   this->GetInputValue(&basinid,basinenum_type);
	for(int i=0;i<12;i++) monthlyfac_b[i]   = monthlyfac[basinid*12+i];

	/*Retrieve input*/
	this->GetInputListOnVertices(varlistinput,varenum_type);

	/*Calculate monthly rate for each month and weight-average it for application over dt*/
	for(int v=0;v<numvertices;v++){
		for(int i=0;i<12;i++){
			if(ratevariable){
				rateinmonth[v*12+i] = varlistinput[v]*monthlyfac_b[i]*12;
				varlist[v]          = varlist[v]+fracdtinmonth[i]*rateinmonth[v*12+i];
			}
			else varlist[v]       = varlist[v]+fracdtinmonth[i]*monthlyfac_b[i]*varlistinput[v];
		}
	}
	/*Update input*/
   this->AddInput(varenum_type,varlist,P1DGEnum);

	/*Clean-up*/
	xDelete<IssmDouble>(fracdtinmonth);
	xDelete<IssmDouble>(rateinmonth);
	xDelete<IssmDouble>(monthlyfac_b);
	xDelete<IssmDouble>(monthlyrate_b);
	xDelete<IssmDouble>(varlist);
	xDelete<IssmDouble>(varlistinput);
}/*}}}*/
void       Element::MonthlyPiecewiseLinearEffectBasin(int nummonthbreaks,IssmDouble* monthlyintercepts,IssmDouble* monthlytrends,IssmDouble* monthlydatebreaks,int enum_type){/*{{{*/

	/*Variable declaration*/
   const int numvertices = this->GetNumberOfVertices();
	int basinid,mindex,basinenum_type,varenum_type,indperiod;
	int numperiods = nummonthbreaks+1;
   IssmDouble time,fracyear,yts; 
   IssmDouble monthsteps[12] = {0.,1./12,2./12,3./12,4./12,5./12,6./12,7./12,8./12,9./12,10./12,11./12};
   IssmDouble* datebreaks_b  = xNew<IssmDouble>(nummonthbreaks);
	IssmDouble* intercepts_b  = xNew<IssmDouble>(numperiods*12);
	IssmDouble* trends_b      = xNew<IssmDouble>(numperiods*12);
	IssmDouble* varlist       = xNew<IssmDouble>(numvertices);

	/*Get field-specific enums*/
   switch(enum_type){
      case(FrontalForcingsRignotarmaEnum):
         basinenum_type = FrontalForcingsBasinIdEnum;
         varenum_type   = ThermalForcingEnum;
         break;
	}

	/*Evaluate the month index*/
	this->parameters->FindParam(&yts,ConstantsYtsEnum);
	this->parameters->FindParam(&time,TimeEnum);
	fracyear = time/yts-floor(time/yts);
	for(int i=1;i<12;i++){
		if(fracyear>=monthsteps[i-1]) mindex = i-1;
	}
	if(fracyear>=monthsteps[11])     mindex = 11;

	/*Get basin-specific parameters of the element*/
   this->GetInputValue(&basinid,basinenum_type);
	if(nummonthbreaks>0){
      for(int i=0;i<nummonthbreaks;i++) datebreaks_b[i] = monthlydatebreaks[basinid*nummonthbreaks+i];
   }
	for(int i=0;i<numperiods;i++){
		intercepts_b[i]  = monthlyintercepts[basinid*12*numperiods+mindex+12*i];
		trends_b[i]      = monthlytrends[basinid*12*numperiods+mindex+12*i];
	}

	/*Compute piecewise-linear function*/
	IssmDouble telapsed_break,piecewiselin;
	if(nummonthbreaks>0){
		/*Find index of time compared to the breakpoints*/
		indperiod = 0;
		for(int i=0;i<nummonthbreaks;i++){
			if(time>datebreaks_b[i]) indperiod = i+1;
		}
		/*Compute intercept+trend terms with parameters of indperiod*/
      if(indperiod==0) telapsed_break = time;
      else             telapsed_break = time-datebreaks_b[indperiod-1];
      piecewiselin = intercepts_b[indperiod]+trends_b[indperiod]*telapsed_break;
	}
   else piecewiselin = intercepts_b[indperiod]+trends_b[indperiod]*time;

	/*Retrieve values non-adjusted for monthly effects*/
   this->GetInputListOnVertices(varlist,varenum_type);

	/*Adjust values using monthly effects*/
	for(int v=0;v<numvertices;v++) varlist[v] = varlist[v]+piecewiselin;

	/*Add input to element*/
   this->AddInput(varenum_type,varlist,P1Enum);

	/*Clean-up*/
   xDelete<IssmDouble>(datebreaks_b);
   xDelete<IssmDouble>(intercepts_b);
   xDelete<IssmDouble>(trends_b);
	xDelete<IssmDouble>(varlist);
}
/*}}}*/
void       Element::BeckmannGoosseFloatingiceMeltingRate(){/*{{{*/

	bool isthermalforcing;       
	int numvertices                 = this->GetNumberOfVertices();
	IssmDouble T_fr,T_forcing,ocean_heat_flux;
	/*Material properties*/
	IssmDouble rho_water            = this->FindParam(MaterialsRhoSeawaterEnum);
	IssmDouble rho_ice              = this->FindParam(MaterialsRhoIceEnum);
	IssmDouble latentheat           = this->FindParam(MaterialsLatentheatEnum);
	IssmDouble mixed_layer_capacity = this->FindParam(MaterialsMixedLayerCapacityEnum);
	IssmDouble thermal_exchange_vel = this->FindParam(MaterialsThermalExchangeVelocityEnum);

	IssmDouble base[MAXVERTICES];
	IssmDouble values[MAXVERTICES];
	IssmDouble oceansalinity [MAXVERTICES];
	IssmDouble oceantemp[MAXVERTICES];
	IssmDouble oceanthermalforcing[MAXVERTICES];
	IssmDouble meltratefactor[MAXVERTICES];

	/*Determine if we use temperature-and-salinity or thermal forcing*/
	this->parameters->FindParam(&isthermalforcing,BasalforcingsIsThermalForcingEnum);
	/*Retrieve Inputs*/
	this->GetInputListOnVertices(base,BaseEnum);
	this->GetInputListOnVertices(meltratefactor,BasalforcingsMeltrateFactorEnum);
	if(isthermalforcing==false){
		this->GetInputListOnVertices(oceansalinity,BasalforcingsOceanSalinityEnum);
		this->GetInputListOnVertices(oceantemp,BasalforcingsOceanTempEnum);
	}
	else{
		this->GetInputListOnVertices(oceanthermalforcing,ThermalForcingEnum);
	}

	Gauss* gauss=this->NewGauss();
	for(int i=0;i<numvertices;i++){
		/*Case of temperature-and-salinity inputs*/
		if(isthermalforcing==false){
			T_fr       = (0.0939-0.057*oceansalinity[i]+7.64e-4*base[i]); //freezing point [degC]
			T_forcing  = oceantemp[i]-T_fr; //thermal forcing [K]
		}
		/*Case of thermal forcing input*/
		else{
			T_forcing  = oceanthermalforcing[i]; 
		}
		/*Compute heat flux from (Beckmann and Goosse, 2003) (>0 means heat flux from ocean to ice)*/
		ocean_heat_flux = meltratefactor[i]*rho_water*mixed_layer_capacity*thermal_exchange_vel*T_forcing; // ocean-to-ice heat flux [W/m^2]
		/*Save melt values [m/s]*/
		values[i]       = ocean_heat_flux/(latentheat*rho_ice);
	}
	this->AddInput(BasalforcingsFloatingiceMeltingRateEnum,values,P1Enum);
	delete gauss;
}/*}}}*/
void       Element::MungsmtpParameterization(void){/*{{{*/
	/*Are we on the base? If not, return*/
	if(!IsOnBase()) return;

	const int NUM_VERTICES 					= this->GetNumberOfVertices();
	const int NUM_VERTICES_MONTHS_PER_YEAR 	= NUM_VERTICES * 12;

	int i;
	int*        vertexlids=xNew<int>(NUM_VERTICES);
	IssmDouble* monthlytemperatures=xNew<IssmDouble>(NUM_VERTICES_MONTHS_PER_YEAR);
	IssmDouble* monthlyprec=xNew<IssmDouble>(NUM_VERTICES_MONTHS_PER_YEAR);
	IssmDouble* TemperaturesPresentday=xNew<IssmDouble>(NUM_VERTICES_MONTHS_PER_YEAR);
	IssmDouble* TemperaturesLgm=xNew<IssmDouble>(NUM_VERTICES_MONTHS_PER_YEAR);
	IssmDouble* PrecipitationsPresentday=xNew<IssmDouble>(NUM_VERTICES_MONTHS_PER_YEAR);
	IssmDouble* PrecipitationsLgm=xNew<IssmDouble>(NUM_VERTICES_MONTHS_PER_YEAR);
	IssmDouble* tmp=xNew<IssmDouble>(NUM_VERTICES);
	IssmDouble TdiffTime,PfacTime;

	/*Recover parameters*/
	IssmDouble time,yts,time_yr;
	this->parameters->FindParam(&time,TimeEnum);
	this->parameters->FindParam(&yts,ConstantsYtsEnum);
	this->GetVerticesLidList(vertexlids);
	time_yr=floor(time/yts)*yts;

	/*Recover present day temperature and precipitation*/
	DatasetInput* dinput1=this->GetDatasetInput(SmbTemperaturesPresentdayEnum);   _assert_(dinput1);
	DatasetInput* dinput2=this->GetDatasetInput(SmbTemperaturesLgmEnum);          _assert_(dinput2);
	DatasetInput* dinput3=this->GetDatasetInput(SmbPrecipitationsPresentdayEnum); _assert_(dinput3);
	DatasetInput* dinput4=this->GetDatasetInput(SmbPrecipitationsLgmEnum);        _assert_(dinput4);

	/*loop over vertices: */
	Gauss* gauss=this->NewGauss();
	for(int month=0;month<12;month++) {
		for(int iv=0;iv<NUM_VERTICES;iv++) {
			gauss->GaussVertex(iv);
			dinput1->GetInputValue(&TemperaturesPresentday[iv*12+month],gauss,month);
			dinput2->GetInputValue(&TemperaturesLgm[iv*12+month],gauss,month);
			dinput3->GetInputValue(&PrecipitationsPresentday[iv*12+month],gauss,month);
			dinput4->GetInputValue(&PrecipitationsLgm[iv*12+month],gauss,month);

			PrecipitationsPresentday[iv*12+month]=PrecipitationsPresentday[iv*12+month]*yts;
			PrecipitationsLgm[iv*12+month]=PrecipitationsLgm[iv*12+month]*yts;
		}
	}

	/*Recover interpolation parameters at time t*/
	this->parameters->FindParam(&TdiffTime,SmbTdiffEnum,time);
	this->parameters->FindParam(&PfacTime,SmbPfacEnum,time);

	/*Compute the temperature and precipitation*/
	for(int iv=0;iv<NUM_VERTICES;iv++){
		ComputeMungsmTemperaturePrecipitation(TdiffTime,PfacTime,
					&PrecipitationsLgm[iv*12],&PrecipitationsPresentday[iv*12],
					&TemperaturesLgm[iv*12], &TemperaturesPresentday[iv*12],
					&monthlytemperatures[iv*12], &monthlyprec[iv*12]);
	}

	/*Update inputs*/
	for (int imonth=0;imonth<12;imonth++) {
		for(i=0;i<NUM_VERTICES;i++) tmp[i]=monthlytemperatures[i*12+imonth];
		switch(this->ObjectEnum()){
			case TriaEnum:  this->inputs->SetTriaDatasetInput(SmbMonthlytemperaturesEnum,imonth,P1Enum,NUM_VERTICES,vertexlids,tmp); break;
			case PentaEnum: this->inputs->SetPentaDatasetInput(SmbMonthlytemperaturesEnum,imonth,P1Enum,NUM_VERTICES,vertexlids,tmp); break;
			default: _error_("Not implemented yet");
		}
		for(i=0;i<NUM_VERTICES;i++) tmp[i]=monthlyprec[i*12+imonth]/yts;
		switch(this->ObjectEnum()){
			case TriaEnum:  this->inputs->SetTriaDatasetInput(SmbPrecipitationEnum,imonth,P1Enum,NUM_VERTICES,vertexlids,tmp); break;
			case PentaEnum: this->inputs->SetPentaDatasetInput(SmbPrecipitationEnum,imonth,P1Enum,NUM_VERTICES,vertexlids,tmp); break;
			default: _error_("Not implemented yet");
		}
	}

	switch(this->ObjectEnum()){
		case TriaEnum: break;
		case PentaEnum:
		case TetraEnum:
							this->DatasetInputExtrude(SmbMonthlytemperaturesEnum,-1);
							this->DatasetInputExtrude(SmbPrecipitationEnum,-1);
							break;
		default: _error_("Not implemented yet");
	}

	/*clean-up*/
	delete gauss;
	xDelete<IssmDouble>(monthlytemperatures);
	xDelete<IssmDouble>(monthlyprec);
	xDelete<IssmDouble>(TemperaturesPresentday);
	xDelete<IssmDouble>(TemperaturesLgm);
	xDelete<IssmDouble>(PrecipitationsPresentday);
	xDelete<IssmDouble>(PrecipitationsLgm);
	xDelete<IssmDouble>(tmp);
	xDelete<int>(vertexlids);

}
/*}}}*/
ElementMatrix* Element::NewElementMatrix(int approximation_enum){/*{{{*/
	return new ElementMatrix(nodes,this->GetNumberOfNodes(),this->parameters,approximation_enum);
}
/*}}}*/
ElementMatrix* Element::NewElementMatrixCoupling(int number_nodes,int approximation_enum){/*{{{*/
	return new ElementMatrix(nodes,number_nodes,this->parameters,approximation_enum);
}
/*}}}*/
ElementVector* Element::NewElementVector(int approximation_enum){/*{{{*/
	return new ElementVector(nodes,this->GetNumberOfNodes(),this->parameters,approximation_enum);
}
/*}}}*/
void       Element::PicoUpdateBoxid(int* max_boxid_basin_list){/*{{{*/

	if(!this->IsIceInElement() || !this->IsAllFloating()) return;

	int        basin_id;
	IssmDouble dist_gl,dist_cf;

	this->GetInputValue(&basin_id,BasalforcingsPicoBasinIdEnum);
	IssmDouble boxid_max=reCast<IssmDouble>(max_boxid_basin_list[basin_id])+1.;

	Input* dist_gl_input=this->GetInput(DistanceToGroundinglineEnum); _assert_(dist_gl_input);
	Input* dist_cf_input=this->GetInput(DistanceToCalvingfrontEnum);  _assert_(dist_cf_input);

	/*Get dist_gl and dist_cf at center of element*/
	Gauss* gauss=this->NewGauss(1); gauss->GaussPoint(0);
	dist_gl_input->GetInputValue(&dist_gl,gauss);
	dist_cf_input->GetInputValue(&dist_cf,gauss);
	delete gauss;

	/*Ensure values are positive for floating ice*/
	dist_gl = fabs(dist_gl);
	dist_cf = fabs(dist_cf);

	/*Compute relative distance to grounding line*/
	IssmDouble rel_dist_gl=dist_gl/(dist_gl+dist_cf);

	/*Assign box numbers based on rel_dist_gl*/
	int boxid = -1;
	for(IssmDouble i=0.;i<boxid_max;i++){
		IssmDouble lowbound  = 1. -sqrt((boxid_max-i   )/boxid_max);
		IssmDouble highbound = 1. -sqrt((boxid_max-i-1.)/boxid_max);
		if(rel_dist_gl>=lowbound && rel_dist_gl<=highbound){
			boxid=reCast<int>(i);
			break;
		}
	}
	if(boxid==-1) _error_("No boxid found for element " << this->Sid() << "!");

	this->SetIntInput(this->inputs,BasalforcingsPicoBoxIdEnum, boxid);

}/*}}}*/
void       Element::PicoUpdateBox(int loopboxid){/*{{{*/

	if(!this->IsIceInElement() || !this->IsAllFloating()) return;

	int boxid;
	this->GetInputValue(&boxid,BasalforcingsPicoBoxIdEnum);
	if(loopboxid!=boxid) return;

	const int NUM_VERTICES = this->GetNumberOfVertices();
	int        basinid, maxbox, num_basins, numnodes, M;
	IssmDouble gamma_T, overturning_coeff, thickness;
	IssmDouble pressure, T_star,p_coeff, q_coeff;
	bool       isplume;

	/*Get variables*/
	IssmDouble rhoi       = this->FindParam(MaterialsRhoIceEnum);
	IssmDouble rhow       = this->FindParam(MaterialsRhoSeawaterEnum);
	IssmDouble earth_grav = this->FindParam(ConstantsGEnum);
	IssmDouble rho_star   = 1033.;             // kg/m^3
	IssmDouble nu         = rhoi/rhow;
	IssmDouble latentheat = this->FindParam(MaterialsLatentheatEnum);
	IssmDouble c_p_ocean  = this->FindParam(MaterialsMixedLayerCapacityEnum);
	IssmDouble lambda     = latentheat/c_p_ocean;
	IssmDouble a          = -0.0572;          // K/psu
	IssmDouble b          = 0.0788 + this->FindParam(MaterialsMeltingpointEnum);  //K
	IssmDouble c          = 7.77e-4;
	IssmDouble alpha      = 7.5e-5;           // 1/K
	IssmDouble Beta       = 7.7e-4;           // K

	/* Get non-box-specific parameters and inputs */
	this->parameters->FindParam(&num_basins, BasalforcingsPicoNumBasinsEnum);
	this->parameters->FindParam(&gamma_T,BasalforcingsPicoGammaTEnum);
	this->parameters->FindParam(&maxbox,BasalforcingsPicoMaxboxcountEnum);
	this->GetInputValue(&basinid,BasalforcingsPicoBasinIdEnum);
	this->parameters->FindParam(&isplume, BasalforcingsPicoIsplumeEnum);
	Input *thickness_input = this->GetInput(ThicknessEnum); 
   _assert_(basinid<=num_basins);
   _assert_(thickness_input);

	IssmDouble* boxareas = NULL;
	this->parameters->FindParam(&boxareas,&M,BasalforcingsPicoBoxAreaEnum);
	_assert_(M==num_basins*maxbox);

	IssmDouble area_boxi        = boxareas[basinid*maxbox+boxid];
	IssmDouble g1               = area_boxi*gamma_T;

	IssmDouble basalmeltrates_shelf[MAXVERTICES];
	IssmDouble potential_pressure_melting_point[MAXVERTICES];
	IssmDouble Tocs[MAXVERTICES];
	IssmDouble Socs[MAXVERTICES];

	/* First box calculations */
	if(boxid==0){
		/* Get box1 parameters and inputs */
		IssmDouble time, toc_farocean, soc_farocean;
      IssmDouble overturnings[MAXVERTICES];
		this->parameters->FindParam(&time,TimeEnum);
		this->parameters->FindParam(&toc_farocean, basinid, time, BasalforcingsPicoFarOceantemperatureEnum);
		this->parameters->FindParam(&soc_farocean, basinid, time, BasalforcingsPicoFarOceansalinityEnum);
		IssmDouble s1 = soc_farocean/(nu*lambda);
		Input *overturningC_input = this->GetInput(BasalforcingsPicoOverturningCoeffEnum); _assert_(overturningC_input);

		/* Start looping on the number of verticies and calculate ocean vars */
		Gauss* gauss=this->NewGauss();
		for(int i=0;i<NUM_VERTICES;i++){
			gauss->GaussVertex(i);
			thickness_input->GetInputValue(&thickness,gauss);
			overturningC_input->GetInputValue(&overturning_coeff,gauss);
			pressure = (rhoi*earth_grav*1e-4)*thickness;
			T_star   = a*soc_farocean+b-c*pressure-toc_farocean;
			p_coeff  = g1/(overturning_coeff*rho_star*(Beta*s1-alpha));
			q_coeff  = T_star*(g1/(overturning_coeff*rho_star*(Beta*s1-alpha)));

			/* To avoid negatives under the square root */
			if((0.25*pow(p_coeff,2)-q_coeff)<0) q_coeff = 0.25*p_coeff*p_coeff;

			Tocs[i] = toc_farocean-(-0.5*p_coeff+sqrt(0.25*pow(p_coeff,2)-q_coeff));
			Socs[i] = soc_farocean-(soc_farocean/(nu*lambda))*(toc_farocean-Tocs[i]);
			potential_pressure_melting_point[i] = a*Socs[i]+b-c*pressure;
			if(!isplume) basalmeltrates_shelf[i] = (-gamma_T/(nu*lambda))*(potential_pressure_melting_point[i]-Tocs[i]);
			overturnings[i] = overturning_coeff*rho_star*(Beta*(soc_farocean-Socs[i])-alpha*(toc_farocean-Tocs[i]));
		}

		if(!isplume) this->AddInput(BasalforcingsFloatingiceMeltingRateEnum,&basalmeltrates_shelf[0],P1DGEnum);
		this->AddInput(BasalforcingsPicoSubShelfOceanTempEnum,&Tocs[0],P1DGEnum);
		this->AddInput(BasalforcingsPicoSubShelfOceanSalinityEnum,&Socs[0],P1DGEnum);
		this->AddInput(BasalforcingsPicoSubShelfOceanOverturningEnum,&overturnings[0],P1DGEnum);

		/*Cleanup and return*/
		delete gauss;
	}

	/* Subsequent box calculations */
	else {
		/* Get subsequent box parameters and inputs */
		IssmDouble* toc_weighted_avg         = NULL;
		IssmDouble* soc_weighted_avg         = NULL;
		IssmDouble* overturning_weighted_avg = NULL;
		this->parameters->FindParam(&toc_weighted_avg,&num_basins,BasalforcingsPicoAverageTemperatureEnum);
		this->parameters->FindParam(&soc_weighted_avg,&num_basins,BasalforcingsPicoAverageSalinityEnum);
		this->parameters->FindParam(&overturning_weighted_avg,&num_basins,BasalforcingsPicoAverageOverturningEnum);
		IssmDouble mean_toc                  = toc_weighted_avg[basinid];
		IssmDouble mean_soc                  = soc_weighted_avg[basinid];
		IssmDouble mean_overturning          = overturning_weighted_avg[basinid];
		IssmDouble g2                        = g1/(nu*lambda);

		/* Start looping on the number of verticies and calculate ocean vars */
		Gauss* gauss=this->NewGauss();
		for(int i=0;i<NUM_VERTICES;i++){
			gauss->GaussVertex(i);
			thickness_input->GetInputValue(&thickness,gauss);
			pressure = (rhoi*earth_grav*1.e-4)*thickness;
			T_star   = a*mean_soc+b-c*pressure-mean_toc;
			Tocs[i]  = mean_toc+T_star*(g1/(mean_overturning+g1-g2*a*mean_soc));
			Socs[i]  = mean_soc-mean_soc*((mean_toc-Tocs[i])/(nu*lambda));
			potential_pressure_melting_point[i] = a*Socs[i]+b-c*pressure;
			if(!isplume) basalmeltrates_shelf[i] = (-gamma_T/(nu*lambda))*(potential_pressure_melting_point[i]-Tocs[i]);
		}

		if(!isplume) this->AddInput(BasalforcingsFloatingiceMeltingRateEnum,basalmeltrates_shelf,P1DGEnum);
		this->AddInput(BasalforcingsPicoSubShelfOceanTempEnum,Tocs,P1DGEnum);
		this->AddInput(BasalforcingsPicoSubShelfOceanSalinityEnum,Socs,P1DGEnum);

		/*Cleanup and return*/
		xDelete<IssmDouble>(toc_weighted_avg);
		xDelete<IssmDouble>(soc_weighted_avg);
		xDelete<IssmDouble>(overturning_weighted_avg);
		delete gauss;
	}

	/*Cleanup and return*/
	xDelete<IssmDouble>(boxareas);

}/*}}}*/
void       Element::PicoComputeBasalMelt(){/*{{{*/

	if(!this->IsIceInElement() || !this->IsAllFloating()) return;

	const int NUM_VERTICES = this->GetNumberOfVertices();

	IssmDouble E0, Cd, CdT, YT, lam1, lam2, lam3, M0, CdTS0, y1, y2, x0;
	IssmDouble Tf_gl, YTS, CdTS, G1, G2, G3, g_alpha, M, l, X_hat, M_hat;
	IssmDouble alpha, zgl, Toc, Soc, z_base, yts, slopex, slopey;

	/*Get variables*/
	E0    = 3.6e-2;        //Entrainment coefficient
	Cd    = 2.5e-3;        //Drag coefficient
	CdT   = 1.1e-3;        //Turbulent heat exchange coefficient
	YT    = CdT/sqrt(Cd);  //Heat exchange coefficient
	lam1  = -5.73e-2;      //Freezing point-salinity coefficient (degrees C)
	lam2  = 8.32e-2;       //Freezing point offset (degrees C)
	lam3  = 7.61e-4;       //Freezing point-depth coefficient (K m-1)
	M0    = 10.;           //Melt-rate parameter (m yr-1 C-2)
	CdTS0 = 6e-4;          //Heat exchange parameter
	y1    = 0.545;         //Heat exchange parameter
	y2    = 3.5e-5;        //Heat exchange parameter
	x0    = 0.56;          //Dimentionless scaling factor

	/*Define arrays*/
	IssmDouble basalmeltrates_shelf[MAXVERTICES];

	/*Polynomial coefficients*/
	IssmDouble p[12];
	p[0]  =  0.1371330075095435;
	p[1]  =  5.527656234709359E1;
	p[2]  = -8.951812433987858E2;
	p[3]  =  8.927093637594877E3;
	p[4]  = -5.563863123811898E4;
	p[5]  =  2.218596970948727E5;
	p[6]  = -5.820015295669482E5;
	p[7]  =  1.015475347943186E6;
	p[8]  = -1.166290429178556E6;
	p[9]  =  8.466870335320488E5;
	p[10] = -3.520598035764990E5;
	p[11] =  6.387953795485420E4;

	/*Get inputs*/
	Input* zgl_input         = this->GetInput(GroundinglineHeightEnum);                     _assert_(zgl_input);
	Input* toc_input         = this->GetInput(BasalforcingsPicoSubShelfOceanTempEnum);      _assert_(toc_input);
	Input* soc_input         = this->GetInput(BasalforcingsPicoSubShelfOceanSalinityEnum);  _assert_(soc_input);
	Input* base_input        = this->GetInput(BaseEnum);                                    _assert_(base_input);
	Input* baseslopex_input  = this->GetInput(BaseSlopeXEnum);                              _assert_(baseslopex_input);
	Input* baseslopey_input  = this->GetInput(BaseSlopeYEnum);                              _assert_(baseslopey_input);
	this->FindParam(&yts, ConstantsYtsEnum);

	/*Loop over the number of vertices in this element*/
	Gauss* gauss=this->NewGauss();
	for(int i=0;i<NUM_VERTICES;i++){
		gauss->GaussVertex(i);

		/*Get inputs*/
		zgl_input->GetInputValue(&zgl,gauss);
		toc_input->GetInputValue(&Toc,gauss); //(K)
		soc_input->GetInputValue(&Soc,gauss);
		base_input->GetInputValue(&z_base,gauss);
		baseslopex_input->GetInputValue(&slopex,gauss);
		baseslopey_input->GetInputValue(&slopey,gauss);

		/*Compute ice shelf base slope (radians)*/
		alpha = atan(sqrt(slopex*slopex + slopey*slopey));
		if(alpha>=M_PI) alpha = M_PI - 0.001;               //ensure sin(alpha) > 0 for meltrate calculations

		/*Make necessary conversions*/
		Toc = Toc-273.15;
		if(zgl>z_base) zgl=z_base;

		/*Low bound for Toc to ensure X_hat is between 0 and 1*/
		if(Toc<lam1*Soc+lam2) Toc=lam1*Soc+lam2;

		/*Compute parameters needed for melt-rate calculation*/
		Tf_gl = lam1*Soc+lam2+lam3*zgl;                                              //Characteristic freezing point
		YTS = YT*(y1+y2*(((Toc-Tf_gl)*E0*sin(alpha))/(lam3*(CdTS0+E0*sin(alpha))))); //Effective heat exchange coefficient
		CdTS = sqrt(Cd)*YTS;                                                         //Heat exchange coefficient
		G1 = sqrt(sin(alpha)/(Cd+E0*sin(alpha)));                                    //Geometric factor
		G2 = sqrt(CdTS/(CdTS+E0*sin(alpha)));                                        //Geometric factor
		G3 = (E0*sin(alpha))/(CdTS+E0*sin(alpha));                                   //Geometric factor
		g_alpha = G1*G2*G3;                                                          //Melt scaling factor
		M = M0*g_alpha*pow((Toc-Tf_gl),2);                                           //Melt-rate scale
		l = ((Toc-Tf_gl)*(x0*CdTS+E0*sin(alpha)))/(lam3*x0*(CdTS+E0*sin(alpha)));    //Length scale
		X_hat = (z_base-zgl)/l;                                                      //Dimentionless coordinate system

		/*Compute polynomial fit*/
		M_hat = 0.;                                                                  //Reset summation variable for each node
		for(int ii=0;ii<12;ii++) {
			M_hat += p[ii]*pow(X_hat,ii);                                               //Polynomial fit
		}

		/*Compute melt-rate*/
		basalmeltrates_shelf[i] = (M*M_hat)/yts;                                     //Basal melt-rate (m/s)
	}

	/*Save computed melt-rate*/
	this->AddInput(BasalforcingsFloatingiceMeltingRateEnum,&basalmeltrates_shelf[0],P1DGEnum);

	/*Cleanup and return*/
	delete gauss;
}/*}}}*/
void       Element::PositiveDegreeDay(IssmDouble* pdds,IssmDouble* pds,IssmDouble signorm,bool ismungsm,bool issetpddfac){/*{{{*/

	const int NUM_VERTICES 		= this->GetNumberOfVertices();
	const int NUM_VERTICES_MONTHS_PER_YEAR = NUM_VERTICES * 12;

	int  		i,vertexlids[MAXVERTICES];
	bool     isinitialized=false;
	IssmDouble* agd=xNew<IssmDouble>(NUM_VERTICES); // surface mass balance
	IssmDouble* melt=xNew<IssmDouble>(NUM_VERTICES); // surface melt
	IssmDouble* accu=xNew<IssmDouble>(NUM_VERTICES); // surface precipitation
	IssmDouble* monthlytemperatures=xNew<IssmDouble>(NUM_VERTICES_MONTHS_PER_YEAR);
	IssmDouble* monthlyprec=xNew<IssmDouble>(NUM_VERTICES_MONTHS_PER_YEAR);
	IssmDouble* yearlytemperatures=xNew<IssmDouble>(NUM_VERTICES); memset(yearlytemperatures, 0., NUM_VERTICES*sizeof(IssmDouble));
	IssmDouble* tmp=xNew<IssmDouble>(NUM_VERTICES);
	IssmDouble* h=xNew<IssmDouble>(NUM_VERTICES);
	IssmDouble* s=xNew<IssmDouble>(NUM_VERTICES);
	IssmDouble* s0p=xNew<IssmDouble>(NUM_VERTICES);
	IssmDouble* s0t=xNew<IssmDouble>(NUM_VERTICES);
	IssmDouble pddsnowfac = -1.;
	IssmDouble pddicefac = -1.;
	IssmDouble accsumM = 0.0;
	IssmDouble accsumSMB = 0.0;
	IssmDouble accsumP = 0.0;
	IssmDouble avgM = 0.0;
	IssmDouble avgSMB = 0.0;
	IssmDouble avgP = 0.0;
	IssmDouble rho_water,rho_ice,desfac,rlaps,rlapslgm;
	IssmDouble PfacTime,TdiffTime,sealevTime;
	IssmDouble mavg=1./12.; //factor for monthly average

	/*Get vertex Lids for later*/
	this->GetVerticesLidList(&vertexlids[0]);

	/*Check for smb initialization*/
	this->GetInputValue(&isinitialized,SmbIsInitializedEnum);

	/*Get material parameters :*/
	rho_water=this->FindParam(MaterialsRhoSeawaterEnum);
	rho_ice=this->FindParam(MaterialsRhoIceEnum);

	/*Get some pdd parameters*/
	desfac=this->FindParam(SmbDesfacEnum);
	rlaps=this->FindParam(SmbRlapsEnum);
	rlapslgm=this->FindParam(SmbRlapslgmEnum);

	IssmDouble time,yts,time_yr,dt;
	this->parameters->FindParam(&time,TimeEnum);
	this->parameters->FindParam(&yts,ConstantsYtsEnum);
	this->parameters->FindParam(&dt,TimesteppingTimeStepEnum);          /*transient core time step*/
	time_yr=floor(time/yts)*yts;

	/*Get inputs*/
	DatasetInput* dinput =this->GetDatasetInput(SmbMonthlytemperaturesEnum); _assert_(dinput);
	DatasetInput* dinput2=this->GetDatasetInput(SmbPrecipitationEnum);       _assert_(dinput2);

	/*loop over vertices: */
	Gauss* gauss=this->NewGauss();
	for(int month=0;month<12;month++) {
		/*Recover monthly temperatures and precipitation and compute the yearly mean temperatures*/

		for(int iv=0;iv<NUM_VERTICES;iv++) {
			gauss->GaussVertex(iv);
			dinput->GetInputValue(&monthlytemperatures[iv*12+month],gauss,month);
			// yearlytemperatures[iv]=yearlytemperatures[iv]+monthlytemperatures[iv*12+month]*mavg; // Has to be in Kelvin
			monthlytemperatures[iv*12+month]=monthlytemperatures[iv*12+month]-273.15; // conversion from Kelvin to celcius for PDD module
			dinput2->GetInputValue(&monthlyprec[iv*12+month],gauss,month);
			monthlyprec[iv*12+month]=monthlyprec[iv*12+month]*yts;
		}
	}

	/*Recover Pfac, Tdiff and sealev at time t:
	 *     This parameters are used to interpolate the temperature
	 *         and precipitaton between PD and LGM when ismungsm==1 */
	if (ismungsm==1){
		this->parameters->FindParam(&TdiffTime,SmbTdiffEnum,time);
		this->parameters->FindParam(&sealevTime,SmbSealevEnum,time);
	}
	else {
		TdiffTime=0;
		sealevTime=0;
	}

	/*Recover pdd factors at time t.
	 *     This parameter is set, if the user wants to define the
	 *     pdd factors regionally, if issetpddfac==1 in the d18opdd method */
	Input* input  = NULL;
	Input* input2 = NULL;
	if(issetpddfac==1){
		input  = this->GetInput(SmbPddfacSnowEnum); _assert_(input);
		input2 = this->GetInput(SmbPddfacIceEnum); _assert_(input2);
	}

	/*Recover info at the vertices: */
	GetInputListOnVertices(&h[0],ThicknessEnum);
	GetInputListOnVertices(&s[0],SurfaceEnum);
	GetInputListOnVertices(&s0p[0],SmbS0pEnum);
	GetInputListOnVertices(&s0t[0],SmbS0tEnum);

	//Get accumulated input
	Input *accsumM_input = NULL;
	Input *accsumP_input = NULL;
	Input *accsumSMB_input = NULL;
	if (isinitialized){
		accsumM_input = this->GetInput(SmbAccumulatedMeltEnum);  _assert_(accsumM_input);
		accsumP_input = this->GetInput(SmbAccumulatedPrecipitationEnum);  _assert_(accsumP_input);
		accsumSMB_input = this->GetInput(SmbAccumulatedMassBalanceEnum);  _assert_(accsumSMB_input);
	}

	/*measure the surface mass balance*/
	for(int iv = 0; iv<NUM_VERTICES; iv++){
		gauss->GaussVertex(iv);
		pddsnowfac=0.;
		pddicefac=0.;
		if(issetpddfac==1){
			input->GetInputValue(&pddsnowfac,gauss);
			input2->GetInputValue(&pddicefac,gauss);
		}
		agd[iv]=PddSurfaceMassBalance(&monthlytemperatures[iv*12], &monthlyprec[iv*12],
					pdds, pds, &melt[iv], &accu[iv], signorm, yts, h[iv], s[iv],
					desfac, s0t[iv], s0p[iv],rlaps,rlapslgm,TdiffTime,sealevTime,
					pddsnowfac,pddicefac,rho_water,rho_ice);
		/*Get yearlytemperatures */
		for(int month=0;month<12;month++) {
			yearlytemperatures[iv]=yearlytemperatures[iv]+(monthlytemperatures[iv*12+month]+273.15)*mavg; // Has to be in Kelvin
		}

	}

	/*Save accumulated inputs {{{*/
	Input *avgM_input = NULL;
	Input *avgP_input = NULL;
	Input *avgSMB_input = NULL;

	if (isinitialized){
		accsumM_input->GetInputAverage(&accsumM);
		accsumP_input->GetInputAverage(&accsumP);
		accsumSMB_input->GetInputAverage(&accsumSMB);
	}

	/*}}}*/

	/*Update inputs*/
	switch(this->ObjectEnum()){
		case TriaEnum:
			this->AddInput(TemperaturePDDEnum,&yearlytemperatures[0],P1Enum);
			this->AddInput(SmbMassBalanceEnum,&agd[0],P1Enum);
			this->AddInput(SmbAccumulationEnum,&accu[0],P1Enum);
			this->AddInput(SmbMeltEnum,&melt[0],P1Enum);
			break;
		case PentaEnum:
			bool isthermal;
         this->parameters->FindParam(&isthermal,TransientIsthermalEnum);
         if(isthermal){
				if(IsOnSurface()){
					/*Here, we want to change the BC of the thermal model, keep
					 * the temperatures as they are for the base of the penta and
					 * yse yearlytemperatures for the top*/
					PentaInput* temp_input = xDynamicCast<PentaInput*>(this->GetInput(TemperatureEnum)); _assert_(temp_input);
					switch(temp_input->GetInputInterpolationType()){
						case P1Enum:
							temp_input->element_values[3] = yearlytemperatures[3];
							temp_input->element_values[4] = yearlytemperatures[4];
							temp_input->element_values[5] = yearlytemperatures[5];
							temp_input->SetInput(P1Enum,NUM_VERTICES,&vertexlids[0],temp_input->element_values);
							break;
						case P1xP2Enum:
						case P1xP3Enum:
						case P1xP4Enum:
						case P1DGEnum:
							temp_input->element_values[3] = yearlytemperatures[3];
							temp_input->element_values[4] = yearlytemperatures[4];
							temp_input->element_values[5] = yearlytemperatures[5];
							temp_input->SetInput(temp_input->GetInputInterpolationType(),this->lid,this->GetNumberOfNodes(temp_input->GetInputInterpolationType()),temp_input->element_values);
							break;
						default:
							_error_("Interpolation "<<EnumToStringx(temp_input->GetInputInterpolationType())<<" not supported yet");
					}

					bool isenthalpy;
					this->parameters->FindParam(&isenthalpy,ThermalIsenthalpyEnum);
					if(isenthalpy){
						/*Convert that to enthalpy for the enthalpy model*/
						PentaInput* enth_input = xDynamicCast<PentaInput*>(this->GetInput(EnthalpyEnum)); _assert_(enth_input);
						switch(enth_input->GetInputInterpolationType()){
							case P1Enum:
								ThermalToEnthalpy(&enth_input->element_values[3],yearlytemperatures[3],0.,0.);
								ThermalToEnthalpy(&enth_input->element_values[4],yearlytemperatures[4],0.,0.);
								ThermalToEnthalpy(&enth_input->element_values[5],yearlytemperatures[5],0.,0.);
								enth_input->SetInput(P1Enum,NUM_VERTICES,&vertexlids[0],enth_input->element_values);
								break;
							case P1xP2Enum:
							case P1xP3Enum:
							case P1xP4Enum:
							case P1DGEnum:
								ThermalToEnthalpy(&enth_input->element_values[3],yearlytemperatures[3],0.,0.);
								ThermalToEnthalpy(&enth_input->element_values[4],yearlytemperatures[4],0.,0.);
								ThermalToEnthalpy(&enth_input->element_values[5],yearlytemperatures[5],0.,0.);
								enth_input->SetInput(enth_input->GetInputInterpolationType(),this->lid,this->GetNumberOfNodes(enth_input->GetInputInterpolationType()),enth_input->element_values);
								break;
							default:
								_error_("Interpolation "<<EnumToStringx(temp_input->GetInputInterpolationType())<<" not supported yet");
						}
					}
				}
			}
			this->AddInput(SmbMassBalanceEnum,&agd[0],P1Enum);
			this->AddInput(TemperaturePDDEnum,&yearlytemperatures[0],P1Enum);
			this->AddInput(SmbAccumulationEnum,&accu[0],P1Enum);
			this->AddInput(SmbMeltEnum,&melt[0],P1Enum);
			this->InputExtrude(TemperaturePDDEnum,-1);
			this->InputExtrude(SmbMassBalanceEnum,-1);
			this->InputExtrude(SmbAccumulationEnum,-1);
			this->InputExtrude(SmbMeltEnum,-1);
			break;
		default: _error_("Not implemented yet");
	}

	avgM_input = this->GetInput(SmbMeltEnum);  _assert_(avgM_input);
	avgP_input = this->GetInput(SmbAccumulationEnum);  _assert_(avgP_input);
	avgSMB_input = this->GetInput(SmbMassBalanceEnum);  _assert_(avgSMB_input);

	avgM_input->GetInputAverage(&avgM);
	avgP_input->GetInputAverage(&avgP);
	avgSMB_input->GetInputAverage(&avgSMB);

	this->SetElementInput(SmbAccumulatedMassBalanceEnum,accsumSMB+avgSMB*dt);
	this->SetElementInput(SmbAccumulatedPrecipitationEnum,accsumP+avgP*dt);
	this->SetElementInput(SmbAccumulatedMeltEnum,accsumM+avgM*dt);
	if (!isinitialized){
		/*Flag the initialization:*/
		this->SetBoolInput(this->inputs,SmbIsInitializedEnum,true);
	}

	/*clean-up*/
	delete gauss;
	xDelete<IssmDouble>(monthlytemperatures);
	xDelete<IssmDouble>(monthlyprec);
	xDelete<IssmDouble>(agd);
	xDelete<IssmDouble>(melt);
	xDelete<IssmDouble>(accu);
	xDelete<IssmDouble>(yearlytemperatures);
	xDelete<IssmDouble>(h);
	xDelete<IssmDouble>(s);
	xDelete<IssmDouble>(s0t);
	xDelete<IssmDouble>(s0p);
	xDelete<IssmDouble>(tmp);
}
/*}}}*/
void       Element::PositiveDegreeDaySicopolis(bool isfirnwarming){/*{{{*/

	/* General FIXMEs: get Tmelting point, pddicefactor, pddsnowfactor, sigma from parameters/user input */

	const int NUM_VERTICES 		= this->GetNumberOfVertices();
	const int NUM_VERTICES_MONTHS_PER_YEAR	= NUM_VERTICES * 12;

	int        	i,vertexlids[MAXVERTICES];;
	IssmDouble* smb=xNew<IssmDouble>(NUM_VERTICES);		// surface mass balance
	IssmDouble* melt=xNew<IssmDouble>(NUM_VERTICES);		// melting comp. of surface mass balance
	IssmDouble* accu=xNew<IssmDouble>(NUM_VERTICES);		// accuumulation comp. of surface mass balance
	IssmDouble* melt_star=xNew<IssmDouble>(NUM_VERTICES);
	IssmDouble* monthlytemperatures=xNew<IssmDouble>(NUM_VERTICES_MONTHS_PER_YEAR);
	IssmDouble* monthlyprec=xNew<IssmDouble>(NUM_VERTICES_MONTHS_PER_YEAR);
	IssmDouble* yearlytemperatures=xNew<IssmDouble>(NUM_VERTICES); memset(yearlytemperatures, 0., NUM_VERTICES*sizeof(IssmDouble));
	IssmDouble* s=xNew<IssmDouble>(NUM_VERTICES);			// actual surface height
	IssmDouble* s0p=xNew<IssmDouble>(NUM_VERTICES);		// reference elevation for precip.
	IssmDouble* s0t=xNew<IssmDouble>(NUM_VERTICES);		// reference elevation for temperature
	IssmDouble* smbcorr=xNew<IssmDouble>(NUM_VERTICES); // surface mass balance correction; will be added after pdd call
	IssmDouble* p_ampl=xNew<IssmDouble>(NUM_VERTICES);	// precip anomaly
	IssmDouble* t_ampl=xNew<IssmDouble>(NUM_VERTICES);	// remperature anomaly
	IssmDouble rho_water,rho_ice,desfac,rlaps;
	IssmDouble pdd_fac_ice,pdd_fac_snow;
	IssmDouble inv_twelve=1./12.;								//factor for monthly average
	IssmDouble time,yts,time_yr;

	/*Get vertex Lids for later*/
	this->GetVerticesLidList(&vertexlids[0]);

	/*Get material parameters :*/
	rho_water=this->FindParam(MaterialsRhoSeawaterEnum);
	rho_ice=this->FindParam(MaterialsRhoIceEnum);

	/*Get parameters for height corrections*/
	desfac=this->FindParam(SmbDesfacEnum);
	rlaps=this->FindParam(SmbRlapsEnum);

	/*Get pdd melt factors*/
	pdd_fac_ice=this->FindParam(PddfacIceEnum);
	pdd_fac_snow=this->FindParam(PddfacSnowEnum);

	/* Get time */
	this->parameters->FindParam(&time,TimeEnum);
	this->parameters->FindParam(&yts,ConstantsYtsEnum);
	time_yr=floor(time/yts)*yts;

	/* Set parameters for finrnwarming */
	IssmDouble MU_0         = 9.7155; //Firn-warming correction, in (d*deg C)/(mm WE)
	IssmDouble mu           = MU_0*(1000.0*86400.0)*(rho_ice/rho_water);   // (d*deg C)/(mm WE) --> (s*deg C)/(m IE)

	/*Get inputs*/
	DatasetInput* dinput =this->GetDatasetInput(SmbMonthlytemperaturesEnum); _assert_(dinput);
	DatasetInput* dinput2=this->GetDatasetInput(SmbPrecipitationEnum);       _assert_(dinput2);

	/*loop over vertices: */
	Gauss* gauss=this->NewGauss();
	for(int month=0;month<12;month++){

		for(int iv=0;iv<NUM_VERTICES;iv++){
			gauss->GaussVertex(iv);
			dinput->GetInputValue(&monthlytemperatures[iv*12+month],gauss,month);
			monthlytemperatures[iv*12+month]=monthlytemperatures[iv*12+month]-273.15; // conversion from Kelvin to celcius for PDD module
			dinput2->GetInputValue(&monthlyprec[iv*12+month],gauss,month);
			monthlyprec[iv*12+month]=monthlyprec[iv*12+month]*yts;
		}
	}

	/*Recover info at the vertices: */
	GetInputListOnVertices(&s[0],SurfaceEnum);
	GetInputListOnVertices(&s0p[0],SmbS0pEnum);
	GetInputListOnVertices(&s0t[0],SmbS0tEnum);
	GetInputListOnVertices(&smbcorr[0],SmbSmbCorrEnum);
	GetInputListOnVertices(&t_ampl[0],SmbTemperaturesAnomalyEnum);
	GetInputListOnVertices(&p_ampl[0],SmbPrecipitationsAnomalyEnum);

	/*measure the surface mass balance*/
	for (int iv = 0; iv<NUM_VERTICES; iv++){
		smb[iv]=PddSurfaceMassBalanceSicopolis(&monthlytemperatures[iv*12], &monthlyprec[iv*12],
					&melt[iv], &accu[iv], &melt_star[iv], &t_ampl[iv], &p_ampl[iv], yts, s[iv],
					desfac, s0t[iv], s0p[iv],rlaps,rho_water,rho_ice,pdd_fac_ice,pdd_fac_snow);

		/* make correction */
		smb[iv] = smb[iv]+smbcorr[iv];
		/*Get yearlytemperatures */
		for(int month=0;month<12;month++) yearlytemperatures[iv]=yearlytemperatures[iv]+((monthlytemperatures[iv*12+month]+273.15)+t_ampl[iv])*inv_twelve; // Has to be in Kelvin

		if(isfirnwarming){
			if(melt_star[iv]>=melt[iv]){
				yearlytemperatures[iv]= yearlytemperatures[iv]+mu*(melt_star[iv]-melt[iv]);
			}
			else{
				yearlytemperatures[iv]= yearlytemperatures[iv];
			}
		}
		if (yearlytemperatures[iv]>273.15) yearlytemperatures[iv]=273.15;
	}

	switch(this->ObjectEnum()){
		case TriaEnum:
			//this->AddInput(TemperatureEnum,&yearlytemperatures[0],P1Enum);
			this->AddInput(TemperaturePDDEnum,&yearlytemperatures[0],P1Enum);
			this->AddInput(SmbMassBalanceEnum,&smb[0],P1Enum);
			this->AddInput(SmbAccumulationEnum,&accu[0],P1Enum);
			this->AddInput(SmbMeltEnum,&melt[0],P1Enum);
			break;
		case PentaEnum:
			bool isthermal;
			this->parameters->FindParam(&isthermal,TransientIsthermalEnum);
			if(isthermal){
				bool isenthalpy;
				this->parameters->FindParam(&isenthalpy,ThermalIsenthalpyEnum);
				if(IsOnSurface()){
					/*Here, we want to change the BC of the thermal model, keep
					 * the temperatures as they are for the base of the penta and
					 * use yearlytemperatures for the top*/

					/*FIXME: look at other function Element::PositiveDegreeDay and propagate change! Just assert for now*/
					PentaInput* temp_input = xDynamicCast<PentaInput*>(this->GetInput(TemperatureEnum)); _assert_(temp_input);
					switch(temp_input->GetInputInterpolationType()){
						case P1Enum:
							temp_input->element_values[3] = yearlytemperatures[3];
							temp_input->element_values[4] = yearlytemperatures[4];
							temp_input->element_values[5] = yearlytemperatures[5];
							temp_input->SetInput(P1Enum,NUM_VERTICES,&vertexlids[0],temp_input->element_values);
							break;
						case P1DGEnum:
						case P1xP2Enum:
						case P1xP3Enum:
						case P1xP4Enum:
							temp_input->element_values[3] = yearlytemperatures[3];
							temp_input->element_values[4] = yearlytemperatures[4];
							temp_input->element_values[5] = yearlytemperatures[5];
							temp_input->SetInput(temp_input->GetInputInterpolationType(),this->lid,this->GetNumberOfNodes(temp_input->GetInputInterpolationType()),temp_input->element_values);
							break;
						default:
							_error_("Interpolation "<<EnumToStringx(temp_input->GetInputInterpolationType())<<" not supported yet");
					}

					if(isenthalpy){
						/*Convert that to enthalpy for the enthalpy model*/
						PentaInput* enth_input = xDynamicCast<PentaInput*>(this->GetInput(EnthalpyEnum)); _assert_(enth_input);
						switch(enth_input->GetInputInterpolationType()){
							case P1Enum:
								ThermalToEnthalpy(&enth_input->element_values[3],yearlytemperatures[3],0.,0.);
								ThermalToEnthalpy(&enth_input->element_values[4],yearlytemperatures[4],0.,0.);
								ThermalToEnthalpy(&enth_input->element_values[5],yearlytemperatures[5],0.,0.);
								enth_input->SetInput(P1Enum,NUM_VERTICES,&vertexlids[0],enth_input->element_values);
								break;
							case P1DGEnum:
							case P1xP2Enum:
							case P1xP3Enum:
							case P1xP4Enum:
								ThermalToEnthalpy(&enth_input->element_values[3],yearlytemperatures[3],0.,0.);
								ThermalToEnthalpy(&enth_input->element_values[4],yearlytemperatures[4],0.,0.);
								ThermalToEnthalpy(&enth_input->element_values[5],yearlytemperatures[5],0.,0.);
								enth_input->SetInput(enth_input->GetInputInterpolationType(),this->lid,this->GetNumberOfNodes(enth_input->GetInputInterpolationType()),enth_input->element_values);
								break;
							default:
								_error_("Interpolation "<<EnumToStringx(temp_input->GetInputInterpolationType())<<" not supported yet");
						}
					}
				}
			}
			this->AddInput(SmbMassBalanceEnum,&smb[0],P1Enum);
			this->AddInput(TemperaturePDDEnum,&yearlytemperatures[0],P1Enum);
			this->AddInput(SmbAccumulationEnum,&accu[0],P1Enum);
			this->AddInput(SmbMeltEnum,&melt[0],P1Enum);
			this->InputExtrude(TemperaturePDDEnum,-1);
			this->InputExtrude(SmbMassBalanceEnum,-1);
			this->InputExtrude(SmbAccumulationEnum,-1);
			this->InputExtrude(SmbMeltEnum,-1);
			break;
		default: _error_("Not implemented yet");
	}

	/*clean-up*/
	delete gauss;
	xDelete<IssmDouble>(monthlytemperatures);
	xDelete<IssmDouble>(monthlyprec);
	xDelete<IssmDouble>(smb);
	xDelete<IssmDouble>(melt);
	xDelete<IssmDouble>(accu);
	xDelete<IssmDouble>(yearlytemperatures);
	xDelete<IssmDouble>(s);
	xDelete<IssmDouble>(s0t);
	xDelete<IssmDouble>(s0p);
	xDelete<IssmDouble>(t_ampl);
	xDelete<IssmDouble>(p_ampl);
	xDelete<IssmDouble>(smbcorr);
	xDelete<IssmDouble>(melt_star);
}
/*}}}*/
void       Element::PositiveDegreeDayGCM(){/*{{{*/

	const int NUM_VERTICES 	= this->GetNumberOfVertices();
	/*Allocate all arrays*/
	IssmDouble* smb         = xNew<IssmDouble>(NUM_VERTICES);
	IssmDouble* accumulation= xNew<IssmDouble>(NUM_VERTICES);
	IssmDouble* ablation    = xNew<IssmDouble>(NUM_VERTICES);
	IssmDouble* downscaleT  = xNew<IssmDouble>(NUM_VERTICES);
	IssmDouble* refreezing	= xNew<IssmDouble>(NUM_VERTICES);
	/*variables*/
	IssmDouble rlaps, ref_surf, surface;
	IssmDouble temp_all_solid,temp_all_liquid,solid_fraction; 
	IssmDouble temperature, precepitation, annual_temperature;
	IssmDouble ddf_snow, ddf_ice, ddf_firn, ddf_debris, efactor;
	/*Load parameters*/
	this->parameters->FindParam(&temp_all_solid, SmbAllSolidTempEnum);
	this->parameters->FindParam(&temp_all_liquid, SmbAllLiquidTempEnum);
	this->parameters->FindParam(&ddf_snow, SmbDdfSnowEnum);
	this->parameters->FindParam(&ddf_ice, SmbDdfIceEnum);
	ddf_firn = (ddf_snow+ddf_ice)*0.5;
	/*Load inputs*/
	Input *T_input			= this->GetInput(SmbTemperatureEnum);			_assert_(T_input);
	Input *P_input			= this->GetInput(SmbPrecipitationEnum);		_assert_(P_input);
	Input *meanT_input	= this->GetInput(SmbMeanTemperatureEnum);		_assert_(meanT_input);
	Input *efactor_input	= this->GetInput(SmbEnhanceFactorEnum);		_assert_(efactor_input);
	Input *ref_s_input	= this->GetInput(SmbGCMRefSurfaceEnum);		_assert_(ref_s_input);
	Input *surface_input = this->GetInput(SurfaceEnum);					_assert_(surface_input);
	Input *lapserate_input	= this->GetInput(SmbGCMLapseratesEnum);		_assert_(lapserate_input);

	Gauss* gauss=this->NewGauss();

	for(int iv=0;iv<NUM_VERTICES;iv++) {
		gauss->GaussVertex(iv);
		/*Step 1: get GCM data from Input*/
		T_input->GetInputValue(&temperature,gauss);
		P_input->GetInputValue(&precepitation,gauss);
		meanT_input->GetInputValue(&annual_temperature,gauss);
		efactor_input->GetInputValue(&efactor,gauss);
		lapserate_input->GetInputValue(&rlaps,gauss);
		ref_s_input->GetInputValue(&ref_surf,gauss);
		surface_input->GetInputValue(&surface,gauss);

		/*Step 2: downscaling temp with elevation and lapse rate*/
		temperature = temperature + rlaps*(ref_surf - surface);
		downscaleT[iv] = temperature;

		/*Step 3: Accumulation*/
		if (temperature <= temp_all_solid){
			solid_fraction = 1.;
		}
		else if (temperature >= temp_all_liquid){
			solid_fraction = 0.;
		}
		else {
			solid_fraction = (temp_all_liquid - temperature) /(temp_all_liquid - temp_all_solid);
		}
		accumulation[iv] = solid_fraction * precepitation; 

		/*Step 3: Ablation*/
		ddf_debris = ddf_ice*efactor;
		ablation[iv] = (ddf_snow+ddf_firn+ddf_ice+ddf_debris)*(max(0., temperature-273.15));

		/*Step 4: Refreezing*/
		refreezing[iv] = min(ablation[iv], max(0., -0.69*(annual_temperature-273.15)+0.0096)*12/(365*3600*24));

		/*Step 5: TODO: add firn=previous years accumulation-ablation*/
		smb[iv] = accumulation[iv]-ablation[iv]+refreezing[iv];
	}
	/*Add input to element and Free memory*/
	this->AddInput(SmbAccumulationEnum,accumulation,P1Enum);
	this->AddInput(SmbAblationEnum,ablation,P1Enum);
	this->AddInput(SmbDownscaleTemperatureEnum,downscaleT,P1Enum);
	this->AddInput(SmbMassBalanceEnum,smb,P1Enum);
	xDelete<IssmDouble>(accumulation);
	xDelete<IssmDouble>(ablation);
	xDelete<IssmDouble>(downscaleT);
	xDelete<IssmDouble>(refreezing);
	xDelete<IssmDouble>(smb);
	delete gauss;
}
/*}}}*/
void       Element::ProjectGridDataToMesh(IssmDouble* griddata,IssmDouble* x_grid,IssmDouble* y_grid,int Nx,int Ny,int input_enum){/*{{{*/

	const int NUM_VERTICES 	= this->GetNumberOfVertices();
	/*Allocate all arrays*/
	IssmDouble* temp        = xNew<IssmDouble>(NUM_VERTICES);
	IssmDouble* xyz_list = NULL;
	IssmDouble x, y;
	int m,n;

	this->GetVerticesCoordinates(&xyz_list);
	for(int iv=0;iv<NUM_VERTICES;iv++) {
		/*Step 1: loop over all vertices, interpolate GCM temperature and precepitation to the mesh */
		x = xyz_list[iv*3+0];
		y = xyz_list[iv*3+1];
		/*Find indices m and n into y_grid and x_grid, for which  y_grid(m)<=y<=y_grid(m+1) and x_grid(n)<=x<=x_grid(n+1)*/
		findindices<IssmDouble>(&n,&m,x_grid,Nx,y_grid,Ny,x,y);
		temp[iv] = bilinearinterp(x_grid,y_grid,griddata,x,y,m,n,Nx);
	}
	/*Add input to element and Free memory*/
	this->AddInput(input_enum,temp,P1Enum);
	xDelete<IssmDouble>(temp);
	xDelete<IssmDouble>(xyz_list);
}
/*}}}*/
void       Element::SmbDebrisEvatt(){/*{{{*/

	const int NUM_VERTICES          = this->GetNumberOfVertices();
	const int NUM_VERTICES_DAYS_PER_YEAR  = NUM_VERTICES * 365; // 365 FIXME

	int             i,vertexlids[MAXVERTICES];;
	IssmDouble* smb=xNew<IssmDouble>(NUM_VERTICES);
	IssmDouble* melt=xNew<IssmDouble>(NUM_VERTICES);
	IssmDouble* summermelt=xNew<IssmDouble>(NUM_VERTICES); 
	IssmDouble* albedo=xNew<IssmDouble>(NUM_VERTICES);
	IssmDouble* summeralbedo=xNew<IssmDouble>(NUM_VERTICES); 
	IssmDouble* accu=xNew<IssmDouble>(NUM_VERTICES);

	// climate inputs
	IssmDouble* temperature=xNew<IssmDouble>(NUM_VERTICES_DAYS_PER_YEAR);
	IssmDouble* precip=xNew<IssmDouble>(NUM_VERTICES_DAYS_PER_YEAR);
	IssmDouble* lw=xNew<IssmDouble>(NUM_VERTICES_DAYS_PER_YEAR);
	IssmDouble* sw=xNew<IssmDouble>(NUM_VERTICES_DAYS_PER_YEAR);
	IssmDouble* wind=xNew<IssmDouble>(NUM_VERTICES_DAYS_PER_YEAR);
	IssmDouble* humidity=xNew<IssmDouble>(NUM_VERTICES_DAYS_PER_YEAR);
	IssmDouble* yearlytemperatures=xNew<IssmDouble>(NUM_VERTICES); memset(yearlytemperatures, 0., NUM_VERTICES*sizeof(IssmDouble));
	IssmDouble* p_ampl=xNew<IssmDouble>(NUM_VERTICES);
	IssmDouble* t_ampl=xNew<IssmDouble>(NUM_VERTICES);
	IssmDouble* lw_ampl=xNew<IssmDouble>(NUM_VERTICES);
	IssmDouble* sw_ampl=xNew<IssmDouble>(NUM_VERTICES);
	IssmDouble* wind_ampl=xNew<IssmDouble>(NUM_VERTICES);
	IssmDouble* humidity_ampl=xNew<IssmDouble>(NUM_VERTICES);

	IssmDouble* surface=xNew<IssmDouble>(NUM_VERTICES);
	IssmDouble* s0t=xNew<IssmDouble>(NUM_VERTICES);
	IssmDouble* snowheight=xNew<IssmDouble>(NUM_VERTICES);
	IssmDouble* debriscover=xNew<IssmDouble>(NUM_VERTICES);
	IssmDouble rho_water,rho_ice,Tf,debris,debris_here;
	IssmDouble qlaps,rlaps,dsgrad,dlgrad,windspeedgrad,humiditygrad,Tm;
	IssmDouble inv_twelve=1./365.;
	IssmDouble time,yts,time_yr,lambda;
	IssmDouble DailyMelt,CleanIceDailyMelt, CumDailyMelt=0,CleanIceMelt,CumDailySummerMelt=0;
	IssmDouble MeanAlbedo=0, MeanSummerAlbedo=0;
	bool isdebris,isAnderson,iscryokarst;
	this->parameters->FindParam(&isdebris,TransientIsdebrisEnum);
	this->parameters->FindParam(&isAnderson,SmbDebrisIsAndersonEnum);
	this->parameters->FindParam(&iscryokarst,SmbDebrisIsCryokarstEnum);
	IssmDouble PhiD=0.,p;
	IssmDouble icealbedo=this->FindParam(SmbIcealbedoEnum);
	IssmDouble snowalbedo=this->FindParam(SmbSnowalbedoEnum);
	IssmDouble debrisalbedo=this->FindParam(SmbDebrisalbedoEnum);
	IssmDouble Lm=this->FindParam(MaterialsLatentheatEnum); 
	IssmDouble D0=this->FindParam(SmbDebrisAndersonD0Enum);
	int step;
	this->FindParam(&step,StepEnum);

	// cryokarst
	int dim=1,domaintype;
	this->parameters->FindParam(&domaintype,DomainTypeEnum);
	if(domaintype!=Domain2DverticalEnum){
		dim=2;
	}
	IssmDouble taud_plus=110e3, taud_minus=60e3;
	IssmDouble taud, slope, gravity, taudx, taudy;
	this->parameters->FindParam(&gravity,ConstantsGEnum);
	IssmDouble* slopex         = xNew<IssmDouble>(NUM_VERTICES);
	IssmDouble* slopey         = xNew<IssmDouble>(NUM_VERTICES);
	IssmDouble* icethickness   = xNew<IssmDouble>(NUM_VERTICES);

	/*Get material parameters :*/
	rho_water=this->FindParam(MaterialsRhoSeawaterEnum);
	rho_ice=this->FindParam(MaterialsRhoIceEnum);
	IssmDouble sconv=(rho_water/rho_ice); 
	Tf=this->FindParam(MaterialsMeltingpointEnum);

	/*Get parameters for height corrections*/
	qlaps=this->FindParam(SmbDesfacEnum); // comment MR; on alpine galciers we dont have the desertification effect
	rlaps=this->FindParam(SmbRlapsEnum);
	dsgrad=this->FindParam(SmbSWgradEnum);
	dlgrad=this->FindParam(SmbLWgradEnum);
	windspeedgrad=this->FindParam(SmbWindspeedgradEnum);
	humiditygrad=this->FindParam(SmbHumiditygradEnum);

	/* Get time */
	this->parameters->FindParam(&time,TimeEnum);
	this->parameters->FindParam(&yts,ConstantsYtsEnum);
	time_yr=floor(time/yts)*yts;

	/*Get inputs*/
	DatasetInput* tempday     =this->GetDatasetInput(SmbMonthlytemperaturesEnum); _assert_(tempday);
	DatasetInput* precipday   =this->GetDatasetInput(SmbPrecipitationEnum);       _assert_(precipday);
	DatasetInput* lwday       =this->GetDatasetInput(SmbMonthlydlradiationEnum); _assert_(lwday);
	DatasetInput* swday       =this->GetDatasetInput(SmbMonthlydsradiationEnum);       _assert_(swday);
	DatasetInput* windday     =this->GetDatasetInput(SmbMonthlywindspeedEnum); _assert_(windday);
	DatasetInput* humidityday =this->GetDatasetInput(SmbMonthlyairhumidityEnum); _assert_(humidityday);

	/*loop over vertices: */
	Gauss* gauss=this->NewGauss();
	for(int month=0;month<365;month++){
		for(int iv=0;iv<NUM_VERTICES;iv++){
			gauss->GaussVertex(iv);
			tempday->GetInputValue(&temperature[iv*365+month],gauss,month);
			temperature[iv*365+month]=temperature[iv*365+month]-Tf; // conversion from Kelvin to celcius for PDD module
			precipday->GetInputValue(&precip[iv*365+month],gauss,month);
			precip[iv*365+month]=precip[iv*365+month]*yts; // from m/s to m/a
			lwday->GetInputValue(&lw[iv*365+month],gauss,month);
			swday->GetInputValue(&sw[iv*365+month],gauss,month);
			windday->GetInputValue(&wind[iv*365+month],gauss,month);
			humidityday->GetInputValue(&humidity[iv*365+month],gauss,month);
		}
	}

	/*Recover info at the vertices: */
	GetInputListOnVertices(&surface[0],SurfaceEnum);
	GetInputListOnVertices(&s0t[0],SmbS0tEnum);
	GetInputListOnVertices(&snowheight[0],SmbSnowheightEnum);
	GetInputListOnVertices(&debriscover[0],DebrisThicknessEnum);
	GetInputListOnVertices(&t_ampl[0],SmbTemperaturesAnomalyEnum);
	GetInputListOnVertices(&p_ampl[0],SmbPrecipitationsAnomalyEnum);
	GetInputListOnVertices(&lw_ampl[0],SmbDsradiationAnomalyEnum);
	GetInputListOnVertices(&sw_ampl[0],SmbDlradiationAnomalyEnum);
	GetInputListOnVertices(&wind_ampl[0],SmbWindspeedAnomalyEnum);
	GetInputListOnVertices(&humidity_ampl[0],SmbAirhumidityAnomalyEnum);
	if(iscryokarst){
		GetInputListOnVertices(&slopex[0],SurfaceSlopeXEnum);
		GetInputListOnVertices(&icethickness[0],ThicknessEnum);
		if(dim==2){
			GetInputListOnVertices(&slopey[0],SurfaceSlopeYEnum);
		}
		taudx=rho_ice*gravity*icethickness[i]*slopex[i];
		if(dim==2) taudy=rho_ice*gravity*icethickness[i]*slopey[i];
		taud=sqrt(taudx*taudx+taudy*taudy);
	}
	IssmDouble Alphaeff,Alphaeff_cleanice;

	/*measure the surface mass balance*/
	for (int iv = 0; iv<NUM_VERTICES; iv++){

		IssmDouble st=(surface[iv]-s0t[iv])/1000.;

		int ismb_end=1;
		if(isdebris & !isAnderson) ismb_end=2;
		for (int ismb=0;ismb<ismb_end;ismb++){
			if(ismb==0){
				// calc a reference smb to identify accum and melt region; debris only develops in ablation area
				debris=0.;
				PhiD=0.;
				if(isAnderson) debris_here=debriscover[iv]; // store debris for later
			}else{
				// debris only develops in ablation area
				/*if((accu[iv]/yts-CleanIceMelt)<(-1e-2)/yts){
				  debris=debriscover[iv];
				  }else{
				  debris=0.;
				  }*/
				debris=0.;
				if(debris<=0.) debris=0.;
				if(isdebris) PhiD=FindParam(DebrisPackingFractionEnum);
				CumDailyMelt=0;
				CumDailySummerMelt=0;
				debris_here=debriscover[iv];
			}

			/* Now run the debris part */

			// Climate inputs
			IssmDouble Tm;          // C air temperature
			IssmDouble In;          // Wm^-2 incoming long wave
			IssmDouble Q;           // Wm^-2 incoming short wave
			IssmDouble Um;          // ms^-1 measured wind speed
			IssmDouble Humidity;    // relative humidity
			IssmDouble P;           // precip

			// other parameters
			IssmDouble Qh=0.006;   // kg m^-3      saturated humidity level // not used
			IssmDouble Qm=0.8*Qh;  // kg m^-3      measured humiditiy level // not used
			IssmDouble Rhoaa=1.22; // kgm^-3       air densitiy
			IssmDouble K=0.585;    // Wm^-1K^-1    thermal conductivity          0.585
			IssmDouble Xr=0.01;    // ms^-1        surface roughness             0.01
			IssmDouble Ustar=0.16; // ms^-1        friction velocity             0.16
			IssmDouble Ca=1000;    // jkg^-1K^-1   specific heat capacity of air
			IssmDouble Lv=2.50E+06;// jkg^-1K^-1   latent heat of evaporation
			IssmDouble Eps=0.95;   //              thermal emissivity
			IssmDouble Sigma=5.67E-08;// Wm^-2K^-4    Stefan Boltzmann constant
			IssmDouble Gamma=180.;    // m^-1         wind speed attenuation        234

			// Calculate effective albedo
			IssmDouble Alphaeff,Alphaeff_cleanice;
			IssmDouble mean_ela,delta=2000;

			// compute cleanice albedo based on previous SMB distribution
			//if(step==1){
			mean_ela=3000; //FIXME
								//}else{
								//        mean_ela=FindParam(SmbMeanElaEnum);
								//}
			Alphaeff_cleanice=icealbedo+(snowalbedo-icealbedo)*(1+tanh(PI*(surface[iv]-mean_ela)/delta))/2;
			Alphaeff=Alphaeff_cleanice; // will be updated below

			accu[iv]=0.;
			for (int iday=0;iday<365;iday++) {

				Tm=temperature[iv*365+iday]-st*rlaps;//+t_ampl[iv];//+(rand()%10-5)/5;
				In=lw[iv*365+iday]-st*dlgrad+lw_ampl[iv];
				Q=sw[iv*365+iday]+st*dsgrad+sw_ampl[iv];
				Humidity=humidity[iv*365+iday]-st*humiditygrad+humidity_ampl[iv];
				Um=wind[iv*365+iday]-st*windspeedgrad+wind_ampl[iv];
				P=(qlaps*st*precip[iv*365+iday]+precip[iv*365+iday]+p_ampl[iv])*sconv/365.; // convert precip from w.e. -> i.e

				/*Partition of precip in solid and liquid parts */
				IssmDouble temp_plus=1; 
				IssmDouble temp_minus=-1.;
				IssmDouble frac_solid;
				if(Tm>=temp_plus){
					frac_solid=0;
				}else if(Tm<=temp_minus){
					frac_solid=1;
				}else{
					frac_solid=1*(1-cos(PI*(temp_plus-Tm)/(temp_plus-temp_minus)))/2;
				}

				/*Get yearly temperatures and accumulation */
				yearlytemperatures[iv]=yearlytemperatures[iv]+((temperature[iv*365+iday]-rlaps*st+Tf+t_ampl[iv]))/365; // Has to be in Kelvin
				accu[iv]=accu[iv]+P*frac_solid;
				if(yearlytemperatures[iv]>Tf) yearlytemperatures[iv]=Tf;

				CleanIceDailyMelt=((In-(Eps*Sigma*(Tf*Tf*Tf*Tf))+
								Q*(1.-Alphaeff)+
								(Rhoaa*Ca*Ustar*Ustar)/(Um-Ustar*(2.-(exp(Gamma*Xr))))*Tm)/((1-PhiD)*rho_ice*Lm)/(1.+
								((Rhoaa*Ca*Ustar*Ustar)/(Um-Ustar*(2.-(exp(Gamma*Xr))))+4.*Eps*Sigma*(Tf*Tf*Tf))/
								K*debris)-(Lv*Ustar*Ustar*((Humidity))*(exp(-Gamma*Xr)))/((1.-PhiD)*
										rho_ice*Lm*Ustar)/(((Um
														-2.*Ustar)*exp(-Gamma*Xr))/Ustar+exp(Gamma*debris)));
				if(CleanIceDailyMelt<0) CleanIceDailyMelt=0.;
				DailyMelt=CleanIceDailyMelt;

				if(ismb==1){

					//snowheight[iv]=snowheight[iv]+(P-CleanIceDailyMelt*yts/365);
					IssmDouble sn_prev;
					sn_prev=snowheight[iv];
					snowheight[iv]=sn_prev+(-CleanIceDailyMelt*yts/365);//P

					if(snowheight[iv]<=0) snowheight[iv]=0.;
					if(snowheight[iv]<=0.0001){
						p=debris_here*PhiD/(2*0.2*0.01); //Eq. 51 from Evatt et al 2015 without source term g*t
						if(p>1.) p=1.;
						if(p>=0.999){
							Alphaeff=debrisalbedo;
						} else {
							Alphaeff=Alphaeff_cleanice+p*(debrisalbedo-Alphaeff_cleanice);
						}
						debris=debris_here;
						DailyMelt=((In-(Eps*Sigma*(Tf*Tf*Tf*Tf))+
										Q*(1.-Alphaeff)+
										(Rhoaa*Ca*Ustar*Ustar)/(Um-Ustar*(2.-(exp(Gamma*Xr))))*Tm)/((1-PhiD)*rho_ice*Lm)/(1.+
										((Rhoaa*Ca*Ustar*Ustar)/(Um-Ustar*(2.-(exp(Gamma*Xr))))+4.*Eps*Sigma*(Tf*Tf*Tf))/
										K*debris)-(Lv*Ustar*Ustar*((Humidity))*(exp(-Gamma*Xr)))/((1.-PhiD)*
												rho_ice*Lm*Ustar)/(((Um-2.*Ustar)*exp(-Gamma*Xr))/Ustar+exp(Gamma*debris)));
						if(DailyMelt<0) DailyMelt=0.;
						MeanSummerAlbedo=MeanSummerAlbedo+Alphaeff;
						CumDailySummerMelt=CumDailySummerMelt+DailyMelt/365;
					}
				}
				CumDailyMelt=CumDailyMelt+DailyMelt/365;
			}
			MeanAlbedo=MeanAlbedo+Alphaeff;
			if(ismb==0) CleanIceMelt=CumDailyMelt;
		}

		if(iscryokarst){
			if(taud>=taud_plus){
				lambda=0;
			}else if(taud>=taud_minus & taud<taud_plus){
				lambda=0.1*(1-cos(PI*(taud_plus-taud)/(taud_plus-taud_minus)))/2;
			}else if(taud<taud_minus){
				lambda=0.1;
			}
		}

		// update values
		melt[iv]=CumDailyMelt; // is already in m/s
		accu[iv]=accu[iv]/yts;
		if(isAnderson){
			smb[iv]=(accu[iv]-melt[iv])*D0/(D0+debris_here);
			if(iscryokarst){ 
				smb[iv]=lambda*(accu[iv]-melt[iv])+(1-lambda)*(accu[iv]-melt[iv])*D0/(D0+debris_here);
			}else{
				smb[iv]=(accu[iv]-melt[iv])*D0/(D0+debris_here);
			}
		}else{
			if(iscryokarst){ 
				smb[iv]=lambda*(accu[iv]-CleanIceMelt)+(1-lambda)*(accu[iv]-melt[iv]);
			}else{
				smb[iv]=(accu[iv]-melt[iv]);
			}
		}
		albedo[iv]=MeanAlbedo;
		summeralbedo[iv]=MeanSummerAlbedo;
		summermelt[iv]=CumDailySummerMelt;
	}

	this->AddInput(SmbMassBalanceEnum,smb,P1Enum);
	this->AddInput(SmbAccumulationEnum,accu,P1Enum);
	this->AddInput(SmbMeltEnum,melt,P1Enum);
	this->AddInput(SmbSummerMeltEnum,summermelt,P1Enum);
	this->AddInput(SmbSnowheightEnum,snowheight,P1Enum);
	this->AddInput(SmbAlbedoEnum,albedo,P1Enum);
	this->AddInput(SmbSummerAlbedoEnum,summeralbedo,P1Enum);
	this->AddInput(TemperaturePDDEnum,yearlytemperatures,P1Enum); // TemperaturePDD is wrong here, but don't want to create new Enum ...

	/*clean-up*/
	xDelete<IssmDouble>(temperature);
	xDelete<IssmDouble>(precip);
	xDelete<IssmDouble>(lw);
	xDelete<IssmDouble>(sw);
	xDelete<IssmDouble>(wind);
	xDelete<IssmDouble>(humidity);
	xDelete<IssmDouble>(smb);
	xDelete<IssmDouble>(surface);
	xDelete<IssmDouble>(melt);
	xDelete<IssmDouble>(summermelt);
	xDelete<IssmDouble>(albedo);
	xDelete<IssmDouble>(summeralbedo);
	xDelete<IssmDouble>(accu);
	xDelete<IssmDouble>(yearlytemperatures);
	xDelete<IssmDouble>(s0t);
	xDelete<IssmDouble>(snowheight);
	xDelete<IssmDouble>(debriscover);
	xDelete<IssmDouble>(t_ampl);
	xDelete<IssmDouble>(p_ampl);
	xDelete<IssmDouble>(lw_ampl);
	xDelete<IssmDouble>(sw_ampl);
	xDelete<IssmDouble>(humidity_ampl);
	xDelete<IssmDouble>(wind_ampl);
	xDelete<IssmDouble>(slopex);
	xDelete<IssmDouble>(slopey);
	xDelete<IssmDouble>(icethickness);
}
/*}}}*/
void       Element::ResultInterpolation(int* pinterpolation,int* pnodesperelement,int* parray_size, int output_enum){/*{{{*/

	/*Some intputs need to be computed, even if they are already in inputs, they might not be up to date!*/
	switch(output_enum){
		case ViscousHeatingEnum: this->ViscousHeatingCreateInput(); break;
		case FrictionAlpha2Enum: this->FrictionAlpha2CreateInput(); break;
		case StressMaxPrincipalEnum: this->StressMaxPrincipalCreateInput(); break;
		case StressTensorxxEnum:
		case StressTensorxyEnum:
		case StressTensorxzEnum:
		case StressTensoryyEnum:
		case StressTensoryzEnum:
		case StressTensorzzEnum: this->ComputeStressTensor(); break;
		case StrainRatexxEnum:
		case StrainRatexyEnum:
		case StrainRatexzEnum:
		case StrainRateyyEnum:
		case StrainRateyzEnum:
		case StrainRatezzEnum:
		case StrainRateeffectiveEnum: this->ComputeStrainRate(); break;
		case DeviatoricStressxxEnum:
		case DeviatoricStressxyEnum:
		case DeviatoricStressxzEnum:
		case DeviatoricStressyyEnum:
		case DeviatoricStressyzEnum:
		case DeviatoricStresszzEnum:
		case DeviatoricStress1Enum:
		case DeviatoricStress2Enum:
		case DeviatoricStresseffectiveEnum: this->ComputeDeviatoricStressTensor(); break;
		case EsaStrainratexxEnum:
		case EsaStrainratexyEnum:
		case EsaStrainrateyyEnum:
		case EsaRotationrateEnum: this->ComputeEsaStrainAndVorticity(); break;
		case SigmaNNEnum: this->ComputeSigmaNN(); break;
		case LambdaSEnum: this->ComputeLambdaS(); break;
		case StressIntensityFactorEnum: this->StressIntensityFactor(); break;
		case CalvingratexEnum:
		case CalvingrateyEnum:
		case CalvingCalvingrateEnum:
												  this->StrainRateparallel();
												  this->StrainRateperpendicular();
												  int calvinglaw;
												  this->FindParam(&calvinglaw,CalvingLawEnum);
												  switch(calvinglaw){
													  case DefaultCalvingEnum:
														  //do nothing
														  break;
													  case CalvingLevermannEnum:
														  this->CalvingRateLevermann();
														  break;
													  case CalvingPollardEnum:
														  this->CalvingPollard();
														  break;
													  case CalvingVonmisesEnum:
													  case CalvingDev2Enum:
														  this->CalvingRateVonmises();
														  break;
													  case CalvingVonmisesADEnum:
														  this->CalvingRateVonmisesAD();
														  break;
													  case CalvingCrevasseDepthEnum:
														  this->CalvingCrevasseDepth();
														  break;
													  case CalvingParameterizationEnum:
														  this->CalvingRateParameterization();
														  break;
													  case CalvingCalvingMIPEnum:
														  this->CalvingRateCalvingMIP();
														  break;
													  case CalvingTestEnum:
														  this->CalvingRateTest();
														  break;
													  default:
														  _error_("Calving law "<<EnumToStringx(calvinglaw)<<" not supported yet");
												  }
												  break;
		case CalvingFluxLevelsetEnum: this->CalvingFluxLevelset(); break;
		case CalvingMeltingFluxLevelsetEnum: this->CalvingMeltingFluxLevelset(); break;
		case StrainRateparallelEnum: this->StrainRateparallel(); break;
		case StrainRateperpendicularEnum: this->StrainRateperpendicular(); break;
		case SurfaceCrevasseEnum: this->CalvingCrevasseDepth(); break;
		case SigmaVMEnum: this->ComputeSigmaVM(); break;
		case PartitioningEnum: this->inputs->SetInput(PartitioningEnum,this->lid,IssmComm::GetRank()); break;
	}

	/*If this input is not already in Inputs, maybe it needs to be computed?*/
	switch(this->inputs->GetInputObjectEnum(output_enum)){
		case TriaInputEnum:
		case PentaInputEnum:
		case ControlInputEnum:
		case TransientInputEnum:{
											Input* input2 = this->GetInput(output_enum);
											if(!input2) _error_("input "<<EnumToStringx(output_enum)<<" not found in element");
											*pinterpolation   = input2->GetResultInterpolation();
											*pnodesperelement = input2->GetResultNumberOfNodes();
											*parray_size      = input2->GetResultArraySize();
										}
									 break;
		case BoolInputEnum:
									 *pinterpolation   = P0Enum;
									 *pnodesperelement = 1;
									 *parray_size      = 1;
									 break;
		case IntInputEnum:
									 *pinterpolation   = P0Enum;
									 *pnodesperelement = 1;
									 *parray_size      = 1;
									 break;
		case ArrayInputEnum:{
									  int M;
									  this->inputs->GetArray(output_enum,this->lid,NULL,&M);
									  *pinterpolation   = P0ArrayEnum;
									  *pnodesperelement = 1;
									  *parray_size      = M;
								  }
								break;
		default:
								_error_("Input type \""<<EnumToStringx(this->inputs->GetInputObjectEnum(output_enum))<<"\" not supported yet (While trying to return "<<EnumToStringx(output_enum)<<")");
	}

	/*Assign output pointer*/

	return;
}/*}}}*/
void       Element::ResultToPatch(IssmDouble* values,int nodesperelement,int output_enum){/*{{{*/

	/*Find input*/
	Input* input=this->GetInput(output_enum);
	if(!input) _error_("input "<<EnumToStringx(output_enum)<<" not found in element");

	/*Cast to ElementInput*/
	if(input->ObjectEnum()!=TriaInputEnum && input->ObjectEnum()!=PentaInputEnum){
		_error_("Input "<<EnumToStringx(output_enum)<<" is not an ElementInput");
	}
	ElementInput* element_input = xDynamicCast<ElementInput*>(input);

	/*Get Number of nodes and make sure that it is the same as the one provided*/
	int numnodes = this->GetNumberOfNodes(element_input->GetInputInterpolationType());
	_assert_(numnodes==nodesperelement);

	/*Fill in arrays*/
	for(int i=0;i<numnodes;i++) values[this->sid*numnodes + i] = element_input->element_values[i];

} /*}}}*/
void       Element::ResultToMatrix(IssmDouble* values,int ncols,int output_enum){/*{{{*/

	IssmDouble* array = NULL;
	int         m;
	this->inputs->GetArray(output_enum,this->lid,&array,&m);
	for(int i=0;i<m;i++) values[this->Sid()*ncols + i] = array[i];
	xDelete<IssmDouble>(array);

} /*}}}*/
void       Element::ResultToVector(Vector<IssmDouble>* vector,int output_enum){/*{{{*/

	IssmDouble values[MAXVERTICES];
	int        connectivity[MAXVERTICES];
	int        sidlist[MAXVERTICES];

	switch(this->inputs->GetInputObjectEnum(output_enum)){
		case TriaInputEnum:
		case PentaInputEnum:
		case ControlInputEnum:
		case TransientInputEnum:{

											Input* input2 = this->GetInput(output_enum);
											if(!input2) _error_("input "<<EnumToStringx(output_enum)<<" not found in element");

											switch(input2->GetResultInterpolation()){
												case P0Enum:{
																	IssmDouble  value;
																	bool        bvalue;
																	Gauss* gauss = this->NewGauss();
																	input2->GetInputValue(&value,gauss);
																	delete gauss;
																	vector->SetValue(this->Sid(),value,INS_VAL);
																	break;
																}
												case P1Enum:{
																	const int NUM_VERTICES = this->GetNumberOfVertices();

																	this->GetVerticesSidList(&sidlist[0]);
																	this->GetVerticesConnectivityList(&connectivity[0]);
																	this->GetInputListOnVertices(&values[0],output_enum);
																	for(int i=0;i<NUM_VERTICES;i++) values[i] = values[i]/reCast<IssmDouble>(connectivity[i]);
																	vector->SetValues(NUM_VERTICES,sidlist,values,ADD_VAL);
																	break;
																}
												default:
															 _error_("interpolation "<<EnumToStringx(input2->GetResultInterpolation())<<" not supported yet");
											}
										}
									 break;
		case BoolInputEnum:
									 bool bvalue;
									 this->GetInputValue(&bvalue,output_enum);
									 vector->SetValue(this->Sid(),reCast<IssmDouble>(bvalue),INS_VAL);
									 break;
		case IntInputEnum:
									 int ivalue;
									 this->GetInputValue(&ivalue,output_enum);
									 vector->SetValue(this->Sid(),reCast<IssmDouble>(ivalue),INS_VAL);
									 break;
		default:
									 _error_("Input type \""<<EnumToStringx(this->inputs->GetInputObjectEnum(output_enum))<<"\" not supported yet");
	}

} /*}}}*/
void       Element::RignotMeltParameterization(){/*{{{*/

	const int numvertices = this->GetNumberOfVertices();
	IssmDouble A, B, alpha, beta;
	IssmDouble bed,qsg,qsg_basin,TF,yts;
	int numbasins,basinid;
	IssmDouble* basin_icefront_area=NULL;

	/* Coefficients */
	A    = 3e-4;
	B    = 0.15;
	alpha = 0.39;
	beta = 1.18;

	/*Get inputs*/
	this->GetInputValue(&basinid,FrontalForcingsBasinIdEnum);
	Input* bed_input = this->GetInput(BedEnum);                     _assert_(bed_input);
	Input* qsg_input = this->GetInput(FrontalForcingsSubglacialDischargeEnum);     _assert_(qsg_input);
	Input* TF_input  = this->GetInput(ThermalForcingEnum);          _assert_(TF_input);

	this->FindParam(&yts, ConstantsYtsEnum);
	this->parameters->FindParam(&numbasins,FrontalForcingsNumberofBasinsEnum);
	this->parameters->FindParam(&basin_icefront_area,&numbasins,FrontalForcingsBasinIcefrontAreaEnum);
	IssmDouble meltrates[MAXVERTICES];  

	/* Start looping on the number of vertices: */
	Gauss* gauss=this->NewGauss();
	for(int iv=0;iv<numvertices;iv++){
		gauss->GaussVertex(iv);

		/* Get variables */
		bed_input->GetInputValue(&bed,gauss);
		qsg_input->GetInputValue(&qsg,gauss);
		TF_input->GetInputValue(&TF,gauss);

		if(basin_icefront_area[basinid]==0.) meltrates[iv]=0.;
		else{
			/* change the unit of qsg (m^3/d -> m/d) with ice front area */
			qsg_basin=qsg/basin_icefront_area[basinid];

			/* calculate melt rates */
			meltrates[iv]=((A*max(-bed,0.)*pow(max(qsg_basin,0.),alpha)+B)*pow(max(TF,0.),beta))/86400; //[m/s]
		}

		if(xIsNan<IssmDouble>(meltrates[iv])) _error_("NaN found in vector");
		if(xIsInf<IssmDouble>(meltrates[iv])) _error_("Inf found in vector");
	}

	/*Add input*/
	this->AddInput(CalvingMeltingrateEnum,&meltrates[0],P1Enum);

	/*Cleanup and return*/
	delete gauss;
	xDelete<IssmDouble>(basin_icefront_area);
}
/*}}}*/
void       Element::SetBoolInput(Inputs* inputs,int enum_in,bool value){/*{{{*/

	_assert_(inputs);
	inputs->SetInput(enum_in,this->lid,value);

}
/*}}}*/
void       Element::SetIntInput(Inputs* inputs,int enum_in,int value){/*{{{*/

	_assert_(inputs);
	inputs->SetInput(enum_in,this->lid,value);

}
/*}}}*/
void       Element::SetwiseNodeConnectivity(int* pd_nz,int* po_nz,Node* node,bool* flags,int* flagsindices,int* flagsindices_counter,int set1_enum,int set2_enum, int analysis_type){/*{{{*/

	/*Intermediaries*/
	const int numnodes = this->GetNumberOfNodes();
	/*Output */
	int d_nz = 0;
	int o_nz = 0;

	/*Loop over all nodes*/
	for(int i=0;i<numnodes;i++){
		int nodelid = this->nodes[i]->Lid();
		if(!flags[nodelid]){
			/*flag current node so that no other element processes it*/
			flags[nodelid]=true;

			flagsindices[flagsindices_counter[0]]=nodelid;
			flagsindices_counter[0]++;

			/*if node is clone, we have an off-diagonal non-zero, else it is a diagonal non-zero*/
			switch(set2_enum){
				case FsetEnum:
					if(nodes[i]->FSize()){
						if(this->nodes[i]->IsClone())
						 o_nz += 1;
						else
						 d_nz += 1;
					}
					break;
				case GsetEnum:
					if(nodes[i]->gsize){
						if(this->nodes[i]->IsClone())
						 o_nz += 1;
						else
						 d_nz += 1;
					}
					break;
				case SsetEnum:
					if(nodes[i]->SSize()){
						if(this->nodes[i]->IsClone())
						 o_nz += 1;
						else
						 d_nz += 1;
					}
					break;
				default: _error_("not supported");
			}
		}
	}

	/*Special case: 2d/3d coupling, the node of this element might be connected
	 *to the basal element*/
	int approximation,numlayers;
	if(analysis_type==StressbalanceAnalysisEnum){
		this->GetInputValue(&approximation,ApproximationEnum);
		if(approximation==SSAHOApproximationEnum || approximation==SSAFSApproximationEnum){
			parameters->FindParam(&numlayers,MeshNumberoflayersEnum);
			o_nz += numlayers*3;
			d_nz += numlayers*3;
		}
	}

	/*Assign output pointers: */
	*pd_nz=d_nz;
	*po_nz=o_nz;
}
/*}}}*/
#ifdef _HAVE_SEMIC_
void       Element::SmbSemic(){/*{{{*/

	/*only compute SMB at the surface: */
	if (!IsOnSurface()) return;

	const int NUM_VERTICES 					= this->GetNumberOfVertices();
	const int NUM_VERTICES_DAYS_PER_YEAR 	= NUM_VERTICES * 365;

	IssmDouble* s=xNew<IssmDouble>(NUM_VERTICES);
	IssmDouble* s0gcm=xNew<IssmDouble>(NUM_VERTICES);
	IssmDouble* st=xNew<IssmDouble>(NUM_VERTICES);

	// daily forcing inputs
	IssmDouble* dailyrainfall=xNew<IssmDouble>(NUM_VERTICES_DAYS_PER_YEAR);
	IssmDouble* dailysnowfall=xNew<IssmDouble>(NUM_VERTICES_DAYS_PER_YEAR);
	IssmDouble* dailydlradiation=xNew<IssmDouble>(NUM_VERTICES_DAYS_PER_YEAR);
	IssmDouble* dailydsradiation=xNew<IssmDouble>(NUM_VERTICES_DAYS_PER_YEAR);
	IssmDouble* dailywindspeed=xNew<IssmDouble>(NUM_VERTICES_DAYS_PER_YEAR);
	IssmDouble* dailypressure=xNew<IssmDouble>(NUM_VERTICES_DAYS_PER_YEAR);
	IssmDouble* dailyairdensity=xNew<IssmDouble>(NUM_VERTICES_DAYS_PER_YEAR);
	IssmDouble* dailyairhumidity=xNew<IssmDouble>(NUM_VERTICES_DAYS_PER_YEAR);
	IssmDouble* dailytemperature=xNew<IssmDouble>(NUM_VERTICES_DAYS_PER_YEAR);
	// daily outputs
	IssmDouble* tsurf_out=xNew<IssmDouble>(NUM_VERTICES); memset(tsurf_out, 0., NUM_VERTICES*sizeof(IssmDouble));
	IssmDouble* smb_out=xNew<IssmDouble>(NUM_VERTICES); memset(smb_out, 0., NUM_VERTICES*sizeof(IssmDouble));
	IssmDouble* saccu_out=xNew<IssmDouble>(NUM_VERTICES); memset(saccu_out, 0., NUM_VERTICES*sizeof(IssmDouble));
	IssmDouble* smelt_out=xNew<IssmDouble>(NUM_VERTICES); memset(smelt_out, 0., NUM_VERTICES*sizeof(IssmDouble));

	IssmDouble rho_water,rho_ice,desfac,rlaps,rdl;
	IssmDouble time,yts,time_yr;

	/* Get time: */
	this->parameters->FindParam(&time,TimeEnum);
	this->parameters->FindParam(&yts,ConstantsYtsEnum);
	time_yr=floor(time/yts)*yts;

	/*Get material parameters :*/
	rho_water=this->FindParam(MaterialsRhoSeawaterEnum);
	rho_ice=this->FindParam(MaterialsRhoIceEnum);
	desfac=this->FindParam(SmbDesfacEnum);
	rlaps=this->FindParam(SmbRlapsEnum);
	rdl=this->FindParam(SmbRdlEnum);

	/* Recover info at the vertices: */
	GetInputListOnVertices(&s[0],SurfaceEnum);
	GetInputListOnVertices(&s0gcm[0],SmbS0gcmEnum);

	/* loop over vertices and days */ //FIXME account for leap years (365 -> 366)
	Gauss* gauss=this->NewGauss();
	for (int iday = 0; iday < 365; iday++){
		/* Retrieve inputs: */
		Input* dailysnowfall_input    = this->GetInput(SmbDailysnowfallEnum,time_yr+(iday+1)/365.*yts); _assert_(dailysnowfall_input);
		Input* dailyrainfall_input    = this->GetInput(SmbDailyrainfallEnum,time_yr+(iday+1)/365.*yts); _assert_(dailyrainfall_input);
		Input* dailydlradiation_input = this->GetInput(SmbDailydlradiationEnum,time_yr+(iday+1)/365.*yts); _assert_(dailydlradiation_input);
		Input* dailydsradiation_input = this->GetInput(SmbDailydsradiationEnum,time_yr+(iday+1)/365.*yts); _assert_(dailydsradiation_input);
		Input* dailywindspeed_input   = this->GetInput(SmbDailywindspeedEnum,time_yr+(iday+1)/365.*yts); _assert_(dailywindspeed_input);
		Input* dailypressure_input    = this->GetInput(SmbDailypressureEnum,time_yr+(iday+1)/365.*yts); _assert_(dailypressure_input);
		Input* dailyairdensity_input  = this->GetInput(SmbDailyairdensityEnum,time_yr+(iday+1)/365.*yts); _assert_(dailyairdensity_input);
		Input* dailyairhumidity_input = this->GetInput(SmbDailyairhumidityEnum,time_yr+(iday+1)/365.*yts); _assert_(dailyairhumidity_input);
		Input* dailytemperature_input = this->GetInput(SmbDailytemperatureEnum,time_yr+(iday+1)/365.*yts); _assert_(dailytemperature_input);

		for(int iv=0;iv<NUM_VERTICES;iv++){
			gauss->GaussVertex(iv);
			/* get forcing */
			dailyrainfall_input->GetInputValue(&dailyrainfall[iv*365+iday],gauss);
			dailysnowfall_input->GetInputValue(&dailysnowfall[iv*365+iday],gauss);
			dailydlradiation_input->GetInputValue(&dailydlradiation[iv*365+iday],gauss);
			dailydsradiation_input->GetInputValue(&dailydsradiation[iv*365+iday],gauss);
			dailywindspeed_input->GetInputValue(&dailywindspeed[iv*365+iday],gauss);
			dailypressure_input->GetInputValue(&dailypressure[iv*365+iday],gauss);
			dailyairdensity_input->GetInputValue(&dailyairdensity[iv*365+iday],gauss);
			dailyairhumidity_input->GetInputValue(&dailyairhumidity[iv*365+iday],gauss);
			dailytemperature_input->GetInputValue(&dailytemperature[iv*365+iday],gauss);

			/* Surface temperature correction */
			st[iv]=(s[iv]-s0gcm[iv])/1000.;
			dailytemperature[iv*365+iday]=dailytemperature[iv*365+iday]-rlaps *st[iv];

			/* Precipitation correction (Vizcaino et al. 2010) */
			if (s0gcm[iv] < 2000.0) {
				dailysnowfall[iv*365+iday] = dailysnowfall[iv*365+iday]*exp(desfac*(max(s[iv],2000.0)-2000.0));
				dailyrainfall[iv*365+iday] = dailyrainfall[iv*365+iday]*exp(desfac*(max(s[iv],2000.0)-2000.0));
			}else{
				dailysnowfall[iv*365+iday] = dailysnowfall[iv*365+iday]*exp(desfac*(max(s[iv],2000.0)-s0gcm[iv]));
				dailyrainfall[iv*365+iday] = dailyrainfall[iv*365+iday]*exp(desfac*(max(s[iv],2000.0)-s0gcm[iv]));
			}

			/* downward longwave radiation correction (Marty et al. 2002) */
			st[iv]=(s[iv]-s0gcm[iv])/1000.;
			dailydlradiation[iv*365+iday]=dailydlradiation[iv*365+iday]+rdl*st[iv];
		}
	}

	for (int iv = 0; iv<NUM_VERTICES; iv++){
		/* call semic */
		run_semic_(&dailysnowfall[iv*365], &dailyrainfall[iv*365], &dailydsradiation[iv*365], &dailydlradiation[iv*365],
					&dailywindspeed[iv*365], &dailypressure[iv*365], &dailyairdensity[iv*365], &dailyairhumidity[iv*365], &dailytemperature[iv*365],
					&tsurf_out[iv], &smb_out[iv], &saccu_out[iv], &smelt_out[iv]);
	}

	switch(this->ObjectEnum()){
		case TriaEnum:
			this->AddInput(TemperatureSEMICEnum,&tsurf_out[0],P1Enum); // TODO add TemperatureSEMICEnum to EnumDefinitions
			this->AddInput(SmbMassBalanceEnum,&smb_out[0],P1Enum);
			this->AddInput(SmbAccumulationEnum,&saccu_out[0],P1Enum);
			this->AddInput(SmbMeltEnum,&smelt_out[0],P1Enum);
			break;
		case PentaEnum:
			// TODO
			break;
		case TetraEnum:
			// TODO
			break;
		default: _error_("Not implemented yet");
	}

	/*clean-up*/
	delete gauss;
	xDelete<IssmDouble>(dailysnowfall);
	xDelete<IssmDouble>(dailyrainfall);
	xDelete<IssmDouble>(dailydlradiation);
	xDelete<IssmDouble>(dailydsradiation);
	xDelete<IssmDouble>(dailywindspeed);
	xDelete<IssmDouble>(dailypressure);
	xDelete<IssmDouble>(dailyairdensity);
	xDelete<IssmDouble>(dailyairhumidity);
	xDelete<IssmDouble>(dailypressure);
	xDelete<IssmDouble>(dailytemperature);
	xDelete<IssmDouble>(smb_out);
	xDelete<IssmDouble>(smelt_out);
	xDelete<IssmDouble>(saccu_out);
	xDelete<IssmDouble>(tsurf_out);
	xDelete<IssmDouble>(s);
	xDelete<IssmDouble>(st);
	xDelete<IssmDouble>(s0gcm);
}
/*}}}*/
void       Element::SmbSemicTransient(){/*{{{*/

	bool isverbose=VerboseSmb();
	if(isverbose && this->Sid()==0){
		_printf0_("smb core: initialize.\n");
	}
	/*only compute SMB at the surface: */
	if (!IsOnSurface()) return;

	const int NUM_VERTICES 					= this->GetNumberOfVertices();

	// daily forcing inputs
	IssmDouble* dailyrainfall   =xNew<IssmDouble>(NUM_VERTICES);
	IssmDouble* dailysnowfall   =xNew<IssmDouble>(NUM_VERTICES);
	IssmDouble* dailydlradiation=xNew<IssmDouble>(NUM_VERTICES);
	IssmDouble* dailydsradiation=xNew<IssmDouble>(NUM_VERTICES);
	IssmDouble* dailywindspeed  =xNew<IssmDouble>(NUM_VERTICES);
	IssmDouble* dailypressure   =xNew<IssmDouble>(NUM_VERTICES);
	IssmDouble* dailyairdensity =xNew<IssmDouble>(NUM_VERTICES);
	IssmDouble* dailyairhumidity=xNew<IssmDouble>(NUM_VERTICES);
	IssmDouble* dailytemperature=xNew<IssmDouble>(NUM_VERTICES);

	// inputs: geometry
	IssmDouble* s=xNew<IssmDouble>(NUM_VERTICES);
	IssmDouble* s0gcm=xNew<IssmDouble>(NUM_VERTICES);
	IssmDouble* st=xNew<IssmDouble>(NUM_VERTICES);

	// inputs
	IssmDouble* tsurf_in        =xNew<IssmDouble>(NUM_VERTICES); 
	IssmDouble* mask_in         =xNew<IssmDouble>(NUM_VERTICES); 
	IssmDouble* Tamp_in         =xNew<IssmDouble>(NUM_VERTICES); 
	IssmDouble* albedo_in       =xNew<IssmDouble>(NUM_VERTICES); 
	IssmDouble* albedo_snow_in  =xNew<IssmDouble>(NUM_VERTICES); 
	IssmDouble* hice_in         =xNew<IssmDouble>(NUM_VERTICES); 
	IssmDouble* hsnow_in        =xNew<IssmDouble>(NUM_VERTICES); 
	IssmDouble* qmr_in          =xNew<IssmDouble>(NUM_VERTICES); 

	/* outputs 
	runoff - runoff calculated melting + rainfall - refreezing 
	subl   - sublimation
	 * */
	IssmDouble* tsurf_out  =xNew<IssmDouble>(NUM_VERTICES); memset(tsurf_out, 0., NUM_VERTICES*sizeof(IssmDouble));
	IssmDouble* smb_out    =xNew<IssmDouble>(NUM_VERTICES); memset(smb_out, 0., NUM_VERTICES*sizeof(IssmDouble));
	IssmDouble* smbi_out   =xNew<IssmDouble>(NUM_VERTICES); memset(smb_out, 0., NUM_VERTICES*sizeof(IssmDouble));
	IssmDouble* smbs_out   =xNew<IssmDouble>(NUM_VERTICES); memset(smb_out, 0., NUM_VERTICES*sizeof(IssmDouble));
	IssmDouble* saccu_out  =xNew<IssmDouble>(NUM_VERTICES); memset(saccu_out, 0., NUM_VERTICES*sizeof(IssmDouble));
	IssmDouble* smelt_out  =xNew<IssmDouble>(NUM_VERTICES); memset(smelt_out, 0., NUM_VERTICES*sizeof(IssmDouble));
	IssmDouble* refr_out  =xNew<IssmDouble>(NUM_VERTICES); memset(refr_out, 0., NUM_VERTICES*sizeof(IssmDouble));
	IssmDouble* albedo_out =xNew<IssmDouble>(NUM_VERTICES); memset(albedo_out, 0., NUM_VERTICES*sizeof(IssmDouble));
	IssmDouble* albedo_snow_out =xNew<IssmDouble>(NUM_VERTICES); memset(albedo_snow_out, 0., NUM_VERTICES*sizeof(IssmDouble));
	IssmDouble* hsnow_out   =xNew<IssmDouble>(NUM_VERTICES); memset(hsnow_out, 0., NUM_VERTICES*sizeof(IssmDouble));
	IssmDouble* hice_out    =xNew<IssmDouble>(NUM_VERTICES); memset(hice_out, 0., NUM_VERTICES*sizeof(IssmDouble));
	IssmDouble* qmr_out     =xNew<IssmDouble>(NUM_VERTICES); memset(qmr_out, 0., NUM_VERTICES*sizeof(IssmDouble)); 
	IssmDouble* runoff_out  =xNew<IssmDouble>(NUM_VERTICES); memset(runoff_out, 0., NUM_VERTICES*sizeof(IssmDouble)); 
	IssmDouble* subl_out  =xNew<IssmDouble>(NUM_VERTICES); memset(subl_out, 0., NUM_VERTICES*sizeof(IssmDouble)); 

	IssmDouble rho_water,rho_ice,desfac,desfacElev,rlaps,rdl;
	IssmDouble alb_smax, alb_smin, albi, albl;
	IssmDouble hcrit, rcrit; // parameters for ? and refreezing.
	int alb_scheme;
	// albedo parameters - slatter
	IssmDouble tmin, tmax;
	// albedo parameters - isba
	IssmDouble mcrit, tau_a, tau_f, wcrit;
	// albedo parameters - alex
	IssmDouble tmid, afac;

	IssmDouble tstart, time,yts,time_yr,dt;

	int isdesertification, isLWDcorrect;

	/* Get time: */
	this->parameters->FindParam(&time,TimeEnum);
	this->parameters->FindParam(&dt,TimesteppingTimeStepEnum);
	this->parameters->FindParam(&yts,ConstantsYtsEnum);
	this->parameters->FindParam(&tstart,TimesteppingStartTimeEnum);
	time_yr=floor(time/yts)*yts;
	//dt = dt * yts;

	/*Get material parameters :*/
	rho_water=this->FindParam(MaterialsRhoFreshwaterEnum);
	rho_ice=this->FindParam(MaterialsRhoIceEnum);
	desfac=this->FindParam(SmbDesfacEnum);
	desfacElev=this->FindParam(SmbDesfacElevEnum);
	rlaps=this->FindParam(SmbRlapsEnum);
	rdl=this->FindParam(SmbRdlEnum);

	this->FindParam(&alb_scheme,SmbAlbedoSchemeEnum);
	this->FindParam(&hcrit,SmbSemicHcritEnum);
	this->FindParam(&rcrit,SmbSemicRcritEnum);
	alb_smax=this->FindParam(SmbAlbedoSnowMaxEnum);
	alb_smin=this->FindParam(SmbAlbedoSnowMinEnum);
	albi=this->FindParam(SmbAlbedoIceEnum);
	albl=this->FindParam(SmbAlbedoLandEnum);

	// albedo parameters
	this->FindParam(&tmid,SmbSemicTmidEnum);
	this->FindParam(&tmin,SmbSemicTminEnum);
	this->FindParam(&tmax,SmbSemicTmaxEnum);
	this->FindParam(&mcrit,SmbSemicMcritEnum);
	this->FindParam(&wcrit,SmbSemicWcritEnum);
	this->FindParam(&tau_a,SmbSemicTauAEnum);
	this->FindParam(&tau_f,SmbSemicTauFEnum);
	this->FindParam(&afac,SmbSemicAfacEnum);

	/* Recover info at the vertices: */
	GetInputListOnVertices(&s[0],SurfaceEnum);
	GetInputListOnVertices(&s0gcm[0],SmbS0gcmEnum);

	/* Get specific parameter options */
	this->FindParam(&isdesertification, SmbSemicIsDesertificationEnum);
	this->FindParam(&isLWDcorrect, SmbSemicIsLWDcorrectEnum);

	if(isverbose && this->Sid()==0){
		_printf0_("smb core: allocate inputs.\n");
		_printf0_("smb core: time_yr  : " << time_yr/yts <<"\n");
		_printf0_("smb core: time     : " << time <<"\n");
		_printf0_("smb core: dt       : " << dt <<"\n");
	}
	/* loop over vertices and days */
	Gauss* gauss=this->NewGauss();
	/* Retrieve inputs: */
	Input* dailysnowfall_input    = this->GetInput(SmbDailysnowfallEnum,time); _assert_(dailysnowfall_input);
	Input* dailyrainfall_input    = this->GetInput(SmbDailyrainfallEnum,time); _assert_(dailyrainfall_input);
	Input* dailydlradiation_input = this->GetInput(SmbDailydlradiationEnum,time); _assert_(dailydlradiation_input);
	Input* dailydsradiation_input = this->GetInput(SmbDailydsradiationEnum,time); _assert_(dailydsradiation_input);
	Input* dailywindspeed_input   = this->GetInput(SmbDailywindspeedEnum,time); _assert_(dailywindspeed_input);
	Input* dailypressure_input    = this->GetInput(SmbDailypressureEnum,time); _assert_(dailypressure_input);
	Input* dailyairdensity_input  = this->GetInput(SmbDailyairdensityEnum,time); _assert_(dailyairdensity_input);
	Input* dailyairhumidity_input = this->GetInput(SmbDailyairhumidityEnum,time); _assert_(dailyairhumidity_input);
	Input* dailytemperature_input = this->GetInput(SmbDailytemperatureEnum,time); _assert_(dailytemperature_input);

	/*temporal Enum depending on time*/
	int enum_temp       =TemperatureSEMICEnum;
	int enum_hice       =SmbHIceEnum;
	int enum_hsnow      =SmbHSnowEnum;
	int enum_albedo     =SmbAlbedoEnum;
	int enum_albedo_snow=SmbAlbedoSnowEnum;
	int enum_qmr        =SmbSemicQmrEnum;
	if (tstart+dt == time) {
		/* Load inital value at first time step*/
		enum_temp=TemperatureEnum;
		enum_hice=SmbHIceInitEnum;
		enum_hsnow=SmbHSnowInitEnum;
		enum_albedo=SmbAlbedoInitEnum;
		enum_albedo_snow=SmbAlbedoSnowInitEnum;
		enum_qmr        =SmbSemicQmrInitEnum;
	} 
	//if(isverbose && this->Sid()==0)_printf0_("smb core: assign temp.\n");
	Input* tsurf_input       = this->GetInput(enum_temp); _assert_(tsurf_in);
	//if(isverbose && this->Sid()==0)_printf0_("smb core: assign mask.\n");
	Input* mask_input        = this->GetInput(SmbMaskEnum); _assert_(mask_input);
	//if(isverbose && this->Sid()==0)_printf0_("smb core: assign Tamp.\n");
	Input* Tamp_input        = this->GetInput(SmbTampEnum); _assert_(Tamp_input);
	//if(isverbose && this->Sid()==0)_printf0_("smb core: assign albedo.\n");
	Input* albedo_input      = this->GetInput(enum_albedo); _assert_(albedo_input);
	Input* albedo_snow_input = this->GetInput(enum_albedo_snow); _assert_(albedo_snow_input);
	Input* hice_input        = this->GetInput(enum_hice); _assert_(hice_input);
	Input* hsnow_input       = this->GetInput(enum_hsnow); _assert_(hsnow_input);
	Input* qmr_input         = this->GetInput(enum_qmr); _assert_(qmr_input);

	if(isverbose && this->Sid()==0)_printf0_("smb core: assign inputs done....\n");
	for(int iv=0;iv<NUM_VERTICES;iv++){
		gauss->GaussVertex(iv);
		/* get forcing */
		dailyrainfall_input->GetInputValue(&dailyrainfall[iv],gauss);
		dailysnowfall_input->GetInputValue(&dailysnowfall[iv],gauss);
		dailydlradiation_input->GetInputValue(&dailydlradiation[iv],gauss);
		dailydsradiation_input->GetInputValue(&dailydsradiation[iv],gauss);
		dailywindspeed_input->GetInputValue(&dailywindspeed[iv],gauss);
		dailypressure_input->GetInputValue(&dailypressure[iv],gauss);
		dailyairdensity_input->GetInputValue(&dailyairdensity[iv],gauss);
		dailyairhumidity_input->GetInputValue(&dailyairhumidity[iv],gauss);
		dailytemperature_input->GetInputValue(&dailytemperature[iv],gauss);
		tsurf_input->GetInputValue(&tsurf_in[iv],gauss);

		/* Get Albedo information */
		albedo_input->GetInputValue(&albedo_in[iv],gauss);
		albedo_snow_input->GetInputValue(&albedo_snow_in[iv],gauss);
		mask_input->GetInputValue(&mask_in[iv],gauss);
		Tamp_input->GetInputValue(&Tamp_in[iv],gauss);

		hsnow_input->GetInputValue(&hsnow_in[iv],gauss);
		hice_input->GetInputValue(&hice_in[iv],gauss);
		qmr_input->GetInputValue(&qmr_in[iv],gauss);

		/* Surface temperature correction */
		st[iv]=(s[iv]-s0gcm[iv])/1000.;
		dailytemperature[iv]=dailytemperature[iv]-rlaps *st[iv];

		/* Precipitation correction (Vizcaino et al. 2010) */
		if (isdesertification == 1){
			if (s0gcm[iv] < desfacElev) {
				dailysnowfall[iv] = dailysnowfall[iv]*exp(desfac*(max(s[iv],desfacElev)-desfacElev));
				dailyrainfall[iv] = dailyrainfall[iv]*exp(desfac*(max(s[iv],desfacElev)-desfacElev));
			}else{
				dailysnowfall[iv] = dailysnowfall[iv]*exp(desfac*(max(s[iv],desfacElev)-s0gcm[iv]));
				dailyrainfall[iv] = dailyrainfall[iv]*exp(desfac*(max(s[iv],desfacElev)-s0gcm[iv]));
			}
		}

		/* downward longwave radiation correction (Marty et al. 2002) */
      /* Unit of "md.smb.rdl" is defined in W m-2 km-1 */
		if (isLWDcorrect == 1){
			st[iv]=(s[iv]-s0gcm[iv])/1000.; /* unit in km */
			dailydlradiation[iv]=dailydlradiation[iv]-rdl*st[iv];
		}
	}
	if(isverbose && this->Sid()==0){
		_printf0_("smb core: assign tsurf_in        :" << tsurf_in[0] << "\n");
		_printf0_("smb core: assign dailytemperature:" << dailytemperature[0] << "\n");
		_printf0_("smb core: assign hsnow           :" << hsnow_in[0] << "\n");
		_printf0_("smb core: assign hice            :" << hice_in[0] << "\n");
		_printf0_("smb core: assign mask            :" << mask_in[0] << "\n");
		_printf0_("smb core: assign Tamp            :" << Tamp_in[0] << "\n");
		_printf0_("smb core: assign albedo          :" << albedo_in[0] << "\n");
		_printf0_("smb core: assign albedo_snow     :" << albedo_snow_in[0] << "\n");
		_printf0_("smb core: assign albedo_scheme   :" << alb_scheme  << "\n");
		_printf0_("smb core: assign qmr             :" << qmr_in[0]  << "\n");
	}

	if(isverbose && this->Sid()==0)_printf0_("smb core: call run_semic_transient module.\n");
	/* call semic */
	int nx=NUM_VERTICES, ntime=1, nloop=1;
	bool semic_verbose=false; //VerboseSmb();
	run_semic_transient_(&nx, &ntime, &nloop,
			dailysnowfall,  dailyrainfall, dailydsradiation, dailydlradiation,
			dailywindspeed, dailypressure, dailyairdensity,  dailyairhumidity, dailytemperature, tsurf_in, qmr_in, 
			&dt,
			&hcrit, &rcrit, 
			mask_in, hice_in, hsnow_in, 
			albedo_in, albedo_snow_in,
			&alb_scheme, &alb_smax, &alb_smin, &albi, &albl,
			Tamp_in,
			&tmin, &tmax, &tmid, &mcrit, &wcrit, &tau_a, &tau_f, &afac, &semic_verbose,
			tsurf_out, smb_out, smbi_out, smbs_out, saccu_out, smelt_out, refr_out, albedo_out, albedo_snow_out, hsnow_out, hice_out, qmr_out, runoff_out, subl_out);

	for (int iv = 0; iv<NUM_VERTICES; iv++){
		/* 
		 unit conversion: water -> ice
		 w.e. : water equivalenet.
		 */
		smb_out[iv]  = smb_out[iv]*rho_water/rho_ice;      // w.e. m/sec -> ice m/yr
		smbi_out[iv] = smbi_out[iv]*rho_water/rho_ice*yts; // w.e. m/sec -> ice m/yr
		smbs_out[iv] = smbs_out[iv]*rho_water/rho_ice*yts; // w.e. m/sec -> ice m/yr
		saccu_out[iv] = saccu_out[iv]*rho_water/rho_ice*yts; // w.e. m/sec -> ice m/yr
		smelt_out[iv] = smelt_out[iv]*rho_water/rho_ice; // w.e. m/sec -> ice m/yr
		refr_out[iv]  = refr_out[iv]*rho_water/rho_ice;  // w.e. m/sec -> ice m/yr
		runoff_out[iv]= runoff_out[iv]*rho_water/rho_ice; // w.e. m/sec -> ice m/yr
	}

	if(isverbose && this->Sid()==0){
		_printf0_("smb core: tsurf_out " << tsurf_out[0] << " " << tsurf_out[1] << " " << tsurf_out[2] << "\n");
		_printf0_("smb core: hice_out  " << hice_out[0] << " " <<  hice_out[1] << " " << hice_out[2] << "\n");
		_printf0_("smb core: hsnow_out " << hsnow_out[0] << "\n");
		_printf0_("smb core: smb_ice   " << smbi_out[0]*yts << "\n");
		_printf0_("smb core: smb_ice   " << albedo_out[0] <<" "<<albedo_out[1] << " " << albedo_out[2] << "\n");
	}

	switch(this->ObjectEnum()){
		case TriaEnum:
			this->AddInput(TemperatureSEMICEnum,  &tsurf_out[0],P1DGEnum);
			// SMBout = SMB_ice + SMB_snow values.
			//this->AddInput(SmbMassBalanceTotalEnum,&smb_out[0],P1DGEnum);
			// water equivalent SMB ice to ice equivalent.
			this->AddInput(SmbMassBalanceEnum,     &smb_out[0],P1DGEnum);
			this->AddInput(SmbMassBalanceIceEnum,  &smbi_out[0],P1DGEnum);
			this->AddInput(SmbMassBalanceSnowEnum, &smbs_out[0],P1DGEnum);
			//this->AddInput(SmbMassBalanceSnowEnum,&smbs_out[0],P1DGEnum);
			// saccu - accumulation of snow.
			this->AddInput(SmbAccumulationEnum,&saccu_out[0],P1DGEnum);
			// smelt 
			this->AddInput(SmbMeltEnum,        &smelt_out[0],P1DGEnum);
			this->AddInput(SmbRefreezeEnum,    &refr_out[0],P1DGEnum);
			this->AddInput(SmbAlbedoEnum,      &albedo_out[0],P1DGEnum);
			this->AddInput(SmbAlbedoSnowEnum,  &albedo_snow_out[0],P1DGEnum);
			this->AddInput(SmbHSnowEnum,       &hsnow_out[0],P1DGEnum);
			this->AddInput(SmbHIceEnum,        &hice_out[0],P1DGEnum);
			this->AddInput(SmbSemicQmrEnum,    &qmr_out[0],P1DGEnum);
			this->AddInput(SmbRunoffEnum,      &runoff_out[0],P1DGEnum);
			this->AddInput(SmbEvaporationEnum, &subl_out[0],P1DGEnum);
			break;
		case PentaEnum:
			// TODO
			break;
		case TetraEnum:
			// TODO
			break;
		default: _error_("Not implemented yet");
	}

	/*clean-up {{{*/
	delete gauss;
	xDelete<IssmDouble>(dailysnowfall);
	xDelete<IssmDouble>(dailyrainfall);
	xDelete<IssmDouble>(dailydlradiation);
	xDelete<IssmDouble>(dailydsradiation);
	xDelete<IssmDouble>(dailywindspeed);
	xDelete<IssmDouble>(dailypressure);
	xDelete<IssmDouble>(dailyairdensity);
	xDelete<IssmDouble>(dailyairhumidity);
	xDelete<IssmDouble>(dailypressure);
	xDelete<IssmDouble>(dailytemperature);

	/*for outputs*/
	xDelete<IssmDouble>(tsurf_out);
	xDelete<IssmDouble>(smb_out);
	xDelete<IssmDouble>(smbi_out);
	xDelete<IssmDouble>(smbs_out);
	xDelete<IssmDouble>(saccu_out);
	xDelete<IssmDouble>(smelt_out);
	xDelete<IssmDouble>(refr_out);
	xDelete<IssmDouble>(runoff_out);
	xDelete<IssmDouble>(subl_out);
	xDelete<IssmDouble>(albedo_out);
	xDelete<IssmDouble>(albedo_snow_out);
	xDelete<IssmDouble>(hsnow_out);
	xDelete<IssmDouble>(hice_out);
	xDelete<IssmDouble>(qmr_out);

	/*for inputs*/
	xDelete<IssmDouble>(hsnow_in);
	xDelete<IssmDouble>(hice_in);
	xDelete<IssmDouble>(mask_in);
	xDelete<IssmDouble>(Tamp_in);
	xDelete<IssmDouble>(albedo_in);
	xDelete<IssmDouble>(albedo_snow_in);
	xDelete<IssmDouble>(tsurf_in);
	xDelete<IssmDouble>(qmr_in);

	/* for inputs:geometry */
	xDelete<IssmDouble>(s);
	xDelete<IssmDouble>(st);
	xDelete<IssmDouble>(s0gcm);
	/*}}}*/
}
/*}}}*/
#endif // _HAVE_SEMIC_
int        Element::Sid(){/*{{{*/

	return this->sid;

}
/*}}}*/
void       Element::SmbGemb(IssmDouble timeinputs, int count, int steps){/*{{{*/

	/*only compute SMB at the surface: */
	if (!IsOnSurface()) return;

	/*Intermediary variables: {{{*/
	bool       isinitialized;
	IssmDouble zTop=0.0;
	IssmDouble dzTop=0.0;
	IssmDouble zMax=0.0;
	IssmDouble zMin=0.0;
	IssmDouble zY=0.0;
	IssmDouble dzMin=0.0;
	IssmDouble Tmean=0.0;
	IssmDouble Vmean=0.0;
	IssmDouble C=0.0;
	IssmDouble Tz,Vz=0.0;
	IssmDouble yts;
	IssmDouble Ta=0.0;
	IssmDouble V=0.0;
	IssmDouble dlw=0.0;
	IssmDouble dsw=0.0;
	IssmDouble dswdiff=0.0;
	IssmDouble P=0.0;
	IssmDouble eAir=0.0;
	IssmDouble pAir=0.0;
	IssmDouble teValue=1.0;
	IssmDouble aValue=0.0;
	IssmDouble dulwrfValue=0.0;
	IssmDouble szaValue=0.0;
	IssmDouble cotValue=0.0;
	IssmDouble ccsnowValue=0.0;
	IssmDouble cciceValue=0.0;
	IssmDouble dt,time,smb_dt;
	IssmDouble currentsurface;
	int        Mappedpoint=0;
	int        aIdx=0;
	int        eIdx=0;
	int        tcIdx=0;
	int        denIdx=0;
	int        dsnowIdx=0;
	int        swIdx=0;
	int        N=0;
	IssmDouble cldFrac,t0wet, t0dry, K;
	IssmDouble lhf=0.0;
	IssmDouble shf=0.0;
	IssmDouble dayEC=0.0;
	IssmDouble initMass=0.0;
   IssmDouble sumR=0.0;
	IssmDouble sumF=0.0;
	IssmDouble sumM=0.0;
	IssmDouble sumMsurf=0.0;
	IssmDouble sumEC=0.0;
	IssmDouble sumP=0.0;
	IssmDouble sumRa=0.0;
	IssmDouble sumW=0.0;
	IssmDouble sumMassAdd=0.0;
	IssmDouble fac=0.0;
	IssmDouble sumMass=0.0;
	IssmDouble sumH=0.0;
	IssmDouble T0m=0.0;
	IssmDouble T10m=0.0;
	IssmDouble T30m=0.0;
	IssmDouble T50m=0.0;
	IssmDouble dMass=0.0;
	IssmDouble accsumR=0.0;
	IssmDouble accsumF=0.0;
	IssmDouble accsumM=0.0;
	IssmDouble accsumSMB=0.0;
	IssmDouble accsumEC=0.0;
	IssmDouble accsumP=0.0;
	IssmDouble accsumRa=0.0;
	bool isgraingrowth,isalbedo,isshortwave,isthermal,isaccumulation,ismelt,isdensification,isturbulentflux,ismappedforcing;
	bool isconstrainsurfaceT=false;
	bool isdeltaLWup=false;
	IssmDouble init_scaling=0.0;
	IssmDouble thermo_scaling=1.0;
	IssmDouble adThresh=1023.0;
	IssmDouble teThresh=10;
	IssmDouble tlapse=0.0;
	IssmDouble dlwlapse=0.0;
	/*}}}*/
	/*Output variables:{{{ */
	IssmDouble* dz=NULL;
	IssmDouble* d = NULL;
	IssmDouble* re = NULL;
	IssmDouble* gdn = NULL;
	IssmDouble* gsp = NULL;
	IssmDouble  EC = 0.0;
	IssmDouble  ulw = 0.0;
	IssmDouble  netSW=0.0;
	IssmDouble  netLW=0.0;
	IssmDouble  meanULW=0.0;
	IssmDouble  meanLHF=0.0;
	IssmDouble  meanSHF=0.0;
	IssmDouble* W = NULL;
	IssmDouble* a = NULL;
	IssmDouble* adiff = NULL;
	IssmDouble* swf=NULL;
	IssmDouble* T = NULL;
	IssmDouble  T_bottom = 0.0;
	IssmDouble  M = 0.0;
	IssmDouble  Msurf = 0.0;
	IssmDouble  R = 0.0;
	IssmDouble  F = 0.0;
	IssmDouble  Ra = 0.0;
	IssmDouble  mAdd = 0.0;
	IssmDouble  dz_add = 0.0;
	IssmDouble* dzini=NULL;
	IssmDouble* dini = NULL;
	IssmDouble* reini = NULL;
	IssmDouble* gdnini = NULL;
	IssmDouble* gspini = NULL;
	IssmDouble* Wini = NULL;
	IssmDouble* aini = NULL;
	IssmDouble* adiffini = NULL;
	IssmDouble* Tini = NULL;
	IssmDouble* ECsub=NULL;
	IssmDouble* SMBsub=NULL;
	IssmDouble* Msub=NULL;
	IssmDouble* Rsub=NULL;
	IssmDouble* Fsub=NULL;
	IssmDouble* Rasub=NULL;
	IssmDouble* Psub=NULL;
	IssmDouble* FACsub=NULL;
	IssmDouble* ECdsub=NULL;
	IssmDouble* SMBdsub=NULL;
	IssmDouble* Mdsub=NULL;
	IssmDouble* Rdsub=NULL;
	IssmDouble* Fdsub=NULL;
	IssmDouble* Radsub=NULL;
	IssmDouble* Pdsub=NULL;
	int         m=0;
	/*}}}*/

	/*Retrieve material properties and parameters:{{{ */
	IssmDouble rho_ice   = FindParam(MaterialsRhoIceEnum);
	IssmDouble rho_water = FindParam(MaterialsRhoFreshwaterEnum);
	IssmDouble aSnow     = parameters->FindParam(SmbASnowEnum);
	IssmDouble aIce      = parameters->FindParam(SmbAIceEnum);
	parameters->FindParam(&time,TimeEnum);                        /*transient core time at which we run the smb core*/
	parameters->FindParam(&dt,TimesteppingTimeStepEnum);          /*transient core time step*/
	parameters->FindParam(&yts,ConstantsYtsEnum);
	parameters->FindParam(&smb_dt,SmbDtEnum);                     /*time period for the smb solution,  usually smaller than the glaciological dt*/
	parameters->FindParam(&aIdx,SmbAIdxEnum);
	parameters->FindParam(&eIdx,SmbEIdxEnum);
	parameters->FindParam(&tcIdx,SmbTcIdxEnum);
	parameters->FindParam(&denIdx,SmbDenIdxEnum);
	parameters->FindParam(&swIdx,SmbSwIdxEnum);
	parameters->FindParam(&dsnowIdx,SmbDsnowIdxEnum);
	parameters->FindParam(&cldFrac,SmbCldFracEnum);
	parameters->FindParam(&t0wet,SmbT0wetEnum);
	parameters->FindParam(&t0dry,SmbT0dryEnum);
	parameters->FindParam(&K,SmbKEnum);
	parameters->FindParam(&isgraingrowth,SmbIsgraingrowthEnum);
	parameters->FindParam(&isalbedo,SmbIsalbedoEnum);
	parameters->FindParam(&isshortwave,SmbIsshortwaveEnum);
	parameters->FindParam(&isthermal,SmbIsthermalEnum);
	parameters->FindParam(&isaccumulation,SmbIsaccumulationEnum);
	parameters->FindParam(&ismelt,SmbIsmeltEnum);
	parameters->FindParam(&isdensification,SmbIsdensificationEnum);
	parameters->FindParam(&isturbulentflux,SmbIsturbulentfluxEnum);
	parameters->FindParam(&isconstrainsurfaceT,SmbIsconstrainsurfaceTEnum);
	parameters->FindParam(&isdeltaLWup,SmbIsdeltaLWupEnum);
	parameters->FindParam(&init_scaling,SmbInitDensityScalingEnum);
	parameters->FindParam(&thermo_scaling,SmbThermoDeltaTScalingEnum);
	parameters->FindParam(&adThresh,SmbAdThreshEnum);
	parameters->FindParam(&teThresh,SmbTeThreshEnum);
	parameters->FindParam(&ismappedforcing,SmbIsmappedforcingEnum);
	/*}}}*/
	/*Retrieve inputs: {{{*/
	Input *zTop_input          = this->GetInput(SmbZTopEnum);         _assert_(zTop_input);
	Input *dzTop_input         = this->GetInput(SmbDzTopEnum);        _assert_(dzTop_input);
	Input *dzMin_input         = this->GetInput(SmbDzMinEnum);        _assert_(dzMin_input);
	Input *zMax_input          = this->GetInput(SmbZMaxEnum);         _assert_(zMax_input);
	Input *zMin_input          = this->GetInput(SmbZMinEnum);         _assert_(zMin_input);
	Input *zY_input            = this->GetInput(SmbZYEnum);           _assert_(zY_input);
	Input *EC_input            = NULL;

	/*Retrieve input values:*/
	Gauss* gauss=this->NewGauss(1); gauss->GaussPoint(0);

	this->GetInputValue(&isinitialized,SmbIsInitializedEnum);
	zTop_input->GetInputValue(&zTop,gauss);
	dzTop_input->GetInputValue(&dzTop,gauss);
	dzMin_input->GetInputValue(&dzMin,gauss);
	zMax_input->GetInputValue(&zMax,gauss);
	zMin_input->GetInputValue(&zMin,gauss);
	zY_input->GetInputValue(&zY,gauss);

	if (!ismappedforcing) {
		Input *Tmean_input         = this->GetInput(SmbTmeanEnum);        _assert_(Tmean_input);
		Input *Vmean_input         = this->GetInput(SmbVmeanEnum);        _assert_(Vmean_input);
		Input *C_input             = this->GetInput(SmbCEnum);            _assert_(C_input);
		Input *Tz_input            = this->GetInput(SmbTzEnum);           _assert_(Tz_input);
		Input *Vz_input            = this->GetInput(SmbVzEnum);           _assert_(Vz_input);

		Tmean_input->GetInputValue(&Tmean,gauss);
		Vmean_input->GetInputValue(&Vmean,gauss);
		C_input->GetInputValue(&C,gauss);
		Tz_input->GetInputValue(&Tz,gauss);
		Vz_input->GetInputValue(&Vz,gauss);

	} else {
		this->GetInputValue(&Mappedpoint,SmbMappedforcingpointEnum);

		IssmDouble* tmean = NULL;
		IssmDouble* vmean = NULL;
		IssmDouble* cmean = NULL;
		IssmDouble* tzval = NULL;
		IssmDouble* vzval = NULL;
	
		parameters->FindParam(&tmean,&N,SmbTmeanParamEnum);
		parameters->FindParam(&vmean,&N,SmbVmeanParamEnum);
		parameters->FindParam(&cmean,&N,SmbCParamEnum);
		parameters->FindParam(&tzval,&N,SmbTzParamEnum);
		parameters->FindParam(&vzval,&N,SmbVzParamEnum);

		Tmean = tmean[Mappedpoint-1];
		Vmean = vmean[Mappedpoint-1];
		C = cmean[Mappedpoint-1];
		Tz = tzval[Mappedpoint-1];
		Vz = vzval[Mappedpoint-1];

		xDelete<IssmDouble>(tmean);
		xDelete<IssmDouble>(vmean);
		xDelete<IssmDouble>(cmean);
		xDelete<IssmDouble>(tzval);
		xDelete<IssmDouble>(vzval);
	}
	/*}}}*/

	/*First, check that the initial structures have been setup in GEMB. If not, initialize profile variables: layer thickness dz, * density d, temperature T, etc. {{{*/
	if(!isinitialized){
		if(VerboseSmb() && this->Sid()==0)_printf0_("smb core: Initializing grid\n");
		//if(this->Sid()==1) for(int i=0;i<m;i++)_printf_("z[" << i << "]=" <<
		//dz[i] << "\n");

		this->inputs->GetArray(SmbDziniEnum,this->lid,&dzini,&m);
		this->inputs->GetArray(SmbDiniEnum,this->lid,&dini,&m);
		this->inputs->GetArray(SmbReiniEnum,this->lid,&reini,&m);
		this->inputs->GetArray(SmbGdniniEnum,this->lid,&gdnini,&m);
		this->inputs->GetArray(SmbGspiniEnum,this->lid,&gspini,&m);
		this->inputs->GetArray(SmbWiniEnum,this->lid,&Wini,&m);
		this->inputs->GetArray(SmbAiniEnum,this->lid,&aini,&m);
		this->inputs->GetArray(SmbAdiffiniEnum,this->lid,&adiffini,&m);
		this->inputs->GetArray(SmbTiniEnum,this->lid,&Tini,&m);
		EC_input = this->GetInput(SmbECiniEnum);  _assert_(EC_input);
		EC_input->GetInputAverage(&EC);

		/*Retrieve the correct value of m (without the zeroes at the end)*/
		this->GetInputValue(&m,SmbSizeiniEnum);

		if(m==2){ //Snow properties are initialized with default values. Vertical grid has to be initialized too
			//            if(VerboseSmb() && this->Sid()==0)_printf0_("Snow properties initialized w DEFAULT values\n");

			/*initialize profile variables:*/
			GembgridInitialize(&dz, &m, zTop, dzTop, zMax, zY);

			d = xNew<IssmDouble>(m); for(int i=0;i<m;i++)d[i]=dini[0]; //ice density [kg m-3]
			re = xNew<IssmDouble>(m); for(int i=0;i<m;i++)re[i]=reini[0];         //set grain size to old snow [mm]
			gdn = xNew<IssmDouble>(m); for(int i=0;i<m;i++)gdn[i]=gdnini[0];         //set grain dentricity to old snow
			gsp = xNew<IssmDouble>(m); for(int i=0;i<m;i++)gsp[i]=gspini[0];         //set grain sphericity to old snow
			W = xNew<IssmDouble>(m); for(int i=0;i<m;i++)W[i]=Wini[0];             //set water content to zero [kg m-2]
			a = xNew<IssmDouble>(m); for(int i=0;i<m;i++)a[i]=aini[0];         //set albedo equal to fresh snow [fraction]
			adiff = xNew<IssmDouble>(m); for(int i=0;i<m;i++)adiff[i]=adiffini[0];         //set diffusive albedo equal to 1 [fraction]
			T = xNew<IssmDouble>(m); for(int i=0;i<m;i++)T[i]=Tmean;         //set initial grid cell temperature to the annual mean temperature [K]
			/*/!\ Default value of T can not be retrived from SMBgemb.m (like other snow properties)
			 *    because don't know Tmean yet when set default values.
			 *    Default value of 0C given in SMBgemb.m is overwritten here with value of Tmean*/

			//fixed lower temperature bounday condition - T is fixed
			T_bottom=T[m-1];
		}
		else{ //Retrieve snow properties from previous run. Need to provide values for all layers
			//            if(VerboseSmb() && this->Sid()==0)_printf0_("Snow properties initialized w RESTART values\n");

			dz = xNew<IssmDouble>(m);for(int i=0;i<m;i++)dz[i]=dzini[i];
			d = xNew<IssmDouble>(m);for(int i=0;i<m;i++)d[i]=dini[i];
			re = xNew<IssmDouble>(m);for(int i=0;i<m;i++)re[i]=reini[i];
			gdn = xNew<IssmDouble>(m);for(int i=0;i<m;i++)gdn[i]=gdnini[i];
			gsp = xNew<IssmDouble>(m);for(int i=0;i<m;i++)gsp[i]=gspini[i];
			W = xNew<IssmDouble>(m);for(int i=0;i<m;i++)W[i]=Wini[i];
			a = xNew<IssmDouble>(m);for(int i=0;i<m;i++)a[i]=aini[i];
			adiff = xNew<IssmDouble>(m);for(int i=0;i<m;i++)adiff[i]=adiffini[i];
			T = xNew<IssmDouble>(m);for(int i=0;i<m;i++)T[i]=Tini[i];

			//fixed lower temperature bounday condition - T is fixed
			_assert_(m>0);
			T_bottom=T[m-1];
		}

		this->SetElementInput(SmbAccumulatedECEnum,0.0);
		this->SetElementInput(SmbAccumulatedMassBalanceEnum,0.0);
		this->SetElementInput(SmbAccumulatedMeltEnum,0.0);
		this->SetElementInput(SmbAccumulatedRunoffEnum,0.0);
		this->SetElementInput(SmbAccumulatedRefreezeEnum,0.0);
		this->SetElementInput(SmbAccumulatedRainEnum,0.0);
		this->SetElementInput(SmbAccumulatedPrecipitationEnum,0.0);
		this->SetElementInput(SmbMSurfEnum,0.0);

		/*Flag the initialization:*/
		this->SetBoolInput(this->inputs,SmbIsInitializedEnum,true);
	}
	else{
		/*Recover inputs: */
		this->inputs->GetArray(SmbDzEnum,this->lid,&dz,&m);
		this->inputs->GetArray(SmbDEnum,this->lid,&d,&m);
		this->inputs->GetArray(SmbReEnum,this->lid,&re,&m);
		this->inputs->GetArray(SmbGdnEnum,this->lid,&gdn,&m);
		this->inputs->GetArray(SmbGspEnum,this->lid,&gsp,&m);
		this->inputs->GetArray(SmbWEnum,this->lid,&W,&m);
		this->inputs->GetArray(SmbAEnum,this->lid,&a,&m);
		this->inputs->GetArray(SmbAdiffEnum,this->lid,&adiff,&m);
		this->inputs->GetArray(SmbTEnum,this->lid,&T,&m);
		EC_input = this->GetInput(SmbECDtEnum);  _assert_(EC_input);
		EC_input->GetInputAverage(&EC);

		//fixed lower temperature bounday condition - T is fixed
		_assert_(m>0);
		T_bottom=T[m-1];
	} /*}}}*/

	// determine initial mass [kg]
	initMass=0; for(int i=0;i<m;i++) initMass += dz[i]*d[i] + W[i];

	// initialize cumulative variables
	sumR = 0; sumF=0; sumM = 0; sumEC = 0; sumP = 0; sumMassAdd = 0; sumMsurf = 0; sumRa = 0;

	//before starting loop, realize that the transient core runs this smb_core at time = time +deltaT.
   //go back to time - deltaT:
   time-=dt;

	if(VerboseSmb() && this->Sid()==0 && IssmComm::GetRank()==0)_printf0_("Time: t=" << setprecision(8) << timeinputs/365.0/24.0/3600.0 << " yr/" << (time+dt)/365.0/24.0/3600.0 << " yr" << setprecision(3) << " Step: " << count << "\n");

	/*Get daily accumulated inputs {{{*/
	if (count>1){
		Input *sumEC_input         = this->GetInput(SmbECEnum);  _assert_(sumEC_input);
		Input *sumM_input          = this->GetInput(SmbMeltEnum);  _assert_(sumM_input);
		Input *sumR_input          = this->GetInput(SmbRunoffEnum);  _assert_(sumR_input);
		Input *sumF_input          = this->GetInput(SmbRefreezeEnum);  _assert_(sumF_input);
		Input *sumRa_input         = this->GetInput(SmbRainEnum);  _assert_(sumRa_input);
		Input *sumP_input          = this->GetInput(SmbPrecipitationEnum);  _assert_(sumP_input);
		Input *ULW_input           = this->GetInput(SmbMeanULWEnum);  _assert_(ULW_input);
		Input *LW_input            = this->GetInput(SmbNetLWEnum);  _assert_(LW_input);
		Input *SW_input            = this->GetInput(SmbNetSWEnum);  _assert_(SW_input);
		Input *LHF_input           = this->GetInput(SmbMeanLHFEnum);  _assert_(LHF_input);
		Input *SHF_input           = this->GetInput(SmbMeanSHFEnum);  _assert_(SHF_input);
		Input *DzAdd_input         = this->GetInput(SmbDzAddEnum);  _assert_(DzAdd_input);
		Input *MassAdd_input       = this->GetInput(SmbMAddEnum);  _assert_(MassAdd_input);
		Input *InitMass_input      = this->GetInput(SmbMInitnum);  _assert_(InitMass_input);
		Input *sumMsurf_input      = this->GetInput(SmbMSurfSumEnum);  _assert_(sumMsurf_input);

		ULW_input->GetInputAverage(&meanULW);
		LW_input->GetInputAverage(&netLW);
		SW_input->GetInputAverage(&netSW);
		LHF_input->GetInputAverage(&meanLHF);
		SHF_input->GetInputAverage(&meanSHF);
		DzAdd_input->GetInputAverage(&dz_add);
		MassAdd_input->GetInputAverage(&sumMassAdd);
		sumMassAdd=sumMassAdd*dt;
		InitMass_input->GetInputAverage(&initMass);
		sumEC_input->GetInputAverage(&sumEC);
		sumEC=sumEC*dt*rho_ice;
		sumM_input->GetInputAverage(&sumM);
		sumM=sumM*dt*rho_ice;
		sumMsurf_input->GetInputAverage(&sumMsurf);
      sumMsurf=sumMsurf*dt*rho_ice;
		sumR_input->GetInputAverage(&sumR);
		sumR=sumR*dt*rho_ice;
		sumF_input->GetInputAverage(&sumF);
		sumF=sumF*dt*rho_ice;
		sumRa_input->GetInputAverage(&sumRa);
		sumRa=sumRa*dt*rho_ice;
		sumP_input->GetInputAverage(&sumP);
		sumP=sumP*dt*rho_ice;

		this->inputs->GetArray(SmbAccumulatedECSubstepEnum,this->lid,&ECsub,&steps);
		this->inputs->GetArray(SmbAccumulatedMeltSubstepEnum,this->lid,&Msub,&steps);
		this->inputs->GetArray(SmbAccumulatedRunoffSubstepEnum,this->lid,&Rsub,&steps);
		this->inputs->GetArray(SmbAccumulatedRefreezeSubstepEnum,this->lid,&Fsub,&steps);
		this->inputs->GetArray(SmbAccumulatedRainSubstepEnum,this->lid,&Rasub,&steps);
		this->inputs->GetArray(SmbAccumulatedPrecipitationSubstepEnum,this->lid,&Psub,&steps);
		this->inputs->GetArray(SmbAccumulatedMassBalanceSubstepEnum,this->lid,&SMBsub,&steps);
		this->inputs->GetArray(SmbFACSubstepEnum,this->lid,&FACsub,&steps);

		this->inputs->GetArray(SmbECSubstepEnum,this->lid,&ECdsub,&steps);
		this->inputs->GetArray(SmbMeltSubstepEnum,this->lid,&Mdsub,&steps);
		this->inputs->GetArray(SmbRunoffSubstepEnum,this->lid,&Rdsub,&steps);
		this->inputs->GetArray(SmbRefreezeSubstepEnum,this->lid,&Fdsub,&steps);
		this->inputs->GetArray(SmbRainSubstepEnum,this->lid,&Radsub,&steps);
		this->inputs->GetArray(SmbPrecipitationSubstepEnum,this->lid,&Pdsub,&steps);
		this->inputs->GetArray(SmbMassBalanceSubstepEnum,this->lid,&SMBdsub,&steps);
	}
	else{
		ECsub=xNewZeroInit<IssmDouble>(steps);
		SMBsub=xNewZeroInit<IssmDouble>(steps);
		Msub=xNewZeroInit<IssmDouble>(steps);
		Rsub=xNewZeroInit<IssmDouble>(steps);
		Fsub=xNewZeroInit<IssmDouble>(steps);
		Rasub=xNewZeroInit<IssmDouble>(steps);
		Psub=xNewZeroInit<IssmDouble>(steps);
		FACsub=xNewZeroInit<IssmDouble>(steps);

		this->inputs->SetArrayInput(SmbAccumulatedECSubstepEnum,this->lid,ECsub,steps);
		this->inputs->SetArrayInput(SmbAccumulatedMeltSubstepEnum,this->lid,Msub,steps);
		this->inputs->SetArrayInput(SmbAccumulatedRunoffSubstepEnum,this->lid,Rsub,steps);
		this->inputs->SetArrayInput(SmbAccumulatedRefreezeSubstepEnum,this->lid,Fsub,steps);
		this->inputs->SetArrayInput(SmbAccumulatedRainSubstepEnum,this->lid,Rasub,steps);
		this->inputs->SetArrayInput(SmbAccumulatedPrecipitationSubstepEnum,this->lid,Psub,steps);
		this->inputs->SetArrayInput(SmbAccumulatedMassBalanceSubstepEnum,this->lid,SMBsub,steps);
		this->inputs->SetArrayInput(SmbFACSubstepEnum,this->lid,FACsub,steps);

		ECdsub=xNewZeroInit<IssmDouble>(steps);
		SMBdsub=xNewZeroInit<IssmDouble>(steps);
		Mdsub=xNewZeroInit<IssmDouble>(steps);
		Rdsub=xNewZeroInit<IssmDouble>(steps);
		Fdsub=xNewZeroInit<IssmDouble>(steps);
		Radsub=xNewZeroInit<IssmDouble>(steps);
		Pdsub=xNewZeroInit<IssmDouble>(steps);

		this->inputs->SetArrayInput(SmbECSubstepEnum,this->lid,ECdsub,steps);
		this->inputs->SetArrayInput(SmbMeltSubstepEnum,this->lid,Mdsub,steps);
		this->inputs->SetArrayInput(SmbRunoffSubstepEnum,this->lid,Rdsub,steps);
		this->inputs->SetArrayInput(SmbRefreezeSubstepEnum,this->lid,Fdsub,steps);
		this->inputs->SetArrayInput(SmbRainSubstepEnum,this->lid,Radsub,steps);
		this->inputs->SetArrayInput(SmbPrecipitationSubstepEnum,this->lid,Pdsub,steps);
		this->inputs->SetArrayInput(SmbMassBalanceSubstepEnum,this->lid,SMBdsub,steps);
	}
	/*}}}*/

	// Get surface melt for albedo calculation
	Input *Msurf_input         = this->GetInput(SmbMSurfEnum);  _assert_(Msurf_input);
	Msurf_input->GetInputAverage(&Msurf);
	Msurf=Msurf*dt*rho_ice;

	Input *teValue_input= this->GetInput(SmbTeValueEnum,timeinputs); _assert_(teValue_input);
	Input *aValue_input= this->GetInput(SmbAValueEnum,timeinputs); _assert_(aValue_input);
	Input *dulwrfValue_input= this->GetInput(SmbDulwrfValueEnum,timeinputs); _assert_(dulwrfValue_input);
	Input *szaValue_input= this->GetInput(SmbSzaValueEnum,timeinputs); _assert_(szaValue_input);
	Input *cotValue_input= this->GetInput(SmbCotValueEnum,timeinputs); _assert_(cotValue_input);
	Input *ccsnowValue_input= this->GetInput(SmbCcsnowValueEnum,timeinputs); _assert_(ccsnowValue_input);
	Input *cciceValue_input= this->GetInput(SmbCciceValueEnum,timeinputs); _assert_(cciceValue_input);

	teValue_input->GetInputValue(&teValue,gauss);  // Emissivity [0-1]
	dulwrfValue_input->GetInputValue(&dulwrfValue,gauss);  // LWup perturbation [W m-2]
	aValue_input->GetInputValue(&aValue,gauss);  // Albedo [0 1]
	szaValue_input->GetInputValue(&szaValue,gauss);  // Solar Zenith Angle [degree]
	cotValue_input->GetInputValue(&cotValue,gauss);  // Cloud Optical Thickness
	ccsnowValue_input->GetInputValue(&ccsnowValue,gauss); //concentration of light absorbing carbon for snow [ppm1]
	cciceValue_input->GetInputValue(&cciceValue,gauss); //concentration of light absorbing carbon for ice [ppm1]

	/*extract daily data:{{{*/
	// Get time forcing inputs
	if (!ismappedforcing) {
		Input *Ta_input  = this->GetInput(SmbTaEnum,timeinputs);    _assert_(Ta_input);
		Input *V_input   = this->GetInput(SmbVEnum,timeinputs);     _assert_(V_input);
		Input *Dlwr_input= this->GetInput(SmbDlwrfEnum,timeinputs); _assert_(Dlwr_input);
		Input *Dswr_input= this->GetInput(SmbDswrfEnum,timeinputs); _assert_(Dswr_input);
		Input *Dswrdiff_input= this->GetInput(SmbDswdiffrfEnum,timeinputs); _assert_(Dswrdiff_input);
		Input *P_input   = this->GetInput(SmbPEnum,timeinputs);     _assert_(P_input);
		Input *eAir_input= this->GetInput(SmbEAirEnum,timeinputs);  _assert_(eAir_input);
		Input *pAir_input= this->GetInput(SmbPAirEnum,timeinputs);  _assert_(pAir_input);

		Ta_input->GetInputValue(&Ta,gauss);//screen level air temperature [K]
		V_input->GetInputValue(&V,gauss);  //wind speed [m s-1]
		Dlwr_input->GetInputValue(&dlw,gauss);   //downward longwave radiation flux [W m-2]
		Dswr_input->GetInputValue(&dsw,gauss);   //downward shortwave radiation flux [W m-2]
		Dswrdiff_input->GetInputValue(&dswdiff,gauss);   //downward shortwave diffuse radiation flux [W m-2]
		P_input->GetInputValue(&P,gauss);        //precipitation [kg m-2]
		eAir_input->GetInputValue(&eAir,gauss);  //screen level vapor pressure [Pa]
		pAir_input->GetInputValue(&pAir,gauss);  // screen level air pressure [Pa]
		//_printf_("Time: " << t << " Ta: " << Ta << " V: " << V << " dlw: " << dlw << " dsw: " << dsw << " P: " << P << " eAir: " << eAir << " pAir: " << pAir << "\n");
	} else {
		IssmDouble Dtol = 1e-11;
		IssmDouble gravity;
		parameters->FindParam(&gravity,ConstantsGEnum);

		int timestepping;
		IssmDouble dt;
		parameters->FindParam(&dt,TimesteppingTimeStepEnum);          /*transient core time step*/
		parameters->FindParam(&timestepping,TimesteppingTypeEnum);

		Input *currentsurface_input = this->GetInput(SurfaceEnum);  _assert_(currentsurface_input);
		currentsurface_input->GetInputAverage(&currentsurface);

		bool isprecipmap=true;
		parameters->FindParam(&isprecipmap,SmbIsprecipforcingremappedEnum);

		IssmDouble* tlapse = NULL;
		parameters->FindParam(&tlapse,&N,SmbLapseTaValueEnum); _assert_(tlapse);

		IssmDouble* dlwlapse = NULL;
		parameters->FindParam(&dlwlapse,&N,SmbLapsedlwrfValueEnum); _assert_(dlwlapse);

		IssmDouble* elevation = NULL;
		parameters->FindParam(&elevation,&N,SmbMappedforcingelevationEnum); _assert_(elevation);

		//Variables for downscaling
		IssmDouble taparam, dlwrfparam, rhparam, eaparam, pparam, prparam;
		parameters->FindParam(&taparam, Mappedpoint-1, timeinputs, timestepping, dt, SmbTaParamEnum);
		parameters->FindParam(&dlwrfparam, Mappedpoint-1, timeinputs, timestepping, dt, SmbDlwrfParamEnum);
		parameters->FindParam(&eaparam, Mappedpoint-1, timeinputs, timestepping, dt, SmbEAirParamEnum);
		parameters->FindParam(&pparam, Mappedpoint-1, timeinputs, timestepping, dt, SmbPAirParamEnum);
		parameters->FindParam(&prparam, Mappedpoint-1, timeinputs, timestepping, dt, SmbPParamEnum);

		//Variables not downscaled
		parameters->FindParam(&V, Mappedpoint-1, timeinputs, timestepping, dt, SmbVParamEnum);
		parameters->FindParam(&dsw, Mappedpoint-1, timeinputs, timestepping, dt, SmbDswrfParamEnum);
		parameters->FindParam(&dswdiff, Mappedpoint-1, timeinputs, timestepping, dt, SmbDswdiffrfParamEnum);

		Ta = taparam + (currentsurface - elevation[Mappedpoint-1])*tlapse[Mappedpoint-1];
		Tmean = Tmean + (currentsurface - elevation[Mappedpoint-1])*tlapse[Mappedpoint-1];
		if (fabs(dlwlapse[Mappedpoint-1]) > Dtol){
			dlw = fmax(dlwrfparam + (currentsurface - elevation[Mappedpoint-1])*dlwlapse[Mappedpoint-1],0.0);
		}else{
			//adjust downward longwave, holding emissivity equal (Glover et al, 1999)
			IssmDouble SB = 5.67e-8; // Stefan-Boltzmann constant (W m-2 K-4)
			IssmDouble effe = 1.;
			effe = dlwrfparam/(SB * pow(taparam,4.0));
			dlw = fmax(effe*SB*pow(Ta,4.0),0.0);
		}

		if ( (fabs(dlwlapse[Mappedpoint-1]) > Dtol) || (fabs(tlapse[Mappedpoint-1]) > Dtol)){
			IssmDouble Rg = 8.314; // gas constant (J mol-1 K-1)
			IssmDouble dAir = 0.0;
			// calculated air density [kg/m3]
			//    dAir = 0.029 * pAir /(R * Ta);
			dAir=0.029 * pparam /(Rg * Ta);
			pAir=pparam-gravity*dAir*(currentsurface - elevation[Mappedpoint-1]);
		}
		else pAir=pparam;

		//Hold relative humidity constant, calculte new saturation vapor pressure,
		// and new saturation specific humidity to scale precipitation
		//https://cran.r-project.org/web/packages/humidity/vignettes/humidity-measures.html
		//es over ice calculation
		//Ding et al., 2019 after Bolton, 1980
		//ea37 = rh37*100*6.112.*exp((17.67*(t237-273.15))./(t237-29.65));
		//ea37s (hPa) = 6.1121*exp(22.587*(t-273.15)/(t-273.15+273.86)); (with respect to ice)
		IssmDouble esparam, es, qsparam, qs;
		esparam=6.112*exp((17.67*(taparam-273.15))/(taparam-29.65));
		es=6.112*exp((17.67*(Ta-273.15))/(Ta-29.65));
		rhparam=eaparam/esparam;
		eAir=fmax(rhparam*es,0.0);
		qs=fmax(0.622*es/(pAir/100 - 0.378*es),0);
		qsparam=fmax(0.622*esparam/(pparam/100 - 0.378*esparam),0);

		if ((isprecipmap) && (qsparam>0)){ 
			P=fmax(prparam*qs/qsparam,0.0);
			C=fmax(C*qs/qsparam,0.0);
		}
		else P=prparam;

		xDelete<IssmDouble>(elevation);
		xDelete<IssmDouble>(tlapse);
		xDelete<IssmDouble>(dlwlapse);
	}
	/*}}}*/

	/*Snow grain metamorphism:*/
	if(isgraingrowth)grainGrowth(&re, &gdn, &gsp, T, dz, d, W, smb_dt, m, aIdx,this->Sid());

	/*Snow, firn and ice albedo:*/
	if(isalbedo)albedo(&a,&adiff,aIdx,re,dz,d,cldFrac,aIce,aSnow,aValue,adThresh,T,W,P,EC,Msurf,ccsnowValue,cciceValue,szaValue,cotValue,t0wet,t0dry,K,smb_dt,rho_ice,m,this->Sid());

	/*Distribution of absorbed short wave radation with depth:*/
	if(isshortwave)shortwave(&swf, swIdx, aIdx, dsw, dswdiff, a[0], adiff[0], d, dz, re,rho_ice,m,this->Sid());

	/*Calculate net shortwave [W m-2]*/
	netSW = netSW + cellsum(swf,m)*smb_dt/dt;

	if(isconstrainsurfaceT){
		if (m>0) T[0]=Ta;
		if (m>1) T[1]=Ta;
	}
	/*Thermal profile computation:*/
	if(isthermal)thermo(&shf, &lhf, &EC, &T, &ulw, re, dz, d, swf, dlw, Ta, V, eAir, pAir, tcIdx, eIdx, teValue, dulwrfValue, teThresh, W[0], smb_dt, dzMin, m, Vz, Tz, thermo_scaling,rho_ice,this->Sid(),isconstrainsurfaceT,isdeltaLWup);

	/*Change in thickness of top cell due to evaporation/condensation  assuming same density as top cell.
	 * need to fix this in case all or more of cell evaporates */
	dz[0] = dz[0] + EC / d[0];

	/*Add snow/rain to top grid cell adjusting cell depth, temperature and density*/
	if(isaccumulation)accumulation(&T, &dz, &d, &W, &a, &adiff, &re, &gdn, &gsp, &Ra, &m, aIdx, dsnowIdx, Tmean, Ta, P, dzMin, aSnow, C, V, Vmean, rho_ice,this->Sid());

	/*Calculate water production, M [kg m-2] resulting from snow/ice temperature exceeding 273.15 deg K
	 * (> 0 deg C), runoff R [kg m-2] and resulting changes in density and determine wet compaction [m]*/
	if(ismelt)melt(&M, &Msurf, &R, &F, &mAdd, &dz_add, &T, &d, &dz, &W, &a, &adiff, &re, &gdn, &gsp, &m, Ra, dzMin, zMax, zMin, zTop, zY, rho_ice,this->Sid());

	/*Allow non-melt densification and determine compaction [m]*/
	if(isdensification)densification(&d,&dz, T, re, denIdx, aIdx, swIdx, adThresh, C, smb_dt, Tmean,rho_ice,m,this->Sid());

	/*Calculate net longwave [W m-2]*/
	meanULW = meanULW + ulw*smb_dt/dt;
	netLW = netLW + (dlw - ulw)*smb_dt/dt;

	/*Verbose some results in debug mode: {{{*/
	if(VerboseSmb() && 0){
		_printf_("smb log: count[" << count << "] m[" << m << "] "
					<< setprecision(16)   << "T[" << cellsum(T,m)  << "] "
					<< "d[" << cellsum(d,m)  << "] "
					<< "dz[" << cellsum(dz,m)  << "] "
					<< "a[" << cellsum(a,m)  << "] "
					<< "W[" << cellsum(W,m)  << "] "
					<< "re[" << cellsum(re,m)  << "] "
					<< "gdn[" << cellsum(gdn,m)  << "] "
					<< "gsp[" << cellsum(gsp,m)  << "] "
					<< "swf[" << netSW << "] "
					<< "lwf[" << netLW << "] "
					<< "a[" << a << "] "
					<< "adiff[" << adiff << "] "
					<< "te[" << teValue << "] "
					<< "\n");
	} /*}}}*/

	meanLHF = meanLHF + lhf*smb_dt/dt;
	meanSHF = meanSHF + shf*smb_dt/dt;

	/*Sum component mass changes [kg m-2]*/
	sumMassAdd = mAdd + sumMassAdd;
	sumM = M + sumM;
	sumMsurf = Msurf + sumMsurf;
	sumR = R + sumR;
	sumW = cellsum(W,m);
	sumP = P +  sumP;
	sumEC = EC + sumEC;  // evap (-)/cond(+)
	sumRa = Ra + sumRa;
	sumF = F + sumF;

	/*Calculate total system mass:*/
	sumMass=0;
	fac=0;
	T0m=0;
	T10m=0;
	T30m=0;
	T50m=0;
	for(int i=0;i<m;i++){
		sumMass += dz[i]*d[i];
		sumH += dz[i];
		if (d[i] > 0) fac += dz[i]*(rho_ice - fmin(d[i],rho_ice));
		if (i==0 || (d[i]<rho_ice && d[i]>0 && sumH <= 50)) T50m = T[i];
		if (i==0 || (d[i]<rho_ice && d[i]>0 && sumH <= 30)) T30m = T[i];
		if (i==0 || (d[i]<rho_ice && d[i]>0 && sumH <= 10)) T10m = T[i];
		if (i==0) T0m = T[i];
	}

	#if defined(_HAVE_AD_)
	/*we want to avoid the round operation at all cost. Not differentiable.*/
	_error_("not implemented yet");
	#else
	dMass = sumMass + sumR + sumW - sumP - sumEC - initMass - sumMassAdd;
	dMass = round(dMass * 100.0)/100.0;

	/*Check mass conservation:*/
	if (dMass != 0.0){
		_printf_("Time: " << setprecision(8) << timeinputs/365.0/24.0/3600.0 << ",total system mass not conserved in MB function \n");
		_printf_("sumMass: " << sumMass << " sumR: " << sumR << " sumW: " << sumW << " sumP: " << sumP << " sumEC: " << sumEC << " initMass: " << initMass << " sumMassAdd: " << sumMassAdd << " \n");
	}
	#endif

	/*Check bottom grid cell T is unchanged:*/
	if(VerboseSmb() && this->Sid()==0 && IssmComm::GetRank()==0){
		if (T[m-1]!=T_bottom) _printf_("T(end)~=T_bottom" << "\n");
	}

	/*Save accumulated inputs {{{*/
   Input *accsumEC_input         = this->GetInput(SmbAccumulatedECEnum);  _assert_(accsumEC_input);
   Input *accsumM_input          = this->GetInput(SmbAccumulatedMeltEnum);  _assert_(accsumM_input);
	Input *accsumR_input          = this->GetInput(SmbAccumulatedRunoffEnum);  _assert_(accsumR_input);
	Input *accsumF_input          = this->GetInput(SmbAccumulatedRefreezeEnum);  _assert_(accsumF_input);
	Input *accsumRa_input         = this->GetInput(SmbAccumulatedRainEnum);  _assert_(accsumRa_input);
	Input *accsumP_input          = this->GetInput(SmbAccumulatedPrecipitationEnum);  _assert_(accsumP_input);
	Input *accsumSMB_input        = this->GetInput(SmbAccumulatedMassBalanceEnum);  _assert_(accsumSMB_input);

   accsumEC_input->GetInputAverage(&accsumEC);
	accsumM_input->GetInputAverage(&accsumM);
	accsumR_input->GetInputAverage(&accsumR);
	accsumF_input->GetInputAverage(&accsumF);
	accsumRa_input->GetInputAverage(&accsumRa);
	accsumP_input->GetInputAverage(&accsumP);
	accsumSMB_input->GetInputAverage(&accsumSMB);

	this->SetElementInput(SmbAccumulatedECEnum,accsumEC+EC/rho_ice);
	this->SetElementInput(SmbAccumulatedMassBalanceEnum,accsumSMB+(P + EC -R)/rho_ice);
	this->SetElementInput(SmbAccumulatedMeltEnum,accsumM+M/rho_ice);
	this->SetElementInput(SmbAccumulatedRunoffEnum,accsumR+R/rho_ice);
	this->SetElementInput(SmbAccumulatedRefreezeEnum,accsumF+F/rho_ice);
	this->SetElementInput(SmbAccumulatedRainEnum,accsumRa+Ra/rho_ice);
	this->SetElementInput(SmbAccumulatedPrecipitationEnum,accsumP+P/rho_ice);

	ECsub[count-1]=accsumEC+EC/rho_ice;
	Msub[count-1]=accsumM+M/rho_ice;
	Fsub[count-1]=accsumF+F/rho_ice;
	Rsub[count-1]=accsumR+R/rho_ice;
	Rasub[count-1]=accsumRa+Ra/rho_ice;
	Psub[count-1]=accsumP+P/rho_ice;
	SMBsub[count-1]=accsumSMB+(P + EC -R)/rho_ice;
	FACsub[count-1]=fac/1000.;

	ECdsub[count-1]=EC/rho_ice;
	Mdsub[count-1]=M/rho_ice;
	Fdsub[count-1]=F/rho_ice;
	Rdsub[count-1]=R/rho_ice;
	Radsub[count-1]=Ra/rho_ice;
	Pdsub[count-1]=P/rho_ice;
	SMBdsub[count-1]=(P + EC -R)/rho_ice;

	this->inputs->SetArrayInput(SmbAccumulatedECSubstepEnum,this->lid,ECsub,steps);
	this->inputs->SetArrayInput(SmbAccumulatedMeltSubstepEnum,this->lid,Msub,steps);
	this->inputs->SetArrayInput(SmbAccumulatedRunoffSubstepEnum,this->lid,Rsub,steps);
	this->inputs->SetArrayInput(SmbAccumulatedRefreezeSubstepEnum,this->lid,Fsub,steps);
	this->inputs->SetArrayInput(SmbAccumulatedRainSubstepEnum,this->lid,Rasub,steps);
	this->inputs->SetArrayInput(SmbAccumulatedPrecipitationSubstepEnum,this->lid,Psub,steps);
	this->inputs->SetArrayInput(SmbAccumulatedMassBalanceSubstepEnum,this->lid,SMBsub,steps);
	this->inputs->SetArrayInput(SmbFACSubstepEnum,this->lid,FACsub,steps);

	this->inputs->SetArrayInput(SmbECSubstepEnum,this->lid,ECdsub,steps);
	this->inputs->SetArrayInput(SmbMeltSubstepEnum,this->lid,Mdsub,steps);
	this->inputs->SetArrayInput(SmbRunoffSubstepEnum,this->lid,Rdsub,steps);
	this->inputs->SetArrayInput(SmbRefreezeSubstepEnum,this->lid,Fdsub,steps);
	this->inputs->SetArrayInput(SmbRainSubstepEnum,this->lid,Radsub,steps);
	this->inputs->SetArrayInput(SmbPrecipitationSubstepEnum,this->lid,Pdsub,steps);
	this->inputs->SetArrayInput(SmbMassBalanceSubstepEnum,this->lid,SMBdsub,steps);
   /*}}}*/

	/*Save generated inputs: */
	this->inputs->SetArrayInput(SmbDzEnum,this->lid,dz,m);
	this->inputs->SetArrayInput(SmbDEnum,this->lid,d,m);
	this->inputs->SetArrayInput(SmbReEnum,this->lid,re,m);
	this->inputs->SetArrayInput(SmbGdnEnum,this->lid,gdn,m);
	this->inputs->SetArrayInput(SmbGspEnum,this->lid,gsp,m);
	this->inputs->SetArrayInput(SmbTEnum,this->lid,T,m);
	this->inputs->SetArrayInput(SmbWEnum,this->lid,W,m);
	this->inputs->SetArrayInput(SmbAEnum,this->lid,a,m);
	this->inputs->SetArrayInput(SmbAdiffEnum,this->lid,adiff,m);
	this->SetElementInput(SmbECEnum,cellsum(ECdsub,steps)/dt);
	this->SetElementInput(SmbMassBalanceEnum,cellsum(SMBdsub,steps)/dt);
	this->SetElementInput(SmbMeltEnum,cellsum(Mdsub,steps)/dt);
	this->SetElementInput(SmbRunoffEnum,cellsum(Rdsub,steps)/dt);
	this->SetElementInput(SmbRefreezeEnum,cellsum(Fdsub,steps)/dt);
	this->SetElementInput(SmbRainEnum,cellsum(Radsub,steps)/dt);
	this->SetElementInput(SmbPrecipitationEnum,cellsum(Pdsub,steps)/dt);
	this->SetElementInput(SmbMeanULWEnum,meanULW);
	this->SetElementInput(SmbNetLWEnum,netLW);
	this->SetElementInput(SmbNetSWEnum,netSW);
	this->SetElementInput(SmbMeanLHFEnum,meanLHF);
	this->SetElementInput(SmbMeanSHFEnum,meanSHF);
	this->SetElementInput(SmbDzAddEnum,dz_add);
	this->SetElementInput(SmbMInitnum,initMass);
	this->SetElementInput(SmbMSurfEnum,Msurf/dt/rho_ice);
	this->SetElementInput(SmbMAddEnum,sumMassAdd/dt);
	this->SetElementInput(SmbMSurfSumEnum,sumMsurf/dt/rho_ice);
	this->SetElementInput(SmbWAddEnum,sumW/dt);
	this->SetElementInput(SmbFACEnum,fac/1000.); // output in meters
	this->SetElementInput(SmbTsEnum,T0m); // output in K at surface
	this->SetElementInput(SmbT10Enum,T10m); // output in K at 10m depth
	this->SetElementInput(SmbT30Enum,T30m); // output in K at 10m depth
	this->SetElementInput(SmbT50Enum,T50m); // output in K ar 10m depth
	this->SetElementInput(SmbECDtEnum,EC);

	/*Free allocations:{{{*/
	if(dz) xDelete<IssmDouble>(dz);
	if(d) xDelete<IssmDouble>(d);
	if(re) xDelete<IssmDouble>(re);
	if(gdn) xDelete<IssmDouble>(gdn);
	if(gsp) xDelete<IssmDouble>(gsp);
	if(W) xDelete<IssmDouble>(W);
	if(a) xDelete<IssmDouble>(a);
	if(adiff) xDelete<IssmDouble>(adiff);
	if(T) xDelete<IssmDouble>(T);
	if(dzini) xDelete<IssmDouble>(dzini);
	if(dini) xDelete<IssmDouble>(dini);
	if(reini) xDelete<IssmDouble>(reini);
	if(gdnini) xDelete<IssmDouble>(gdnini);
	if(gspini) xDelete<IssmDouble>(gspini);
	if(Wini) xDelete<IssmDouble>(Wini);
	if(aini) xDelete<IssmDouble>(aini);
	if(adiffini) xDelete<IssmDouble>(adiffini);
	if(Tini) xDelete<IssmDouble>(Tini);
	if(swf) xDelete<IssmDouble>(swf);

	if(ECsub) xDelete<IssmDouble>(ECsub);
	if(SMBsub) xDelete<IssmDouble>(SMBsub);
	if(Msub) xDelete<IssmDouble>(Msub);
	if(Rsub) xDelete<IssmDouble>(Rsub);
	if(Fsub) xDelete<IssmDouble>(Fsub);
	if(Rasub) xDelete<IssmDouble>(Rasub);
	if(Psub) xDelete<IssmDouble>(Psub);
	if(FACsub) xDelete<IssmDouble>(FACsub);
	if(ECdsub) xDelete<IssmDouble>(ECdsub);
	if(SMBdsub) xDelete<IssmDouble>(SMBdsub);
	if(Mdsub) xDelete<IssmDouble>(Mdsub);
	if(Rdsub) xDelete<IssmDouble>(Rdsub);
	if(Fdsub) xDelete<IssmDouble>(Fdsub);
	if(Radsub) xDelete<IssmDouble>(Radsub);
	if(Pdsub) xDelete<IssmDouble>(Pdsub);

	delete gauss;
	/*}}}*/
}
/*}}}*/
void       Element::SubglacialWaterPressure(int output_enum){/*{{{*/

	bool        ispwHydroArma;
	int         numvertices = this->GetNumberOfVertices();
	IssmDouble  p_water[MAXVERTICES];
	IssmDouble *perturbationvalues   = xNew<IssmDouble>(numvertices);
	Gauss      *gauss     = this->NewGauss();
	Friction   *friction  = new Friction(this);
   /*Calculate subglacial water pressure*/
   for(int i=0;i<numvertices;i++){
         gauss->GaussVertex(i);
         p_water[i] = friction->SubglacialWaterPressure(gauss);
   }
   /*Add perturbation in water pressure if HydrologyIsWaterPressureArmaEnum is true*/
   this->parameters->FindParam(&ispwHydroArma,HydrologyIsWaterPressureArmaEnum);
   if(ispwHydroArma){
      this->GetInputListOnVertices(perturbationvalues,WaterPressureArmaPerturbationEnum);
		for(int i=0;i<numvertices;i++) p_water[i] = p_water[i]+perturbationvalues[i];
   }
   /*Save*/
   this->AddInput(output_enum,p_water,P1DGEnum);
   /*Clean-up*/
   delete gauss;
   delete friction;
   xDelete<IssmDouble>(perturbationvalues);

}/*}}}*/
void       Element::StrainRateESA(IssmDouble* epsilon,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input){/*{{{*/

	/*Intermediaries*/
	IssmDouble dvx[3];
	IssmDouble dvy[3];

	/*Check that both inputs have been found*/
	if(!vx_input || !vy_input){
		_error_("Input missing. Here are the input pointers we have for vx: " << vx_input << ", vy: " << vy_input << "\n");
	}

	/*Get strain rate assuming that epsilon has been allocated*/
	vx_input->GetInputDerivativeValue(&dvx[0],xyz_list,gauss);
	vy_input->GetInputDerivativeValue(&dvy[0],xyz_list,gauss);
	epsilon[0] = dvx[0];	// normal strain rate x-direction
	epsilon[1] = dvy[1]; // normal strain rate y-direction
	epsilon[2] = 0.5*(dvx[1] + dvy[0]); // shear strain rate
	epsilon[3] = 0.5*(dvx[1] - dvy[0]); // rotation rate

}/*}}}*/
void       Element::StrainRateFS(IssmDouble* epsilon,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input,Input* vz_input){/*{{{*/
	/*Compute the 3d Strain Rate (6 components):
	 *
	 * epsilon=[exx eyy ezz exy exz eyz]
	 */

	/*Intermediaries*/
	IssmDouble dvx[3];
	IssmDouble dvy[3];
	IssmDouble dvz[3];

	/*Check that both inputs have been found*/
	if (!vx_input || !vy_input || !vz_input){
		_error_("Input missing. Here are the input pointers we have for vx: " << vx_input << ", vy: " << vy_input << ", vz: " << vz_input << "\n");
	}

	/*Get strain rate assuming that epsilon has been allocated*/
	vx_input->GetInputDerivativeValue(&dvx[0],xyz_list,gauss);
	vy_input->GetInputDerivativeValue(&dvy[0],xyz_list,gauss);
	vz_input->GetInputDerivativeValue(&dvz[0],xyz_list,gauss);

	epsilon[0] = dvx[0];
	epsilon[1] = dvy[1];
	epsilon[2] = dvz[2];
	epsilon[3] = 0.5*(dvx[1] + dvy[0]);
	epsilon[4] = 0.5*(dvx[2] + dvz[0]);
	epsilon[5] = 0.5*(dvy[2] + dvz[1]);

}/*}}}*/
void       Element::StrainRateHO(IssmDouble* epsilon,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input){/*{{{*/
	/*Compute the 3d Blatter/HOStrain Rate (5 components):
	 *
	 * epsilon=[exx eyy exy exz eyz]
	 *
	 * with exz=1/2 du/dz
	 *      eyz=1/2 dv/dz
	 *
	 * the contribution of vz is neglected
	 */

	/*Intermediaries*/
	IssmDouble dvx[3];
	IssmDouble dvy[3];

	/*Check that both inputs have been found*/
	if (!vx_input || !vy_input){
		_error_("Input missing. Here are the input pointers we have for vx: " << vx_input << ", vy: " << vy_input << "\n");
	}

	/*Get strain rate assuming that epsilon has been allocated*/
	vx_input->GetInputDerivativeValue(&dvx[0],xyz_list,gauss);
	vy_input->GetInputDerivativeValue(&dvy[0],xyz_list,gauss);
	epsilon[0] = dvx[0];
	epsilon[1] = dvy[1];
	epsilon[2] = 0.5*(dvx[1] + dvy[0]);
	epsilon[3] = 0.5*dvx[2];
	epsilon[4] = 0.5*dvy[2];

}/*}}}*/
void       Element::StrainRateHO2dvertical(IssmDouble* epsilon,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input){/*{{{*/
	/*Compute the 2d Blatter/HOStrain Rate (2 components):
	 *
	 * epsilon=[exx exz]
	 *
	 * with exz=1/2 du/dz
	 *
	 * the contribution of vz is neglected
	 */

	/*Intermediaries*/
	IssmDouble dvx[3];

	/*Check that both inputs have been found*/
	if (!vx_input){
		_error_("Input missing. Here are the input pointers we have for vx: " << vx_input <<"\n");
	}

	/*Get strain rate assuming that epsilon has been allocated*/
	vx_input->GetInputDerivativeValue(&dvx[0],xyz_list,gauss);
	epsilon[0] = dvx[0];
	epsilon[1] = 0.5*dvx[1];

}/*}}}*/
void       Element::StrainRateMOLHO(IssmDouble* epsilon,IssmDouble* xyz_list,Gauss* gauss,Input* vxbase_input,Input* vybase_input,Input* vxshear_input,Input* vyshear_input,Input* thickness_input,Input* n_input,IssmDouble zeta){/*{{{*/
	/*Compute the 2d Blatter/MOLHO Strain Rate (5 components) for a given vertical coordinate (zeta):
	 *
	 * epsilon=[exx eyy exy exz eyz]
	 *
	 * with exz=1/2 du/dz
	 *      eyz=1/2 dv/dz
	 *
	 * the contribution of vz is neglected
	 *
	 * zeta = (surface - z)/thickness 
	 * i.e., zeta=0 <=> ice surface
	 * i.e., zeta=1 <=> ice base
	 */

	/*Intermediaries*/
	IssmDouble dvxb[2];  // [dvx/dx, dvx/dy] at the base
	IssmDouble dvyb[2];  // [dvy/dx, dvy/dy] at the base
	IssmDouble dvxsh[2]; // [dvx/dx, dvx/dy] at the surface (only the shear)
	IssmDouble dvysh[2]; // [dvy/dx, dvy/dy] at the surface (only the shear)
	IssmDouble vxsh,vysh; // shear vx and shear at the surface
	IssmDouble thickness,n;

	/*Check that both inputs have been found*/
	if (!vxbase_input || !vybase_input || !vxshear_input || !vyshear_input || !thickness_input || !n_input){
		_error_("Input missing. Here are the input pointers we have for vxb: " << vxbase_input << ", vyb: " << vybase_input << ", vxshear: " << vxshear_input << ", vyshear: " << vyshear_input << ", thickness: " << thickness_input << ", n: " << n_input<< "\n");
	}

	/*Get strain rate assuming that epsilon has been allocated*/
	vxbase_input->GetInputDerivativeValue(&dvxb[0],xyz_list,gauss);
	vybase_input->GetInputDerivativeValue(&dvyb[0],xyz_list,gauss);
	vxshear_input->GetInputDerivativeValue(&dvxsh[0],xyz_list,gauss);
	vyshear_input->GetInputDerivativeValue(&dvysh[0],xyz_list,gauss);
	thickness_input->GetInputValue(&thickness,gauss); _assert_(thickness>0);
   vxshear_input->GetInputValue(&vxsh,gauss);
   vyshear_input->GetInputValue(&vysh,gauss);
   n_input->GetInputValue(&n,gauss);

	epsilon[0] = dvxb[0] + dvxsh[0]*( 1 - pow(zeta,n+1) );
	epsilon[1] = dvyb[1] + dvysh[1]*( 1 - pow(zeta,n+1) );
	epsilon[2] = 0.5*( dvxb[1] + dvxsh[1]*( 1 - pow(zeta,n+1) ) + dvyb[0] + dvysh[0]*( 1 - pow(zeta,n+1) ) ); 
	epsilon[3] = 0.5*(vxsh/thickness)*(n+1)*pow(zeta,n);
	epsilon[4] = 0.5*(vysh/thickness)*(n+1)*pow(zeta,n);

}/*}}}*/
void       Element::StrainRateSSA(IssmDouble* epsilon,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input){/*{{{*/

	/*Intermediaries*/
	IssmDouble dvx[3];
	IssmDouble dvy[3];

	/*Check that both inputs have been found*/
	if(!vx_input || !vy_input){
		_error_("Input missing. Here are the input pointers we have for vx: " << vx_input << ", vy: " << vy_input << "\n");
	}

	/*Get strain rate assuming that epsilon has been allocated*/
	vx_input->GetInputDerivativeValue(&dvx[0],xyz_list,gauss);
	vy_input->GetInputDerivativeValue(&dvy[0],xyz_list,gauss);
	epsilon[0] = dvx[0];
	epsilon[1] = dvy[1];
	epsilon[2] = 0.5*(dvx[1] + dvy[0]);

}/*}}}*/
void       Element::StrainRateSSA1d(IssmDouble* epsilon,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input){/*{{{*/

	/*Intermediaries*/
	IssmDouble dvx[3];

	/*Check that both inputs have been found*/
	if (!vx_input){
		_error_("Input missing. Here are the input pointers we have for vx: " << vx_input << "\n");
	}

	/*Get strain rate assuming that epsilon has been allocated*/
	vx_input->GetInputDerivativeValue(&dvx[0],xyz_list,gauss);
	*epsilon = dvx[0];

}/*}}}*/
void       Element::StressMaxPrincipalCreateInput(void){/*{{{*/

	/*Intermediaries*/
	IssmDouble *xyz_list = NULL;
	IssmDouble  sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz;
	IssmDouble  a,b,c,d,x[3],max;
	int         dim,numroots;

	/*First: get stress tensor*/
	this->ComputeStressTensor();

	/*Get domain dimension*/
	this->FindParam(&dim,DomainDimensionEnum);

	/*Fetch number vertices and allocate memory*/
	const int NUM_VERTICES  = this->GetNumberOfVertices();
	IssmDouble maxprincipal[MAXVERTICES];

	/*Retrieve all inputs and parameters*/
	this->GetVerticesCoordinatesBase(&xyz_list);
	Input* sigma_xx_input  = this->GetInput(StressTensorxxEnum); _assert_(sigma_xx_input);
	Input* sigma_yy_input  = this->GetInput(StressTensoryyEnum); _assert_(sigma_yy_input);
	Input* sigma_xy_input  = this->GetInput(StressTensorxyEnum); _assert_(sigma_xy_input);
	Input* sigma_xz_input  = NULL;
	Input* sigma_yz_input  = NULL;
	Input* sigma_zz_input  = NULL;
	if(dim==3){
		sigma_xz_input  = this->GetInput(StressTensorxzEnum); _assert_(sigma_xz_input);
		sigma_yz_input  = this->GetInput(StressTensoryzEnum); _assert_(sigma_yz_input);
		sigma_zz_input  = this->GetInput(StressTensorzzEnum); _assert_(sigma_zz_input);
	}

	/*loop over vertices: */
	Gauss* gauss=this->NewGauss();
	for(int iv=0;iv<NUM_VERTICES;iv++){
		gauss->GaussVertex(iv);

		sigma_xx_input->GetInputValue(&sigma_xx,gauss);
		sigma_yy_input->GetInputValue(&sigma_yy,gauss);
		sigma_xy_input->GetInputValue(&sigma_xy,gauss);
		if(dim==3){
			sigma_xz_input->GetInputValue(&sigma_xz,gauss);
			sigma_yz_input->GetInputValue(&sigma_yz,gauss);
			sigma_zz_input->GetInputValue(&sigma_zz,gauss);
		}

		if(dim==2){
			a = 0.;
			b = 1.;
			c = -sigma_yy -sigma_xx;
			d = sigma_xx*sigma_yy - sigma_xy*sigma_xy;
		}
		else{
			a = -1.;
			b = sigma_xx+sigma_yy+sigma_zz;
			c = -sigma_xx*sigma_yy -sigma_xx*sigma_zz -sigma_yy*sigma_zz + sigma_xy*sigma_xy +sigma_xz*sigma_xz +sigma_yz*sigma_yz;
			d = sigma_xx*sigma_yy*sigma_zz - sigma_xx*sigma_yz*sigma_yz -sigma_yy*sigma_xz*sigma_xz - sigma_zz*sigma_xy*sigma_xy + 2.*sigma_xy*sigma_xz*sigma_yz;
		}

		/*Get roots of polynomials*/
		cubic(a,b,c,d,x,&numroots);

		/*Initialize maximum eigne value*/
		if(numroots>0){
			max = fabs(x[0]);
		}
		else{
			_error_("No eigen value found");
		}

		/*Get max*/
		for(int i=1;i<numroots;i++){
			if(fabs(x[i])>max) max = fabs(x[i]);
		}

		maxprincipal[iv]=max;
	}

	/*Create input*/
	this->AddInput(StressMaxPrincipalEnum,&maxprincipal[0],P1Enum);

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	delete gauss;
}
/*}}}*/
IssmDouble Element::TotalFloatingBmb(IssmDouble* mask, bool scaled){/*{{{*/

	/*Retrieve values of the mask defining the element: */
	for(int i=0;i<this->GetNumberOfVertices();i++){
		if(mask[this->vertices[i]->Sid()]<=0.){
			return 0.;
		}
	}

	/*Return: */
	return this->TotalFloatingBmb(scaled);
}
/*}}}*/
IssmDouble Element::TotalGroundedBmb(IssmDouble* mask, bool scaled){/*{{{*/

	/*Retrieve values of the mask defining the element: */
	for(int i=0;i<this->GetNumberOfVertices();i++){
		if(mask[this->vertices[i]->Sid()]<=0.){
			return 0.;
		}
	}

	/*Return: */
	return this->TotalGroundedBmb(scaled);
}
/*}}}*/
IssmDouble Element::TotalSmb(IssmDouble* mask, bool scaled){/*{{{*/

	/*Retrieve values of the mask defining the element: */
	for(int i=0;i<this->GetNumberOfVertices();i++){
		if(mask[this->vertices[i]->Sid()]<=0.){
			return 0.;
		}
	}

	/*Return: */
	return this->TotalSmb(scaled);
}
/*}}}*/
IssmDouble Element::TotalSmbMelt(IssmDouble* mask, bool scaled){/*{{{*/

	/*Retrieve values of the mask defining the element: */
	for(int i=0;i<this->GetNumberOfVertices();i++){
		if(mask[this->vertices[i]->Sid()]<=0.){
			return 0.;
		}
	}

	/*Return: */
	return this->TotalSmbMelt(scaled);
}
/*}}}*/
IssmDouble Element::TotalSmbRefreeze(IssmDouble* mask, bool scaled){/*{{{*/

	/*Retrieve values of the mask defining the element: */
	for(int i=0;i<this->GetNumberOfVertices();i++){
		if(mask[this->vertices[i]->Sid()]<=0.){
			return 0.;
		}
	}

	/*Return: */
	return this->TotalSmbRefreeze(scaled);
}
/*}}}*/
void       Element::TransformInvStiffnessMatrixCoord(ElementMatrix* Ke,int transformenum){/*{{{*/

	/*All nodes have the same Coordinate System*/
	int numnodes  = this->GetNumberOfNodes();

	/*Do we need to actually do anything? (only if we have a rotated coordinate system)*/
	bool isrotation = false;
	for(int i=0;i<numnodes;i++){
		if(this->nodes[i]->isrotated){
			isrotation=true;
			break;
		}
	}
	if(!isrotation) return;

	/*Rotate stiffness*/
	int* cs_array = xNew<int>(numnodes);
	for(int i=0;i<numnodes;i++) cs_array[i]=transformenum;
	TransformInvStiffnessMatrixCoord(Ke,this->nodes,numnodes,cs_array);

	/*Clean-up*/
	xDelete<int>(cs_array);
}/*}}}*/
void       Element::TransformInvStiffnessMatrixCoord(ElementMatrix* Ke,Node** nodes_list,int numnodes,int* cs_array){/*{{{*/

	int         i,j;
	int         numdofs   = 0;
	IssmDouble *transform = NULL;
	IssmDouble *values    = NULL;

	/*Get total number of dofs*/
	for(i=0;i<numnodes;i++){
		switch(cs_array[i]){
			case PressureEnum: numdofs+=1; break;
			case XYEnum:       numdofs+=2; break;
			case XYZEnum:      numdofs+=3; break;
			default: _error_("Coordinate system " << EnumToStringx(cs_array[i]) << " not supported yet");
		}
	}

	/*Copy current stiffness matrix*/
	values=xNew<IssmDouble>(Ke->nrows*Ke->nrows);
	for(i=0;i<Ke->nrows;i++) for(j=0;j<Ke->nrows;j++) values[i*Ke->nrows+j]=Ke->values[i*Ke->nrows+j];

	/*Get Coordinate Systems transform matrix*/
	CoordinateSystemTransform(&transform,nodes_list,numnodes,cs_array);

	/*Transform matrix: R*Ke*R^T */
	TripleMultiply(transform,numdofs,numdofs,0,
				values,Ke->nrows,Ke->nrows,0,
				transform,numdofs,numdofs,1,
				&Ke->values[0],0);

	/*Free Matrix*/
	xDelete<IssmDouble>(transform);
	xDelete<IssmDouble>(values);
}/*}}}*/
void       Element::TransformLoadVectorCoord(ElementVector* pe,int transformenum){/*{{{*/

	/*All nodes have the same Coordinate System*/
	int  numnodes = this->GetNumberOfNodes();
	int* cs_array = xNew<int>(numnodes);
	for(int i=0;i<numnodes;i++) cs_array[i]=transformenum;

	/*Call core*/
	this->TransformLoadVectorCoord(pe,this->nodes,numnodes,cs_array);

	/*Clean-up*/
	xDelete<int>(cs_array);
}/*}}}*/
void       Element::TransformLoadVectorCoord(ElementVector* pe,int* cs_array){/*{{{*/

	this->TransformLoadVectorCoord(pe,this->nodes,this->GetNumberOfNodes(),cs_array);

}/*}}}*/
void       Element::TransformLoadVectorCoord(ElementVector* pe,Node** nodes_list,int numnodes,int* cs_array){/*{{{*/

	int         i;
	int         numdofs   = 0;
	IssmDouble *transform = NULL;
	IssmDouble *values    = NULL;

	/*Get total number of dofs*/
	for(i=0;i<numnodes;i++){
		switch(cs_array[i]){
			case PressureEnum: numdofs+=1; break;
			case XYEnum:       numdofs+=2; break;
			case XYZEnum:      numdofs+=3; break;
			default: _error_("Coordinate system " << EnumToStringx(cs_array[i]) << " not supported yet");
		}
	}

	/*Copy current load vector*/
	values=xNew<IssmDouble>(pe->nrows);
	for(i=0;i<pe->nrows;i++) values[i]=pe->values[i];

	/*Get Coordinate Systems transform matrix*/
	CoordinateSystemTransform(&transform,nodes_list,numnodes,cs_array);

	/*Transform matrix: R^T*pe */
	MatrixMultiply(transform,numdofs,numdofs,1,
				values,pe->nrows,1,0,
				&pe->values[0],0);

	/*Free Matrices*/
	xDelete<IssmDouble>(transform);
	xDelete<IssmDouble>(values);
}/*}}}*/
void       Element::TransformSolutionCoord(IssmDouble* values,int transformenum){/*{{{*/

	/*All nodes have the same Coordinate System*/
	int  numnodes = this->GetNumberOfNodes();

	/*Do we need to actually do anything? (only if we have a rotated coordinate system)*/
	bool isrotation = false;
	for(int i=0;i<numnodes;i++){
		if(this->nodes[i]->isrotated){
			isrotation=true;
			break;
		}
	}
	if(!isrotation) return;

	/*Rotate solution*/
	int* cs_array = xNew<int>(numnodes);
	for(int i=0;i<numnodes;i++) cs_array[i]=transformenum;
	this->TransformSolutionCoord(values,this->nodes,numnodes,cs_array);

	/*Clean-up*/
	xDelete<int>(cs_array);
}/*}}}*/
void       Element::TransformSolutionCoord(IssmDouble* values,int* transformenum_list){/*{{{*/
	this->TransformSolutionCoord(values,this->nodes,this->GetNumberOfNodes(),transformenum_list);
}/*}}}*/
void       Element::TransformSolutionCoord(IssmDouble* values,int numnodes,int transformenum){/*{{{*/

	/*Do we need to actually do anything? (only if we have a rotated coordinate system)*/
	bool isrotation = false;
	for(int i=0;i<numnodes;i++){
		if(this->nodes[i]->isrotated){
			isrotation=true;
			break;
		}
	}
	if(!isrotation) return;

	/*All nodes have the same Coordinate System*/
	int* cs_array = xNew<int>(numnodes);
	for(int i=0;i<numnodes;i++) cs_array[i]=transformenum;

	/*Call core*/
	this->TransformSolutionCoord(values,this->nodes,numnodes,cs_array);

	/*Clean-up*/
	xDelete<int>(cs_array);
}/*}}}*/
void       Element::TransformSolutionCoord(IssmDouble* solution,int numnodes,int* cs_array){/*{{{*/
	this->TransformSolutionCoord(solution,this->nodes,numnodes,cs_array);
}/*}}}*/
void       Element::TransformSolutionCoord(IssmDouble* values,Node** nodes_list,int numnodes,int transformenum){/*{{{*/
	/*NOT NEEDED*/
	_error_("NOT NEEDED??");
	/*All nodes have the same Coordinate System*/
	int* cs_array = xNew<int>(numnodes);
	for(int i=0;i<numnodes;i++) cs_array[i]=transformenum;

	/*Call core*/
	this->TransformSolutionCoord(values,nodes_list,numnodes,cs_array);

	/*Clean-up*/
	xDelete<int>(cs_array);
}/*}}}*/
void       Element::TransformSolutionCoord(IssmDouble* solution,Node** nodes_list,int numnodes,int* cs_array){/*{{{*/

	int         i;
	int         numdofs   = 0;
	IssmDouble *transform = NULL;
	IssmDouble *values    = NULL;

	/*Get total number of dofs*/
	for(i=0;i<numnodes;i++){
		switch(cs_array[i]){
			case PressureEnum: numdofs+=1; break;
			case XYEnum:       numdofs+=2; break;
			case XYZEnum:      numdofs+=3; break;
			default: _error_("Coordinate system " << EnumToStringx(cs_array[i]) << " not supported yet");
		}
	}

	/*Copy current solution vector*/
	values=xNew<IssmDouble>(numdofs);
	for(int i=0;i<numdofs;i++) values[i]=solution[i];

	/*Get Coordinate Systems transform matrix*/
	CoordinateSystemTransform(&transform,nodes_list,numnodes,cs_array);

	/*Transform matrix: R*U */
	MatrixMultiply(transform,numdofs,numdofs,0,
				values,numdofs,1,0,
				&solution[0],0);

	/*Free Matrices*/
	xDelete<IssmDouble>(transform);
	xDelete<IssmDouble>(values);
}/*}}}*/
void       Element::TransformStiffnessMatrixCoord(ElementMatrix* Ke,int transformenum){/*{{{*/

	/*All nodes have the same Coordinate System*/
	int  numnodes = this->GetNumberOfNodes();

	/*Do we need to actually do anything? (only if we have a rotated coordinate system)*/
	bool isrotation = false;
	for(int i=0;i<numnodes;i++){
		if(this->nodes[i]->isrotated){
			isrotation=true;
			break;
		}
	}
	if(!isrotation) return;

	/*Call core*/
	int* cs_array = xNew<int>(numnodes);
	for(int i=0;i<numnodes;i++) cs_array[i]=transformenum;
	this->TransformStiffnessMatrixCoord(Ke,this->nodes,numnodes,cs_array);

	/*Clean-up*/
	xDelete<int>(cs_array);
}/*}}}*/
void       Element::TransformStiffnessMatrixCoord(ElementMatrix* Ke,int* transformenum_list){/*{{{*/
	this->TransformStiffnessMatrixCoord(Ke,this->nodes,this->GetNumberOfNodes(),transformenum_list);
}/*}}}*/
void       Element::TransformStiffnessMatrixCoord(ElementMatrix* Ke,Node** nodes_list,int numnodes,int* cs_array){/*{{{*/

	int         numdofs = 0;
	IssmDouble *transform = NULL;
	IssmDouble *values    = NULL;

	/*Get total number of dofs*/
	for(int i=0;i<numnodes;i++){
		switch(cs_array[i]){
			case PressureEnum: numdofs+=1; break;
			case XYEnum:       numdofs+=2; break;
			case XYZEnum:      numdofs+=3; break;
			default: _error_("Coordinate system " << EnumToStringx(cs_array[i]) << " not supported yet");
		}
	}

	/*Copy current stiffness matrix*/
	values=xNew<IssmDouble>(Ke->nrows*Ke->nrows);
	for(int i=0;i<Ke->nrows*Ke->nrows;i++) values[i]=Ke->values[i];

	/*Get Coordinate Systems transform matrix*/
	CoordinateSystemTransform(&transform,nodes_list,numnodes,cs_array);

	/*Transform matrix: R^T*Ke*R */
	TripleMultiply(transform,numdofs,numdofs,1,
				values,Ke->nrows,Ke->nrows,0,
				transform,numdofs,numdofs,0,
				&Ke->values[0],0);

	/*Free Matrix*/
	xDelete<IssmDouble>(transform);
	xDelete<IssmDouble>(values);
}/*}}}*/
void       Element::ViscousHeatingCreateInput(void){/*{{{*/

	/*Intermediaries*/
	IssmDouble *xyz_list = NULL;

	/*Fetch number vertices and allocate memory*/
	const int NUM_VERTICES = this->GetNumberOfVertices();
	IssmDouble viscousheating[MAXVERTICES];

	/*Retrieve all inputs and parameters*/
	this->GetVerticesCoordinates(&xyz_list);
	Input* vx_input = this->GetInput(VxEnum); _assert_(vx_input);
	Input* vy_input = this->GetInput(VyEnum); _assert_(vy_input);
	Input* vz_input = this->GetInput(VzEnum); _assert_(vz_input);

	/*loop over vertices: */
	Gauss* gauss=this->NewGauss();
	for(int iv=0;iv<NUM_VERTICES;iv++){
		gauss->GaussVertex(iv);
		this->ViscousHeating(&viscousheating[iv],xyz_list,gauss,vx_input,vy_input,vz_input);
	}

	/*Create PentaVertex input, which will hold the basal friction:*/
	this->AddInput(ViscousHeatingEnum,&viscousheating[0],P1Enum);

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	delete gauss;
}
/*}}}*/

/*Enthalpy*/
void       Element::ThermalToEnthalpy(IssmDouble * penthalpy,IssmDouble temperature,IssmDouble waterfraction,IssmDouble pressure){/*{{{*/

	/*Ouput*/
	IssmDouble enthalpy;

	/*Get necessary parameters*/
	IssmDouble latentheat,referencetemperature,heatcapacity;
	parameters->FindParam(&latentheat,MaterialsLatentheatEnum);
	parameters->FindParam(&referencetemperature,ConstantsReferencetemperatureEnum);
	parameters->FindParam(&heatcapacity,MaterialsHeatcapacityEnum);

	if(temperature<TMeltingPoint(pressure)){
		enthalpy=heatcapacity*(temperature-referencetemperature);
	}
	else{
		enthalpy=PureIceEnthalpy(pressure)+latentheat*waterfraction;
	}

	/*Assign output pointers:*/
	*penthalpy=enthalpy;
}
/*}}}*/
IssmDouble Element::TMeltingPoint(IssmDouble pressure){/*{{{*/

	/*Get necessary parameters*/
	IssmDouble beta,meltingpoint;
	parameters->FindParam(&beta,MaterialsBetaEnum);
	parameters->FindParam(&meltingpoint,MaterialsMeltingpointEnum);

	return meltingpoint-beta*pressure;
}
/*}}}*/
void       Element::EnthalpyToThermal(IssmDouble* ptemperature,IssmDouble* pwaterfraction,IssmDouble enthalpy,IssmDouble pressure){/*{{{*/

	/*Ouput*/
	IssmDouble temperature,waterfraction;

	/*Get necessary parameters*/
	IssmDouble latentheat,referencetemperature,heatcapacity;
	parameters->FindParam(&latentheat,MaterialsLatentheatEnum);
	parameters->FindParam(&referencetemperature,ConstantsReferencetemperatureEnum);
	parameters->FindParam(&heatcapacity,MaterialsHeatcapacityEnum);

	if(enthalpy<PureIceEnthalpy(pressure)){
		temperature=referencetemperature+enthalpy/heatcapacity;
		waterfraction=0.;
	}
	else{
		temperature=TMeltingPoint(pressure);
		waterfraction=(enthalpy-PureIceEnthalpy(pressure))/latentheat;
	}

	/*Assign output pointers:*/
	*pwaterfraction=waterfraction;
	*ptemperature=temperature;
}
/*}}}*/
IssmDouble Element::EnthalpyDiffusionParameter(IssmDouble enthalpy,IssmDouble pressure){/*{{{*/

	/*Get necessary parameters*/
	IssmDouble heatcapacity,thermalconductivity,temperateiceconductivity;
	parameters->FindParam(&heatcapacity,MaterialsHeatcapacityEnum);
	parameters->FindParam(&thermalconductivity,MaterialsThermalconductivityEnum);
	parameters->FindParam(&temperateiceconductivity,MaterialsTemperateiceconductivityEnum);

	if(enthalpy<PureIceEnthalpy(pressure)){
		return thermalconductivity/heatcapacity;
	}
	else{
		return temperateiceconductivity/heatcapacity;
	}
}
/*}}}*/
IssmDouble Element::EnthalpyDiffusionParameterVolume(int numvertices,IssmDouble* enthalpy,IssmDouble* pressure){/*{{{*/

	IssmDouble  lambda;                 // fraction of cold ice
	IssmDouble  kappa,kappa_c,kappa_t;  //enthalpy conductivities
	IssmDouble  Hc,Ht;
	IssmDouble PIE[MAXVERTICES];
	IssmDouble dHpmp[MAXVERTICES];

	for(int iv=0; iv<numvertices; iv++){
		PIE[iv]=PureIceEnthalpy(pressure[iv]);
		dHpmp[iv]=enthalpy[iv]-PIE[iv];
	}

	bool allequalsign=true;
	if(dHpmp[0]<0){
		for(int iv=1; iv<numvertices;iv++) allequalsign=(allequalsign && (dHpmp[iv]<0));
	}
	else{
		for(int iv=1; iv<numvertices;iv++) allequalsign=(allequalsign && (dHpmp[iv]>=0));
	}

	if(allequalsign){
		kappa=EnthalpyDiffusionParameter(enthalpy[0], pressure[0]);
	}
	else {
		/* return harmonic mean of thermal conductivities, weighted by fraction of cold/temperate ice,
			cf Patankar 1980, pp44 */
		kappa_c=EnthalpyDiffusionParameter(PureIceEnthalpy(0.)-1.,0.);
		kappa_t=EnthalpyDiffusionParameter(PureIceEnthalpy(0.)+1.,0.);
		Hc=0.; Ht=0.;
		for(int iv=0; iv<numvertices;iv++){
			if(enthalpy[iv]<PIE[iv])
			 Hc+=(PIE[iv]-enthalpy[iv]);
			else
			 Ht+=(enthalpy[iv]-PIE[iv]);
		}
		_assert_((Hc+Ht)>0.);
		lambda = Hc/(Hc+Ht);
		kappa  = 1./(lambda/kappa_c + (1.-lambda)/kappa_t);
	}

	/*Clean up and return*/
	return kappa;
}
/*}}}*/
IssmDouble Element::PureIceEnthalpy(IssmDouble pressure){/*{{{*/

	/*Get necessary parameters*/
	IssmDouble referencetemperature,heatcapacity;
	parameters->FindParam(&referencetemperature,ConstantsReferencetemperatureEnum);
	parameters->FindParam(&heatcapacity,MaterialsHeatcapacityEnum);

	return heatcapacity*(TMeltingPoint(pressure)-referencetemperature);
}
/*}}}*/

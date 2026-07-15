/*!\file: controlnudging_core.cpp
 * \brief: core of the control solution
 */

#include "./cores.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"

void ComputeRMSEs(FemModel* femmodel, IssmDouble deltat, IssmDouble* pJ){/*{{{*/

	/*output: */
	IssmDouble RMSE_H = 0.;
	IssmDouble RMSE_H_sum = 0.;
	IssmDouble RMSE_dHdt = 0.;
	IssmDouble RMSE_dHdt_sum = 0.;
	IssmDouble RMSE_vel = 0.;
	IssmDouble RMSE_vel_sum = 0.;
	IssmDouble S = 0.;
	IssmDouble S_sum = 1.;

	IssmDouble  weight,Jdet,H,Hobs,Hold;
	IssmDouble  vx, vy, vxobs, vyobs;
	IssmDouble* xyz_list = NULL;

	/*Retrieve parameters*/
   IssmDouble yts = femmodel->parameters->FindParam(ConstantsYtsEnum);

	/*Compute Misfit: */
	for(Object* & object : femmodel->elements->objects){
		Element* element = xDynamicCast<Element*>(object);

		/*If on water, return 0: */
		if(!element->IsOnSurface() || !element->IsIceInElement()) continue;

		/*Spawn surface element*/
		element = element->SpawnTopElement();

		/* Get node coordinates*/
		element->GetVerticesCoordinates(&xyz_list);
		S += element->SurfaceArea();

		/*Retrieve all inputs we will be needing: */
		Input* vx_input   = element->GetInput(VxEnum);             _assert_(vx_input);
		Input* vy_input   = element->GetInput(VyEnum);             _assert_(vy_input);
		Input* vxobs_input= element->GetInput(InversionVxObsEnum); _assert_(vxobs_input);
		Input* vyobs_input= element->GetInput(InversionVyObsEnum); _assert_(vyobs_input);
		Input* H_input    = element->GetInput(ThicknessEnum);      _assert_(H_input);
		Input* Hobs_input = element->GetInput(InversionThicknessObsEnum);        _assert_(Hobs_input);
		Input* Hold_input = element->GetInput(ThicknessPreviousNudgingStepEnum); _assert_(Hold_input);

		/* Start  looping on the number of gaussian points: */
		Gauss* gauss=element->NewGauss(1);
		while(gauss->next()){

			/* Get Jacobian determinant: */
			element->JacobianDeterminant(&Jdet,xyz_list,gauss);

			/*Get all parameters at gaussian point*/
			H_input->GetInputValue(&H, gauss);
			Hobs_input->GetInputValue(&Hobs, gauss);
			Hold_input->GetInputValue(&Hold, gauss);
			vx_input->GetInputValue(&vx, gauss);
			vxobs_input->GetInputValue(&vxobs, gauss);
			vy_input->GetInputValue(&vy, gauss);
			vyobs_input->GetInputValue(&vyobs, gauss);

			RMSE_H    += pow(H-Hobs, 2)*Jdet*gauss->weight;
			RMSE_dHdt += pow( (H - Hold)/deltat*yts, 2)*Jdet*gauss->weight;
			RMSE_vel  += pow( (sqrt(vx*vx+vy*vy) - sqrt(vxobs*vxobs+vyobs*vyobs))*yts, 2)*Jdet*gauss->weight;
		}

		/*clean up and Return: */
		xDelete<IssmDouble>(xyz_list);
		delete gauss;
		if(element->IsSpawnedElement()){element->DeleteMaterials(); delete element;};
	}

	/*Sum all J from all cpus of the cluster:*/
	ISSM_MPI_Reduce (&RMSE_H,   &RMSE_H_sum,    1, ISSM_MPI_DOUBLE, ISSM_MPI_SUM, 0, IssmComm::GetComm());
	ISSM_MPI_Reduce (&RMSE_dHdt,&RMSE_dHdt_sum, 1, ISSM_MPI_DOUBLE, ISSM_MPI_SUM, 0, IssmComm::GetComm());
	ISSM_MPI_Reduce (&RMSE_vel, &RMSE_vel_sum,  1, ISSM_MPI_DOUBLE, ISSM_MPI_SUM, 0, IssmComm::GetComm());
	ISSM_MPI_Reduce (&S,        &S_sum,         1, ISSM_MPI_DOUBLE, ISSM_MPI_SUM, 0, IssmComm::GetComm());
	_printf0_("   → RMSE H   : " << sqrt(RMSE_H_sum/S_sum)   << " m\n");
	_printf0_("   → RMSE dHdt: " << sqrt(RMSE_dHdt_sum/S_sum)<< " m/yr\n");
	_printf0_("   → RMSE v   : " << sqrt(RMSE_vel_sum/S_sum) << " m/yr\n");
	pJ[0] = sqrt(RMSE_H_sum/S_sum);
	pJ[1] = sqrt(RMSE_dHdt_sum/S_sum);
	pJ[2] = sqrt(RMSE_vel_sum/S_sum);

	/*Assign output pointers: */
}/*}}}*/
void controlnudging_core(FemModel* femmodel){

   /*Intermediaries*/
	IssmDouble dCdt1, dCdt2, dCdt3, dHdt;
	IssmDouble dMeltdt1, dMeltdt2, dMeltdt3;
   IssmDouble time;
	int        maxiter;

   /*User-defined nudging parameters*/
   femmodel->parameters->FindParam(&maxiter, InversionMaxiterEnum);
   IssmDouble C0           = femmodel->parameters->FindParam(InversionC0Enum);
	IssmDouble melt0        = femmodel->parameters->FindParam(InversionMelt0Enum);
   IssmDouble tau_C        = femmodel->parameters->FindParam(InversionTauCEnum);
	IssmDouble tau_melt     = femmodel->parameters->FindParam(InversionTauMeltEnum);
   IssmDouble max_inc_C    = femmodel->parameters->FindParam(InversionMaxIncrementCEnum);
   IssmDouble max_inc_melt = femmodel->parameters->FindParam(InversionMaxIncrementMeltEnum);
   IssmDouble H0_C         = femmodel->parameters->FindParam(InversionH0CEnum);
   IssmDouble H0_melt      = femmodel->parameters->FindParam(InversionH0MeltEnum);
   IssmDouble r_C          = femmodel->parameters->FindParam(InversionRelaxationCEnum);
   IssmDouble r_melt       = femmodel->parameters->FindParam(InversionRelaxationMeltEnum);
   IssmDouble yts          = femmodel->parameters->FindParam(ConstantsYtsEnum);
   IssmDouble tmax         = femmodel->parameters->FindParam(TimesteppingFinalTimeEnum);
   IssmDouble tmin         = femmodel->parameters->FindParam(TimesteppingStartTimeEnum);
   IssmDouble deltat       = (tmax - tmin);

   /*Fields before/after*/
	IssmDouble *C       = NULL;
	IssmDouble *Cinit   = NULL;
	IssmDouble *Cmin    = NULL;
	IssmDouble *Cmax    = NULL;
	IssmDouble *Melt    = NULL;
	IssmDouble *Meltmin = NULL;
	IssmDouble *Meltmax = NULL;
	IssmDouble *H_obs   = NULL;
	IssmDouble *H_old   = NULL;
	IssmDouble *H       = NULL;
	IssmDouble *V       = NULL;
	IssmDouble *O_ls    = NULL;
	
	/*Cost functions*/
	IssmDouble* J = xNewZeroInit<IssmDouble>(maxiter*3);

   /*Get Fields once and for all*/
   int numvertices = femmodel->vertices->NumberOfVertices();
	GetVectorFromInputsx(&Melt, femmodel, BasalforcingsPerturbationMeltingRateEnum, VertexSIdEnum);
	GetVectorFromInputsx(&C,    femmodel, FrictionCoefficientEnum, VertexSIdEnum);
   GetVectorFromInputsx(&Cinit,femmodel, FrictionCoefficientEnum, VertexSIdEnum);
	GetVectorFromInputsx(&Cmin, femmodel, InversionMinCEnum, VertexSIdEnum);
	GetVectorFromInputsx(&Cmax, femmodel, InversionMaxCEnum, VertexSIdEnum);
	GetVectorFromInputsx(&Meltmin,femmodel, InversionMinMeltEnum, VertexSIdEnum);
	GetVectorFromInputsx(&Meltmax,femmodel, InversionMaxMeltEnum, VertexSIdEnum);
   GetVectorFromInputsx(&H_obs,  femmodel, ThicknessEnum, VertexSIdEnum);

	/*NEW*/
	InputDuplicatex(femmodel, ThicknessEnum, InversionThicknessObsEnum);

   femmodel->parameters->SetParam(true, DoNotSaveResultsEnum);
   for(int m=0;m<maxiter;m++){
		int size = femmodel->Size();
      _printf0_("\n=== NUDGING STEP "<< m+1 <<"/"<< maxiter << " === SIZE = "<<size<<"\n");

      /*Get ice thickness before we run a transient step*/
      GetVectorFromInputsx(&H_old, femmodel, ThicknessEnum, VertexSIdEnum);
		InputDuplicatex(femmodel, ThicknessEnum, ThicknessPreviousNudgingStepEnum);

      /*Solve another transient simulation*/
      femmodel->parameters->FindParam(&time,TimeEnum);
      femmodel->parameters->SetParam(time+deltat,TimesteppingFinalTimeEnum);
		if(m==maxiter-1){
			femmodel->parameters->SetParam(false, DoNotSaveResultsEnum);
		}
      transient_core(femmodel);

      /*Extract results*/
      GetVectorFromInputsx(&H,   femmodel, ThicknessEnum, VertexSIdEnum);
		GetVectorFromInputsx(&V,   femmodel, VelEnum, VertexSIdEnum);
      GetVectorFromInputsx(&O_ls,femmodel, MaskOceanLevelsetEnum, VertexSIdEnum);

      /*Update friction coefficient accordingly*/
      for(int i=0;i<numvertices;i++){

         /*Compute thickness change for this vertex*/
         IssmDouble dH_now   = H[i] - H_obs[i];
         IssmDouble dHdt_now = (H[i] - H_old[i])/(deltat);

         /*1. : thickness error — push C to reduce H deviation
          *     Sign: if H > H_obs (too thick), decrease C (less friction → faster ice → larger flux, lower H)*/
         dCdt1 = -C0*(dH_now/H0_C) / tau_C;

         /*2. : tendency — damp ongoing thinning/thickening
          *     Sign: if dH/dt > 0 (thickening), if it is thickening and already too thick, decrease c faster.
          *     If it is thinning and the thickness is already to large, the ice is already moving in the correct
          *     direction, so make C decrease a little bit less */
         dCdt2 = -C0*(dHdt_now/H0_C);

         /*3. Relaxation term*/
         IssmDouble C_log  = log10(max(C[i],  1.));
         IssmDouble C_log0 = log10(max(Cinit[i], 1.));
         dCdt3 = C0*(r_C / tau_C)*(C_log0 - C_log);

			/* Do not nudge C where velocity ratio is very low AND ice is too thin
			 * These cells need upstream replenishment, not local friction changes*/
			//if(V[i]<15./yts && dHdt<-50/yts && O_ls[i]>0.) {
			//	dCdt2 = 0.;
			//	dCdt3 = 0.;
			//}

			/*Compute total dC in log space*/
			IssmDouble dC_log = deltat*(dCdt1 + dCdt2 + dCdt3);
         if(dC_log> max_inc_C) dC_log = max_inc_C;
         if(dC_log<-max_inc_C) dC_log = -max_inc_C;

         /*Correction: prevent C from decreasing where ice is too thin and
          * still thinning — the model needs more friction here, not less*/
         //if( ((H[i] - H_obs[i])< -20.) && (dHdt<-1.) ){
         //   dC_log = max(dC_log, 0.);
         //}

         /*Update friction coefficient now*/
         C[i] = pow(10., C_log + dC_log);
			if(C[i] > Cmax[i]) C[i] = Cmax[i];
			if(C[i] < Cmin[i]) C[i] = Cmin[i];

			/*Melt update if floating*/
			//if(O_ls[i]<0.){
				/*1. : thickness error — push melt to reduce H deviation
				 *     Sign: if H > H_obs (too thick), increase melt */
				dMeltdt1 = +melt0*(dH_now/H0_melt) / tau_melt;

				/*2. : tendency — damp ongoing thinning/thickening
				 *     Sign: if dH/dt > 0 (thickening), if it is thickening and already too thick, increase melt*/
				dMeltdt2 = +melt0*(dHdt_now/H0_melt);

				/*3. Relaxation term*/
				dMeltdt3 = -melt0*(r_melt / tau_melt)*(Melt[i]);

				/*Compute total dMelt by combining all 3 contributions*/
				IssmDouble dMelt = deltat*(dMeltdt1 + dMeltdt2 + dMeltdt3);
				if(dMelt> max_inc_melt) dMelt = max_inc_melt;
				if(dMelt<-max_inc_melt) dMelt = -max_inc_melt;

				/*Update Melt*/
				Melt[i] += dMelt;
				if(Melt[i] > Meltmax[i]) Melt[i] = Meltmax[i];
				if(Melt[i] < Meltmin[i]) Melt[i] = Meltmin[i];
			//}
      }

      InputUpdateFromVectorx(femmodel, C,    FrictionCoefficientEnum, VertexSIdEnum);
		InputUpdateFromVectorx(femmodel, Melt, BasalforcingsPerturbationMeltingRateEnum, VertexSIdEnum);

		/*Print statistics*/
		ComputeRMSEs(femmodel, deltat,&J[3*m]);

      xDelete<IssmDouble>(H_old);
      xDelete<IssmDouble>(H);
		xDelete<IssmDouble>(V);
      xDelete<IssmDouble>(O_ls);
   }

	/*Add C/melt/J to results*/
	femmodel->results->AddObject(new GenericExternalResult<IssmDouble*>(femmodel->results->Size()+1,
					FrictionCoefficientEnum,C, numvertices, 1, 0, 0));
	femmodel->results->AddObject(new GenericExternalResult<IssmDouble*>(femmodel->results->Size()+1,
					BasalforcingsPerturbationMeltingRateEnum,Melt, numvertices, 1, 0, 0));
	femmodel->results->AddObject(new GenericExternalResult<IssmDouble*>(femmodel->results->Size()+1,JEnum,J, maxiter, 3,0,0));

   /*Clean up and return*/
	xDelete<IssmDouble>(C);
   xDelete<IssmDouble>(Cinit);
	xDelete<IssmDouble>(Melt);
	xDelete<IssmDouble>(Meltmin);
   xDelete<IssmDouble>(Meltmax);
	xDelete<IssmDouble>(Cmin);
	xDelete<IssmDouble>(Cmax);
   xDelete<IssmDouble>(H_obs);
	xDelete<IssmDouble>(J);
}

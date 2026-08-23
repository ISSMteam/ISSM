/*!\file: controlnudging_core.cpp
 * \brief: core of the control solution
 */

#include "./cores.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"

void controlnudging_core(FemModel* femmodel){

   /*Intermediaries*/
	IssmDouble dCdt1, dCdt2, dCdt3, dHdt;
	IssmDouble dMeltdt1, dMeltdt2, dMeltdt3;
	IssmDouble RMSE_dHdt, RMSE_H, RMSE_dHdt_sum, RMSE_H_sum;
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
	IssmDouble *Cinit      = NULL;
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

   /*Get Fields once and for all*/
   int numvertices       = femmodel->vertices->NumberOfVertices();
	int numvertices_local = femmodel->vertices->NumberOfVerticesLocalAll();
	GetVectorFromInputsx(&Melt, femmodel, BasalforcingsPerturbationMeltingRateEnum, VertexLIdEnum);
   GetVectorFromInputsx(&Cinit, femmodel, FrictionCoefficientEnum, VertexLIdEnum);
	GetVectorFromInputsx(&Cmin,femmodel, InversionMinCEnum, VertexLIdEnum);
	GetVectorFromInputsx(&Cmax,femmodel, InversionMaxCEnum, VertexLIdEnum);
	GetVectorFromInputsx(&Meltmin,femmodel, InversionMinMeltEnum, VertexLIdEnum);
	GetVectorFromInputsx(&Meltmax,femmodel, InversionMaxMeltEnum, VertexLIdEnum);
   GetVectorFromInputsx(&H_obs, femmodel, ThicknessEnum, VertexLIdEnum);

   femmodel->parameters->SetParam(false,SaveResultsEnum);
   for(int m=0;m<maxiter;m++){
      _printf0_("\n=== NUDGING STEP "<< m+1 <<"/"<< maxiter << " ===\n");

      /*Get ice thickness before we run a transient step*/
      GetVectorFromInputsx(&H_old, femmodel, ThicknessEnum, VertexLIdEnum);

      /*Solve another transient simulation*/
      femmodel->parameters->FindParam(&time,TimeEnum);
      femmodel->parameters->SetParam(time+deltat,TimesteppingFinalTimeEnum);
      transient_core(femmodel);

      /*Extract results*/
		xDelete<IssmDouble>(C);
      GetVectorFromInputsx(&H,   femmodel, ThicknessEnum, VertexLIdEnum);
      GetVectorFromInputsx(&C,   femmodel, FrictionCoefficientEnum, VertexLIdEnum);
		GetVectorFromInputsx(&V,   femmodel, VelEnum, VertexLIdEnum);
      GetVectorFromInputsx(&O_ls,femmodel, MaskOceanLevelsetEnum, VertexLIdEnum);

      /*Update friction coefficient accordingly*/
		RMSE_H    = 0.;
		RMSE_dHdt = 0.;
      for(int index=0;index<numvertices_local;index++){

			/*Convert index to lid since this is how we look at vectors*/
			Vertex* vertex = xDynamicCast<Vertex*>(femmodel->vertices->GetObjectByOffset(index));
			int i = vertex->Lid();

         /*Compute thickness change for this vertex*/
         IssmDouble dH_now   = H[i] - H_obs[i];
         IssmDouble dHdt_now = (H[i] - H_old[i])/(deltat);


			if(!vertex->clone){
				RMSE_H    += dH_now*dH_now;
				RMSE_dHdt += pow(dHdt_now*yts, 2); //Convert to m/yr
			}

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

      InputUpdateFromVectorx(femmodel, C,    FrictionCoefficientEnum, VertexLIdEnum);
		InputUpdateFromVectorx(femmodel, Melt, BasalforcingsPerturbationMeltingRateEnum, VertexLIdEnum);

		/*Print statistics*/
		ISSM_MPI_Reduce(&RMSE_H,    &RMSE_H_sum, 1, ISSM_MPI_DOUBLE, ISSM_MPI_SUM, 0, IssmComm::GetComm() );
		ISSM_MPI_Reduce(&RMSE_dHdt, &RMSE_dHdt_sum, 1, ISSM_MPI_DOUBLE, ISSM_MPI_SUM, 0, IssmComm::GetComm() );
		_printf0_("   → RMSE H   : " << sqrt(RMSE_H_sum/numvertices)    << " m\n");
		_printf0_("   → RMSE dHdt: " << sqrt(RMSE_dHdt_sum/numvertices) << " m/yr\n");

      xDelete<IssmDouble>(H_old);
      xDelete<IssmDouble>(H);
      xDelete<IssmDouble>(O_ls);
		xDelete<IssmDouble>(V);
   }

	/*Add C to results*/
	femmodel->results->AddObject(new GenericExternalResult<IssmDouble*>(femmodel->results->Size()+1,
					FrictionCoefficientEnum,C, numvertices, 1, 0, 0));
	femmodel->results->AddObject(new GenericExternalResult<IssmDouble*>(femmodel->results->Size()+1,
					BasalforcingsPerturbationMeltingRateEnum,Melt, numvertices, 1, 0, 0));

   /*Clean up and return*/
	xDelete<IssmDouble>(C);
   xDelete<IssmDouble>(Cinit);
	xDelete<IssmDouble>(Melt);
	xDelete<IssmDouble>(Meltmin);
   xDelete<IssmDouble>(Meltmax);
	xDelete<IssmDouble>(Cmin);
	xDelete<IssmDouble>(Cmax);
   xDelete<IssmDouble>(H_obs);
}

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
   int numvertices=femmodel->vertices->NumberOfVertices();

   /*Hard-coded nudging parameters*/
   IssmDouble timeadjust = 1/6;   // yr per nudging step (2 months)
   IssmDouble H0         = 100;   // m — thickness error scale (smaller = more sensitive), 100 m used in van den Akker et al (2025)
   IssmDouble r          = 0.4;   // relaxation strength toward C_inv (0  = none, 1 = strong) , 0.5 used in van den Akker et al (2025)
   int        nsteps     = 12000; // number of nudging steps

   /*User-defined nudging parameters*/
   IssmDouble tau_C  = femmodel->parameters->FindParam(InversionTauCEnum);
   IssmDouble yts    = femmodel->parameters->FindParam(ConstantsYtsEnum);
   IssmDouble tmax   = femmodel->parameters->FindParam(TimesteppingFinalTimeEnum);
   IssmDouble tmin   = femmodel->parameters->FindParam(TimesteppingStartTimeEnum);
   IssmDouble deltat = (tmax - tmin)*yts; // convert to years

   /*Fields before/after*/
   IssmDouble* C0   = NULL;
   IssmDouble* C    = NULL;
   IssmDouble* Hobs = NULL;
   IssmDouble* H    = NULL;
	IssmDouble* V    = NULL;
   IssmDouble* O_ls = NULL;

   /*Get Fields once and for all*/
   GetVectorFromInputsx(&C0, femmodel, FrictionCoefficientEnum, VertexSIdEnum);
   GetVectorFromInputsx(&Hobs, femmodel, ThicknessEnum, VertexSIdEnum);

   femmodel->parameters->SetParam(false,SaveResultsEnum);
   for(int m=0;m<nsteps;m++){
      _printf0_("\n=== NUDGING STEP "<< m <<"/"<< nsteps << " ===\n");

      /*we need to make sure we do not modify femmodel at each iteration, make a copy*/
      FemModel* femmodel_temp = femmodel->copy();

      /*Solve another transient simulation*/
      transient_core(femmodel);

      /*Extract results*/
      GetVectorFromInputsx(&H, femmodel, ThicknessEnum, VertexSIdEnum);
      GetVectorFromInputsx(&C, femmodel, FrictionCoefficientEnum, VertexSIdEnum);
		GetVectorFromInputsx(&V, femmodel, VelEnum, VertexSIdEnum);
      GetVectorFromInputsx(&O_ls, femmodel, MaskOceanLevelsetEnum, VertexSIdEnum);

      /*Update friction coefficient accordingly*/
      for(int i=0;i<numvertices;i++){

         /*1. : thickness error — push C to reduce H deviation
          *     Sign: if H > H_obs (too thick), decrease C (less friction → faster ice → larger flux, lower H)*/
         dCdt1 = -1*((H[i] - Hobs[i])/ H0) / tau_C;

         /*2. : tendency — damp ongoing thinning/thickening
          *     Sign: if dH/dt > 0 (thickening), if it is thickening and already too thick, decrease c faster.
          *     If it is thinning and the thickness is already to large, the ice is already moving in the correct
          *     direction, so make C decrease a little bit less */
			dHdt = (H[i] - Hobs[i])/deltat;
         dCdt2 = -1*(dHdt/H0) /tau_C;

         /*3. Relaxation term*/
         IssmDouble C_loginv  = log10(max(C[i],  1.));
         IssmDouble C_loginv0 = log10(max(C0[i], 1.));
         dCdt3 = 1*(r / tau_C)*(C_loginv-C_loginv0);

			/* Do not nudge C where velocity ratio is very low AND ice is too thin
			 * These cells need upstream replenishment, not local friction changes*/
			if(V[i]<15./yts && dHdt<-50/yts && O_ls[i]>0.) {
				dCdt2 = 0.;
				dCdt3 = 0.;
			}

			/*Compute total dC in log space*/
			IssmDouble dC_log = deltat*(dCdt1 + dCdt2 + dCdt3);

         /*Correction: prevent C from decreasing where ice is too thin and
          * still thinning — the model needs more friction here, not less*/
         if( ((H[i] - Hobs[i])< -20.) && (dHdt<-1.) ){
            dC_log = max(dC_log, 0.);
         }

         /*Update friction coefficient now*/
         C[i] = pow(10., C_loginv + dC_log);
      }

      InputUpdateFromVectorx(femmodel, C, FrictionCoefficientEnum, VertexSIdEnum);
      _error_("not implemented yet");

      xDelete<IssmDouble>(H);
      xDelete<IssmDouble>(C);
      xDelete<IssmDouble>(O_ls);
		xDelete<IssmDouble>(V);
      delete femmodel_temp;
   }

   /*Clean up and return*/
   xDelete<IssmDouble>(C0);
   xDelete<IssmDouble>(Hobs);
}

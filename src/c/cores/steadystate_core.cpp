/*!\file: steadystate_core.cpp
 * \brief: core of the steadystate solution 
 */ 

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif
#include "./cores.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"
#include "../solutionsequences/solutionsequences.h"

/*Local prototypes*/
bool steadystateconvergence(Vector<IssmDouble>* tg,Vector<IssmDouble>* tg_old,Vector<IssmDouble>* ug,Vector<IssmDouble>* ug_old,IssmDouble reltol);

void steadystate_core(FemModel* femmodel){ //{{{

	/*intermediary: */
	int step; 
	Vector<IssmDouble>* ug     = NULL;
	Vector<IssmDouble>* ug_old = NULL;
	Vector<IssmDouble>* tg     = NULL;
	Vector<IssmDouble>* tg_old = NULL;

	/*parameters: */
	bool        save_results,isenthalpy;
	int         maxiter;
	IssmDouble  reltol;
	int         numoutputs   = 0;
	char** requested_outputs = NULL;

	/* recover parameters:*/
	femmodel->parameters->FindParam(&save_results,SaveResultsEnum);
	femmodel->parameters->FindParam(&maxiter,SteadystateMaxiterEnum);
	femmodel->parameters->FindParam(&isenthalpy,ThermalIsenthalpyEnum);
	femmodel->parameters->FindParam(&reltol,SteadystateReltolEnum);
	femmodel->parameters->SetParam(false,SaveResultsEnum);
	femmodel->parameters->FindParam(&numoutputs,SteadystateNumRequestedOutputsEnum);
	if(numoutputs) femmodel->parameters->FindParam(&requested_outputs,&numoutputs,SteadystateRequestedOutputsEnum);

	/*intialize counters: */
	step=1;

	for(;;){

		/* Compute first velocity, then temperature due to high sensitivity of temperature to velocity. */
		if(VerboseSolution()) _printf0_("\n==================================================\n");
		if(VerboseSolution()) _printf0_("   computing velocity and temperature for step: " << step << "\n");
		if(VerboseSolution()) _printf0_("====================================================\n");

		if(VerboseSolution()) _printf0_("\n   -- computing new velocity -- \n\n");
		stressbalance_core(femmodel);
		GetSolutionFromInputsx(&ug,femmodel);

		if(VerboseSolution()) _printf0_("\n   -- computing new temperature --\n\n");
		thermal_core(femmodel);
		if(!isenthalpy)femmodel->SetCurrentConfiguration(ThermalAnalysisEnum);/*Could be MeltingAnalysis...*/
		GetSolutionFromInputsx(&tg,femmodel);

		if(step>1){
			if(VerboseSolution()) _printf0_("   checking steadystate convergence\n");
			if(steadystateconvergence(tg,tg_old,ug,ug_old,reltol)) break;
		}
		if(step>=maxiter){
			if(VerboseSolution()) _printf0_("   maximum number steadystate iterations " << maxiter << " reached\n");
			break;
		}

		/*update results and increase counter*/
		delete tg_old;tg_old=tg;
		delete ug_old;ug_old=ug;
		step++;
	}

	if(save_results){
		femmodel->RequestedOutputsx(&femmodel->results,requested_outputs,numoutputs);
	}

	/*Free resources:*/
	delete tg_old;
	delete ug_old;
	delete tg;
	delete ug;	
	if(numoutputs){for(int i=0;i<numoutputs;i++){xDelete<char>(requested_outputs[i]);} xDelete<char*>(requested_outputs);}
}//}}}
bool steadystateconvergence(Vector<IssmDouble>* tg,Vector<IssmDouble>* tg_old,Vector<IssmDouble>* ug,Vector<IssmDouble>* ug_old,IssmDouble reltol){//{{{

	/*Output*/
	bool converged = true;

	/*Intermediary*/
	Vector<IssmDouble>* dug    = NULL;
	Vector<IssmDouble>* dtg    = NULL;
	IssmDouble          ndt,nt;
	IssmDouble          ndu,nu;

	/*compute norm(du)/norm(u)*/
	dug=ug_old->Duplicate(); ug_old->Copy(dug); dug->AYPX(ug,-1.0);
	ndu=dug->Norm(NORM_TWO); nu=ug_old->Norm(NORM_TWO);
	if (xIsNan<IssmDouble>(ndu) || xIsNan<IssmDouble>(nu)) _error_("convergence criterion is NaN!");
	if((ndu/nu)<reltol){
		if(VerboseConvergence()) _printf0_("\n"<<setw(50)<<left<<"   Velocity convergence: norm(du)/norm(u)"<<ndu/nu*100<<" < "<<reltol*100<<" %\n");
	}
	else{ 
		if(VerboseConvergence()) _printf0_("\n"<<setw(50)<<left<<"   Velocity convergence: norm(du)/norm(u)"<<ndu/nu*100<<" > "<<reltol*100<<" %\n");
		converged=false;
	}

	/*compute norm(dt)/norm(t)*/
	dtg=tg_old->Duplicate(); tg_old->Copy(dtg); dtg->AYPX(tg,-1.0);
	ndt=dtg->Norm(NORM_TWO); nt=tg_old->Norm(NORM_TWO);
	if (xIsNan<IssmDouble>(ndt) || xIsNan<IssmDouble>(nt)) _error_("convergence criterion is NaN!");
	if((ndt/nt)<reltol){
		if(VerboseConvergence()) _printf0_(setw(50)<<left<<"   Temperature convergence: norm(dt)/norm(t)"<<ndt/nt*100<<" < "<<reltol*100<<" %\n");
	}
	else{ 
		if(VerboseConvergence()) _printf0_(setw(50)<<left<<"   Temperature convergence: norm(dt)/norm(t)"<<ndt/nt*100<<" > "<<reltol*100<<" %\n");
		converged=false;
	}

	/*clean up and return*/
	delete dtg;
	delete dug;
	return converged;
}//}}}

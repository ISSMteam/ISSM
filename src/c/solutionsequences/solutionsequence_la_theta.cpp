/*!\file: solutionsequence_la_theta.cpp
 * \brief: numerical core of la_theta solutions
 */ 

#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"

void solutionsequence_la_theta(FemModel* femmodel){

	/*intermediary: */
	Matrix<IssmDouble>*  Kff    = NULL;
	Matrix<IssmDouble>*  Kfs    = NULL;
	Vector<IssmDouble>*  ug_old = NULL;
	Vector<IssmDouble>*  ug     = NULL;
	Vector<IssmDouble>*  uf     = NULL;
	Vector<IssmDouble>*  pf     = NULL;
	Vector<IssmDouble>*  df     = NULL;
	Vector<IssmDouble>*  ys     = NULL;
	IssmDouble           eps_rel,r,theta; // 0<theta<.5   -> .15<theta<.45
	int                  max_nonlinear_iterations;

	/*Create analysis*/
	StressbalanceAnalysis* analysis = new StressbalanceAnalysis();

	/*Recover parameters: */
	femmodel->parameters->FindParam(&eps_rel,StressbalanceReltolEnum);
	femmodel->parameters->FindParam(&max_nonlinear_iterations,StressbalanceMaxiterEnum);
	femmodel->parameters->FindParam(&r,AugmentedLagrangianREnum);
	femmodel->parameters->FindParam(&theta,AugmentedLagrangianThetaEnum);

	/*Update constraints and initialize d and tau if necessary*/
	femmodel->UpdateConstraintsx();
	analysis->InitializeXTH(femmodel->elements,femmodel->parameters);

	/*Convergence criterion*/
	int  count = 0;
	GetSolutionFromInputsx(&ug,femmodel);
	Vector<IssmDouble>* vx     = NULL;
	Vector<IssmDouble>* vx_old = NULL;
	GetVectorFromInputsx(&vx,femmodel,VxEnum,VertexPIdEnum);

	while(true){
		count++;

		/*save pointer to old velocity*/
		delete ug_old;ug_old=ug;
		delete vx_old;vx_old=vx;

		/*Calculate d*/
		if(theta>0.){
			analysis->InputUpdateFromSolutionFSXTH_d(femmodel->elements,femmodel->parameters);
		}

		/*Solve KU=F*/
		SystemMatricesx(&Kff,&Kfs,&pf,&df,NULL,femmodel);
		CreateNodalConstraintsx(&ys,femmodel->nodes);
		Reduceloadx(pf, Kfs, ys); delete Kfs;

		femmodel->profiler->Start(SOLVER);
		Solverx(&uf, Kff, pf, NULL, df, femmodel->parameters); 
		femmodel->profiler->Stop(SOLVER);

		delete Kff; delete pf; delete df;
		Mergesolutionfromftogx(&ug, uf,ys,femmodel->nodes,femmodel->parameters);delete uf; delete ys;

		/*Update solution*/
		InputUpdateFromSolutionx(femmodel,ug); 

		/*Update d and tau accordingly*/
		analysis->InputUpdateFromSolutionFSXTH_d(  femmodel->elements,femmodel->parameters);
		analysis->InputUpdateFromSolutionFSXTH_tau(femmodel->elements,femmodel->parameters);
		GetVectorFromInputsx(&vx,femmodel,VxEnum,VertexPIdEnum);

		/*Check for convergence*/
		//Vector<IssmDouble>* dug=ug_old->Duplicate(); ug_old->Copy(dug); dug->AYPX(ug,-1.0);
		//IssmDouble ndu=dug->Norm(NORM_TWO);   delete dug;
		//IssmDouble nu =ug_old->Norm(NORM_TWO);
		Vector<IssmDouble>* dvx=vx_old->Duplicate(); vx_old->Copy(dvx); dvx->AYPX(vx,-1.0);
		IssmDouble ndu=dvx->Norm(NORM_TWO);   delete dvx;
		IssmDouble nu =vx_old->Norm(NORM_TWO);
		if (xIsNan<IssmDouble>(ndu) || xIsNan<IssmDouble>(nu)) _error_("convergence criterion is NaN!");
		if((ndu/nu)<eps_rel){
			if(VerboseConvergence()) _printf0_(setw(50) << left << "   Convergence criterion: norm(du)/norm(u)" << ndu/nu*100 << " < " << eps_rel*100 << " %\n");
			break;
		}
		else{ 
			if(VerboseConvergence()) _printf0_(setw(50) << left << "   Convergence criterion: norm(du)/norm(u)" << ndu/nu*100 << " > " << eps_rel*100 << " %\n");
		}

		if(count>=max_nonlinear_iterations){
			_printf0_("   maximum number of nonlinear iterations (" << max_nonlinear_iterations << ") exceeded\n"); 
			break;
		}
	}

	if(VerboseConvergence()) _printf0_("\n   total number of iterations: " << count-1 << "\n");

	delete ug;  
	delete ug_old;  
	delete vx;  
	delete vx_old;  
	delete analysis;
}

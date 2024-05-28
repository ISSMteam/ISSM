/*!\file: solutionsequence_la.cpp
 * \brief: numerical core of la solutions
 */ 

#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"

void solutionsequence_la(FemModel* femmodel){

	/*intermediary: */
	Matrix<IssmDouble>*  Kff     = NULL;
	Matrix<IssmDouble>*  Kfs     = NULL;
	Vector<IssmDouble>*  ug      = NULL;
	Vector<IssmDouble>*  uf      = NULL;
	Vector<IssmDouble>*  pf      = NULL;
	Vector<IssmDouble>*  df      = NULL;
	Vector<IssmDouble>*  ys      = NULL;
	Vector<IssmDouble>*  pug     = NULL;
	Vector<IssmDouble>*  pug_old = NULL;
	IssmDouble           eps_rel,r,theta; // 0<theta<.5   -> .15<theta<.45
	int                  max_nonlinear_iterations;
	bool                 vel_converged      = false;
	bool                 pressure_converged = false;

	/*Create analysis*/
	StressbalanceAnalysis* stressanalysis   = new StressbalanceAnalysis();
	UzawaPressureAnalysis* pressureanalysis = new UzawaPressureAnalysis();
	femmodel->SetCurrentConfiguration(StressbalanceAnalysisEnum);

	/*Recover parameters: */
	femmodel->parameters->FindParam(&eps_rel,StressbalanceReltolEnum);
	femmodel->parameters->FindParam(&max_nonlinear_iterations,StressbalanceMaxiterEnum);

	/*Update constraints*/
	femmodel->UpdateConstraintsx();

	/*Convergence criterion*/
	int  count       = 0;
	int  count_local = 0;
	Vector<IssmDouble>* vel           = NULL;
	Vector<IssmDouble>* vel_old       = NULL;
	Vector<IssmDouble>* vel_old_local = NULL;
	GetVectorFromInputsx(&vel,femmodel,VxEnum,VertexPIdEnum);
	GetVectorFromInputsx(&pug,femmodel,PressureEnum,VertexPIdEnum);

	while(true){
		count++;

		/*save pointer to old velocity*/
		delete vel_old;vel_old=vel->Duplicate(); vel->Copy(vel_old);
		delete pug_old;pug_old=pug; 
		pressure_converged=false; vel_converged=false;

		while(true){
			count_local++;
			delete vel_old_local;vel_old_local=vel;
			/*Solve KU=F*/
			femmodel->SetCurrentConfiguration(StressbalanceAnalysisEnum);
			SystemMatricesx(&Kff,&Kfs,&pf,&df,NULL,femmodel);
			CreateNodalConstraintsx(&ys,femmodel->nodes);
			Reduceloadx(pf, Kfs, ys); delete Kfs;

			femmodel->profiler->Start(SOLVER);
			Solverx(&uf, Kff, pf, NULL, df, femmodel->parameters); 
			femmodel->profiler->Stop(SOLVER);

			delete Kff; delete pf; delete df;
			Mergesolutionfromftogx(&ug, uf,ys,femmodel->nodes,femmodel->parameters);delete uf; delete ys;

			/*Update solution*/
			InputUpdateFromSolutionx(femmodel,ug); delete ug;
			GetVectorFromInputsx(&vel,femmodel,VxEnum,VertexPIdEnum);
			/*Check for convergence*/
			Vector<IssmDouble>* dvel_local=vel_old_local->Duplicate(); vel_old_local->Copy(dvel_local); dvel_local->AYPX(vel,-1.0);
			IssmDouble ndu=dvel_local->Norm(NORM_TWO);   delete dvel_local;
			IssmDouble nu =vel_old_local->Norm(NORM_TWO);
			if (xIsNan<IssmDouble>(ndu) || xIsNan<IssmDouble>(nu)) _error_("convergence criterion is NaN!");
			if((ndu/nu)<eps_rel){
				if(VerboseConvergence()) _printf0_(setw(50) << left << "   Convergence criterion: norm(du)/norm(u)" << ndu/nu*100 << " < " << eps_rel*100 << " %\n");
				break;
			}
			else{ 
				if(VerboseConvergence()) _printf0_(setw(50) << left << "   Convergence criterion: norm(du)/norm(u)" << ndu/nu*100 << " > " << eps_rel*100 << " %\n");
			}
			if(count_local>=max_nonlinear_iterations){
				_printf0_("   maximum number of nonlinear iterations (" << max_nonlinear_iterations << ") exceeded\n"); 
				break;
			}
		}
		count_local=0;

		femmodel->SetCurrentConfiguration(UzawaPressureAnalysisEnum);
		SystemMatricesx(&Kff,&Kfs,&pf,&df,NULL,femmodel);
		CreateNodalConstraintsx(&ys,femmodel->nodes);
		Reduceloadx(pf, Kfs, ys); delete Kfs;

		femmodel->profiler->Start(SOLVER);
		Solverx(&uf, Kff, pf, NULL, df, femmodel->parameters); 
		femmodel->profiler->Stop(SOLVER);

		delete Kff; delete pf; delete df;
		Mergesolutionfromftogx(&pug, uf,ys,femmodel->nodes,femmodel->parameters); delete uf; delete ys;

		/*Update solution*/
		InputUpdateFromSolutionx(femmodel,pug); delete pug;
		GetVectorFromInputsx(&pug,femmodel,PressureEnum,VertexPIdEnum);

		/*Check for convergence*/
		Vector<IssmDouble>* dvel=vel_old_local->Duplicate(); vel_old_local->Copy(dvel); dvel->AYPX(vel,-1.0);
		IssmDouble ndu=dvel->Norm(NORM_TWO);   delete dvel;
		IssmDouble nu =vel_old_local->Norm(NORM_TWO);
		if (xIsNan<IssmDouble>(ndu) || xIsNan<IssmDouble>(nu)) _error_("convergence criterion is NaN!");
		if((ndu/nu)<eps_rel){
			if(VerboseConvergence()) _printf0_(setw(50) << left << "   Convergence criterion: norm(du)/norm(u)" << ndu/nu*100 << " < " << eps_rel*100 << " %\n");
			vel_converged=true;
		}
		else{ 
			if(VerboseConvergence()) _printf0_(setw(50) << left << "   Convergence criterion: norm(du)/norm(u)" << ndu/nu*100 << " > " << eps_rel*100 << " %\n");
			vel_converged=false;
		}
		Vector<IssmDouble>* dup=pug_old->Duplicate(); pug_old->Copy(dup); dup->AYPX(pug,-1.0);
		IssmDouble ndp=dup->Norm(NORM_TWO);   delete dup;
		IssmDouble np =pug_old->Norm(NORM_TWO);
		if (xIsNan<IssmDouble>(ndp) || xIsNan<IssmDouble>(np)) _error_("convergence criterion is NaN!");
		if((ndp/np)<eps_rel){
			if(VerboseConvergence()) _printf0_(setw(50) << left << "   Convergence criterion: norm(dp)/norm(p)" << ndp/np*100 << " < " << eps_rel*100 << " %\n");
			pressure_converged=true;
		}
		else{ 
			if(VerboseConvergence()) _printf0_(setw(50) << left << "   Convergence criterion: norm(dp/)/norm(p)" << ndp/np*100 << " > " << eps_rel*100 << " %\n");
			pressure_converged=false;
		}

		if(pressure_converged && vel_converged) break;
		if(count>=max_nonlinear_iterations){
			_printf0_("   maximum number of nonlinear iterations (" << max_nonlinear_iterations << ") exceeded\n"); 
			break;
		}
	}

	if(VerboseConvergence()) _printf0_("\n   total number of iterations: " << count-1 << "\n");

	delete pug;  
	delete pug_old;  
	delete vel;  
	delete vel_old;  
	delete vel_old_local;  
	delete stressanalysis;
	delete pressureanalysis;
}

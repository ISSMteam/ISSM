/*!\file: solutionsequence_nonlinear.cpp
 * \brief: core of a non-linear solution, using fixed-point method 
 */ 

#include "./solutionsequences.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"

void solutionsequence_newton(FemModel* femmodel){

	/*intermediary: */
	bool   converged;
	int    count,newton;
	IssmDouble kmax;
	Matrix<IssmDouble>* Kff = NULL;
	Matrix<IssmDouble>* Kfs    = NULL;
	Matrix<IssmDouble>* Jff = NULL;
	Vector<IssmDouble>* ug  = NULL;
	Vector<IssmDouble>* old_ug = NULL;
	Vector<IssmDouble>* uf  = NULL;
	Vector<IssmDouble>* old_uf = NULL;
	Vector<IssmDouble>* duf = NULL;
	Vector<IssmDouble>* pf  = NULL;
	Vector<IssmDouble>* pJf    = NULL;
	Vector<IssmDouble>* df  = NULL;
	Vector<IssmDouble>* ys  = NULL;

	/*parameters:*/
	int max_nonlinear_iterations;
	IssmDouble eps_res,eps_rel,eps_abs;

	/*Recover parameters: */
	femmodel->parameters->FindParam(&max_nonlinear_iterations,StressbalanceMaxiterEnum);
	femmodel->parameters->FindParam(&newton,StressbalanceIsnewtonEnum);
	femmodel->parameters->FindParam(&eps_res,StressbalanceRestolEnum);
	femmodel->parameters->FindParam(&eps_rel,StressbalanceReltolEnum);
	femmodel->parameters->FindParam(&eps_abs,StressbalanceAbstolEnum);
	femmodel->UpdateConstraintsx();

	count=0;
	converged=false;

	/*Start non-linear iteration using input velocity: */
	GetSolutionFromInputsx(&ug,femmodel);
	Reducevectorgtofx(&uf,ug,femmodel->nodes,femmodel->parameters);

	//Update once again the solution to make sure that vx and vxold are similar (for next step in transient or steadystate)
	InputUpdateFromConstantx(femmodel,converged,ConvergedEnum);
	InputUpdateFromSolutionx(femmodel,ug);

	for(;;){

		delete old_ug;old_ug=ug;
		delete old_uf;old_uf=uf;

		/*Solver forward model*/
		if(count==0 || newton==2){
			SystemMatricesx(&Kff,&Kfs,&pf,&df,NULL,femmodel);
			CreateNodalConstraintsx(&ys,femmodel->nodes);
			Reduceloadx(pf,Kfs,ys);delete Kfs;
			femmodel->profiler->Start(SOLVER);
			Solverx(&uf,Kff,pf,old_uf,df,femmodel->parameters);delete df; delete Kff; delete pf;
			femmodel->profiler->Stop(SOLVER);
			Mergesolutionfromftogx(&ug,uf,ys,femmodel->nodes,femmodel->parameters);delete ys;
			InputUpdateFromSolutionx(femmodel,ug);
			delete old_ug;old_ug=ug;
			delete old_uf;old_uf=uf;
		}
		uf=old_uf->Duplicate(); old_uf->Copy(uf);

		/*Prepare next iteration using Newton's method*/
		SystemMatricesx(&Kff,&Kfs,&pf,&df,&kmax,femmodel);delete df;
		CreateNodalConstraintsx(&ys,femmodel->nodes);
		Reduceloadx(pf,Kfs,ys);delete Kfs;

		pJf=pf->Duplicate();
		Kff->MatMult(uf,pJf);
		pJf->Scale(-1.0); pJf->AXPY(pf,+1.0);

		CreateJacobianMatrixx(&Jff,femmodel,kmax);
		femmodel->profiler->Start(SOLVER);
		Solverx(&duf,Jff,pJf,NULL,NULL,femmodel->parameters); delete Jff; delete pJf;
		femmodel->profiler->Stop(SOLVER);
		uf->AXPY(duf, 1.0); delete duf;
		Mergesolutionfromftogx(&ug,uf,ys,femmodel->nodes,femmodel->parameters);delete ys;
		InputUpdateFromSolutionx(femmodel,ug);
		count++;

		/*Check convergence*/
		convergence(&converged,Kff,pf,uf,old_uf,eps_res,eps_rel,eps_abs); 
		delete Kff; delete pf;
		if(converged==true) break;
		if(count>=max_nonlinear_iterations){
			_printf0_("   maximum number of Newton iterations (" << max_nonlinear_iterations << ") exceeded\n"); 
			break;
		}
	}

	if(VerboseConvergence()) _printf0_("\n   total number of iterations: " << count << "\n");
	femmodel->results->AddResult(new GenericExternalResult<int>(femmodel->results->Size()+1,StressbalanceConvergenceNumStepsEnum,count));

	/*clean-up*/
	delete uf;
	delete ug;
	delete old_ug;
	delete old_uf;
}

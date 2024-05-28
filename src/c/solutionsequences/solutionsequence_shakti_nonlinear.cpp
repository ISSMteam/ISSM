/*!\file: solutionsequence_nonlinear.cpp
 * \brief: core of a non-linear solution, using fixed-point method 
 */ 

#include "./solutionsequences.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"

void solutionsequence_shakti_nonlinear(FemModel* femmodel){

	/*intermediary: */
	Matrix<IssmDouble>* Kff = NULL;
	Matrix<IssmDouble>* Kfs = NULL;
	Vector<IssmDouble>* ug  = NULL;
	Vector<IssmDouble>* uf  = NULL;
	Vector<IssmDouble>* old_uf = NULL;
	Vector<IssmDouble>* pf  = NULL;
	Vector<IssmDouble>* df  = NULL;
	Vector<IssmDouble>* ys  = NULL;

	/*parameters:*/
	int max_nonlinear_iterations;
	IssmDouble eps_res,eps_rel,eps_abs;

	/*Recover parameters: */
	femmodel->parameters->FindParam(&max_nonlinear_iterations,StressbalanceMaxiterEnum);
	femmodel->parameters->FindParam(&eps_res,StressbalanceRestolEnum);
	femmodel->parameters->FindParam(&eps_rel,StressbalanceReltolEnum);
	femmodel->parameters->FindParam(&eps_abs,StressbalanceAbstolEnum);
	femmodel->UpdateConstraintsx();

	int  count=0;
	bool converged=false;

	/*Start non-linear iteration using input velocity: */
	GetSolutionFromInputsx(&ug,femmodel);
	Reducevectorgtofx(&uf, ug, femmodel->nodes,femmodel->parameters);

	for(;;){
		/*save pointer to old solution*/
		delete old_uf;old_uf=uf;
		delete ug;

		SystemMatricesx(&Kff,&Kfs,&pf,&df,NULL,femmodel);
		CreateNodalConstraintsx(&ys,femmodel->nodes);
		Reduceloadx(pf, Kfs, ys); delete Kfs;
		femmodel->profiler->Start(SOLVER);
		Solverx(&uf, Kff, pf, old_uf, df, femmodel->parameters);
		femmodel->profiler->Stop(SOLVER);
		Mergesolutionfromftogx(&ug, uf,ys,femmodel->nodes,femmodel->parameters);delete ys;

		convergence(&converged,Kff,pf,uf,old_uf,eps_res,eps_rel,eps_abs); delete Kff; delete pf; delete df;
		InputUpdateFromSolutionx(femmodel,ug);

		/*Increase count: */
		count++;
		if(count>=max_nonlinear_iterations){
			_printf0_("   maximum number of nonlinear iterations (" << max_nonlinear_iterations << ") exceeded\n"); 
			converged = true;
		}
		if(converged) break;
	}

	if(VerboseConvergence()) _printf0_("\n   total number of iterations: " << count << "\n");

	/*clean-up*/
	delete uf;
	delete ug;
	delete old_uf;
}

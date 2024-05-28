/*!\file: solutionsequence_nonlinear.cpp
 * \brief: core of a non-linear solution, using fixed-point method 
 */ 

#include "./solutionsequences.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"

void solutionsequence_glads_nonlinear(FemModel* femmodel){

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
	HydrologyGlaDSAnalysis* analysis = new HydrologyGlaDSAnalysis();

	/*Recover parameters (FIXME: from Stress balance for now :( )*/
	femmodel->parameters->FindParam(&max_nonlinear_iterations,StressbalanceMaxiterEnum);
	femmodel->parameters->FindParam(&eps_res,StressbalanceRestolEnum);
	femmodel->parameters->FindParam(&eps_rel,StressbalanceReltolEnum);
	femmodel->parameters->FindParam(&eps_abs,StressbalanceAbstolEnum);
	femmodel->UpdateConstraintsx();

	int  count_out=0;
	bool converged_out=false;
	int  count_in=0;
	bool converged_in=false;

	/*Start non-linear iteration using input velocity: */
	GetSolutionFromInputsx(&ug,femmodel);
	Reducevectorgtofx(&uf, ug, femmodel->nodes,femmodel->parameters);

	while(!converged_out){

		count_in=0;
		converged_in=false;

		while(!converged_in){
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

			convergence(&converged_in,Kff,pf,uf,old_uf,eps_res,eps_rel,eps_abs); delete Kff; delete pf; delete df;
			InputUpdateFromSolutionx(femmodel,ug);

			/*Increase count: */
			count_in++;
			if(count_in>=max_nonlinear_iterations && !converged_in){
				//_printf0_("   maximum number of nonlinear iterations of inner loop (" << max_nonlinear_iterations << ") exceeded\n"); 
				converged_in = true;
			}
		}
		if(VerboseConvergence()) _printf0_(setw(50) << left << "   Inner loop converged in "<<count_in<<" iterations\n");

		if(VerboseConvergence()) _printf0_("   updating sheet thickness and channels cross section\n");
		analysis->UpdateSheetThickness(femmodel);
		analysis->UpdateChannelCrossSection(femmodel);

		/*Converged if inner loop converged in one solution*/
		if(count_in==1) converged_out = true;

		/*Increase count: */
		count_out++;
		if(count_out>=max_nonlinear_iterations && !converged_out){
			_printf0_("   maximum number of nonlinear iterations of outer loop (" << max_nonlinear_iterations << ") exceeded\n"); 
			converged_out = true;
		}
	}

	if(VerboseConvergence()) _printf0_("\n   total number of iterations: " << count_out<<"x"<<count_in<<"\n");

	/*clean-up*/
	delete uf;
	delete ug;
	delete old_uf;
	delete analysis;
}

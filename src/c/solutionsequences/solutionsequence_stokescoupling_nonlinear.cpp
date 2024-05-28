/*!\file: solutionsequence_FScoupling_nonlinear.cpp
 * \brief: core of the coupling between FS and SSAHO
 */ 

#include "./solutionsequences.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"

void solutionsequence_FScoupling_nonlinear(FemModel* femmodel,bool conserve_loads){

	/*intermediary: */
	Matrix<IssmDouble> *Kff_horiz    = NULL;
	Matrix<IssmDouble> *Kfs_horiz    = NULL;
	Vector<IssmDouble> *ug_horiz     = NULL;
	Vector<IssmDouble> *uf_horiz     = NULL;
	Vector<IssmDouble> *old_uf_horiz = NULL;
	Vector<IssmDouble> *pf_horiz     = NULL;
	Vector<IssmDouble> *df_horiz     = NULL;
	Matrix<IssmDouble> *Kff_vert     = NULL;
	Matrix<IssmDouble> *Kfs_vert     = NULL;
	Vector<IssmDouble> *ug_vert      = NULL;
	Vector<IssmDouble> *uf_vert      = NULL;
	Vector<IssmDouble> *pf_vert      = NULL;
	Vector<IssmDouble> *df_vert      = NULL;
	Vector<IssmDouble> *ys           = NULL;
	bool converged;
	int  count;

	/*parameters:*/
	int  min_mechanical_constraints;
	int  max_nonlinear_iterations;
	IssmDouble eps_res,eps_rel,eps_abs;

	/*Recover parameters: */
	femmodel->parameters->FindParam(&min_mechanical_constraints,StressbalanceRiftPenaltyThresholdEnum);
	femmodel->parameters->FindParam(&max_nonlinear_iterations,StressbalanceMaxiterEnum);
	femmodel->parameters->FindParam(&eps_res,StressbalanceRestolEnum);
	femmodel->parameters->FindParam(&eps_rel,StressbalanceReltolEnum);
	femmodel->parameters->FindParam(&eps_abs,StressbalanceAbstolEnum);
	femmodel->UpdateConstraintsx();

	count=1;
	converged=false;

	/*First get ug_horiz:*/
	femmodel->SetCurrentConfiguration(StressbalanceAnalysisEnum);
	GetSolutionFromInputsx(&ug_horiz,femmodel);
	Reducevectorgtofx(&uf_horiz, ug_horiz, femmodel->nodes,femmodel->parameters);

	for(;;){

		/*First stressbalance horiz:*/
		femmodel->SetCurrentConfiguration(StressbalanceAnalysisEnum);

		//Update once again the solution to make sure that vx and vxold are similar (for next step in transient or steadystate)
		InputUpdateFromSolutionx(femmodel,ug_horiz);
		delete ug_horiz;

		//save pointer to old velocity
		delete old_uf_horiz; old_uf_horiz=uf_horiz;

		/*solve: */
		SystemMatricesx(&Kff_horiz, &Kfs_horiz, &pf_horiz, &df_horiz, NULL,femmodel);
		CreateNodalConstraintsx(&ys,femmodel->nodes);
		Reduceloadx(pf_horiz, Kfs_horiz, ys); delete Kfs_horiz;
		femmodel->profiler->Start(SOLVER);
		Solverx(&uf_horiz, Kff_horiz, pf_horiz, old_uf_horiz, df_horiz,femmodel->parameters);
		femmodel->profiler->Stop(SOLVER);
		Mergesolutionfromftogx(&ug_horiz, uf_horiz,ys,femmodel->nodes,femmodel->parameters); delete ys;
		InputUpdateFromSolutionx(femmodel,ug_horiz);

		convergence(&converged,Kff_horiz,pf_horiz,uf_horiz,old_uf_horiz,eps_res,eps_rel,eps_abs); delete Kff_horiz; delete pf_horiz; delete df_horiz;

		/*Second compute vertical velocity: */
		femmodel->SetCurrentConfiguration(StressbalanceVerticalAnalysisEnum);

		/*solve: */
		SystemMatricesx(&Kff_vert, &Kfs_vert, &pf_vert, &df_vert,NULL,femmodel);
		CreateNodalConstraintsx(&ys,femmodel->nodes);
		Reduceloadx(pf_vert, Kfs_vert, ys); delete Kfs_vert;
		femmodel->profiler->Start(SOLVER);
		Solverx(&uf_vert, Kff_vert, pf_vert, NULL, df_vert,femmodel->parameters); delete Kff_vert; delete pf_vert; delete df_vert;
		femmodel->profiler->Stop(SOLVER);
		Mergesolutionfromftogx(&ug_vert, uf_vert,ys,femmodel->nodes,femmodel->parameters);
		delete uf_vert; 
		delete ys; 
		InputUpdateFromSolutionx(femmodel,ug_vert);
		delete ug_vert;

		/*Increase count: */
		count++;
		if(converged==true)break;
		if(count>=max_nonlinear_iterations){
			_printf0_("   maximum number of iterations (" << max_nonlinear_iterations << ") exceeded\n"); 
			break;
		}
	}

	/*clean-up*/
	delete old_uf_horiz;
	delete uf_horiz;
	delete ug_horiz;
}

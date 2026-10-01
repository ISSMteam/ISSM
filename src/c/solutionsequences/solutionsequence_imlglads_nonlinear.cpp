/*!\file: solutionsequence_nonlinear.cpp
 * \brief: core of a non-linear solution, using fixed-point method 
 */ 
#include "../cores/cores.h"
#include "./solutionsequences.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"
#include <float.h> /* defines DBL_EPSILON*/
#define h_reltol      1e-3
#define h_abstol      1e-1
/*local routine to check qr/lh convergence:{{{*/
bool lakelhconvergence(Vector<IssmDouble>* lh_new,Vector<IssmDouble>* lh_old);


void solutionsequence_imlglads_nonlinear(FemModel* femmodel){/*{{{*/

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
	bool islakes;
	HydrologyIMLGlaDSAnalysis* analysis = new HydrologyIMLGlaDSAnalysis();

	/*Recover parameters*/
	femmodel->parameters->FindParam(&max_nonlinear_iterations,HydrologyMaxiterEnum);
	femmodel->parameters->FindParam(&eps_res,HydrologyRestolEnum);
	femmodel->parameters->FindParam(&eps_rel,HydrologyReltolEnum);
	femmodel->parameters->FindParam(&eps_abs,HydrologyAbstolEnum);
	femmodel->parameters->FindParam(&islakes,HydrologyLakeFlagEnum);
	femmodel->UpdateConstraintsx();

	/*lh convergence criterion*/
	Vector<IssmDouble>* lh_old = NULL;
	if(islakes){
		GetVectorFromInputsx(&lh_old,femmodel,HydrologyLakeHeightOldEnum,VertexSIdEnum);
	}

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
		if(islakes){
			/*Update lake height based on Qr i) from previous timestep, then ii) from each iteration of this outer loop*/
			if(VerboseConvergence()) _printf0_("   updating lake height\n");
			analysis->UpdateLakeDepth(femmodel);
			/*reset boundary conditions*/
			if(VerboseConvergence()) _printf0_("   updating phi at the lake outlet\n");
			ResetBoundaryConditions(femmodel,HydrologyIMLGlaDSAnalysisEnum);
		}

		/*begin inner loop*/
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

		if(VerboseConvergence()) _printf0_("   updating sheet thickness\n");
		analysis->UpdateSheetThickness(femmodel);
		if(!islakes){
			if(VerboseConvergence()) _printf0_("   updating channel cross section\n");
			analysis->UpdateChannelCrossSection(femmodel);
			/*Converged if inner loop converged in one solution*/
			if(count_in==1) converged_out = true;
		}
		else {
			converged_out = true; /*initially set to true, change if any lakes do not converge*/
			if(VerboseConvergence()) _printf0_("   updating channel cross section and lake outlet flux\n");
			analysis->UpdateChannelCrossSection(femmodel);
			analysis->UpdateLakeOutletDischarge(femmodel);
			
			/*check convergence of lake height*/
			Vector<IssmDouble>* lh_new = NULL;
			GetVectorFromInputsx(&lh_new,femmodel,HydrologyLakeHeightEnum,VertexSIdEnum);
			if(!lakelhconvergence(lh_new,lh_old)){
				converged_out = false;
			}
			/*Update lh_old to lh_new for next iteration*/
			delete lh_old;
			lh_old = lh_new;
			lh_new = NULL;
		}
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
	delete lh_old;
	delete analysis;
}/*}}}*/

bool lakelhconvergence(Vector<IssmDouble>* lh_new, Vector<IssmDouble>* lh_old){/*{{{*/
	
	/*Compute difference*/
	Vector<IssmDouble>* diff = lh_new->Duplicate();
	lh_new->Copy(diff);
	diff->AXPY(lh_old,-1.0);
	
	/*Compute norms*/
	IssmDouble diff_inf = diff->Norm(NORM_INF); // max pointwise change
	IssmDouble nlh		= lh_new->Norm(NORM_TWO); // for optional relative scaling

	bool converged = false;

	if(!xIsNan<IssmDouble>(h_abstol)){
		/* Absolute criterion: max change anywhere is small*/
		converged = (diff_inf < h_abstol);
	}

	if(!converged && !xIsNan<IssmDouble>(h_reltol)){
		/* Relative criterion: max change anywhere is small relative to overall magnitude*/
		IssmDouble rel = diff_inf / (nlh + DBL_EPSILON);
		converged = (rel < h_reltol);
	}

	if(VerboseConvergence()){
		_printf0_(setw(50) << left
        << "              max abs lh diff: " << diff_inf << " m (tol: " << h_abstol << " m)"
        << "  rel lh diff: " << diff_inf / (nlh + DBL_EPSILON) * 100 << " % (tol: " << h_reltol * 100 << " %)\n");
	}

    // Edge case: both essentially zero
    if(nlh < 1e-10 && lh_old->Norm(NORM_TWO) < 1e-10) converged = true;

	delete diff;
	return converged;
}/*}}}*/




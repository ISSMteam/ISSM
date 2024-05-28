/*!\file: solutionsequence_nonlinear.cpp
 * \brief: core of a non-linear solution, using fixed-point method 
 */ 

#include "./solutionsequences.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"

void solutionsequence_nonlinear(FemModel* femmodel,bool conserve_loads){

	/*intermediary: */
	Matrix<IssmDouble>* Kff = NULL;
	Matrix<IssmDouble>* Kfs = NULL;
	Vector<IssmDouble>* ug  = NULL;
	Vector<IssmDouble>* uf  = NULL;
	Vector<IssmDouble>* old_uf = NULL;
	Vector<IssmDouble>* pf  = NULL;
	Vector<IssmDouble>* df  = NULL;
	Vector<IssmDouble>* ys  = NULL;

	int constraints_converged;
	int num_unstable_constraints;

	/*parameters:*/
	int min_mechanical_constraints;
	int max_nonlinear_iterations;
	int configuration_type;
	IssmDouble eps_res,eps_rel,eps_abs;

	/*Recover parameters: */
	femmodel->parameters->FindParam(&min_mechanical_constraints,StressbalanceRiftPenaltyThresholdEnum);
	femmodel->parameters->FindParam(&max_nonlinear_iterations,StressbalanceMaxiterEnum);
	femmodel->parameters->FindParam(&eps_res,StressbalanceRestolEnum);
	femmodel->parameters->FindParam(&eps_rel,StressbalanceReltolEnum);
	femmodel->parameters->FindParam(&eps_abs,StressbalanceAbstolEnum);
	femmodel->parameters->FindParam(&configuration_type,ConfigurationTypeEnum);
	femmodel->UpdateConstraintsx();

	/*Were loads requested as output? : */
	Loads* savedloads=NULL;
	if(conserve_loads){
		savedloads = static_cast<Loads*>(femmodel->loads->Copy());
	}

	int  count=0;
	bool converged=false;

	/*Start non-linear iteration using input velocity: */
	GetSolutionFromInputsx(&ug,femmodel);
	Reducevectorgtofx(&uf, ug, femmodel->nodes,femmodel->parameters);

	/*Update once again the solution to make sure that vx and vxold are similar (for next step in transient or steadystate)*/
	InputUpdateFromConstantx(femmodel,converged,ConvergedEnum);
	InputUpdateFromSolutionx(femmodel,ug);

	/*allocate the matrices once and reuse them per iteration*/
	if(femmodel->loads->numrifts == 0){
		AllocateSystemMatricesx(&Kff,&Kfs,&df,&pf,femmodel);
	}

	for(;;){

		//save pointer to old velocity
		delete old_uf;old_uf=uf;
		delete ug;

		if(femmodel->loads->numrifts){
			AllocateSystemMatricesx(&Kff,&Kfs,&df,&pf,femmodel);
		}
		SystemMatricesx(&Kff,&Kfs,&pf,&df,NULL,femmodel, true);
		CreateNodalConstraintsx(&ys,femmodel->nodes);
		Reduceloadx(pf, Kfs, ys);
		femmodel->profiler->Start(SOLVER);
		Solverx(&uf, Kff, pf, old_uf, df, femmodel->parameters);
		femmodel->profiler->Stop(SOLVER);

		Mergesolutionfromftogx(&ug, uf,ys,femmodel->nodes,femmodel->parameters);delete ys;

		convergence(&converged,Kff,pf,uf,old_uf,eps_res,eps_rel,eps_abs);
		InputUpdateFromConstantx(femmodel,converged,ConvergedEnum);
		InputUpdateFromSolutionx(femmodel,ug);

		/*Clean up if rifts*/
		if(femmodel->loads->numrifts){
			delete Kfs; delete Kff; delete pf; delete df;
		}

		ConstraintsStatex(&constraints_converged,&num_unstable_constraints,femmodel);
		if(VerboseConvergence()) _printf0_("   number of unstable constraints: " << num_unstable_constraints << "\n");

		//rift convergence
		if (!constraints_converged) {
			if (converged){
				if (num_unstable_constraints <= min_mechanical_constraints) converged=true;
				else converged=false;
			}
		}

		/*Increase count: */
		count++;
		if(converged==true){
			femmodel->results->AddResult(new GenericExternalResult<int>(femmodel->results->Size()+1,StressbalanceConvergenceNumStepsEnum,count));
			break;
		}
		if(count>=max_nonlinear_iterations){
			_printf0_("   maximum number of nonlinear iterations (" << max_nonlinear_iterations << ") exceeded\n"); 
			converged=true;
			femmodel->results->AddResult(new GenericExternalResult<int>(femmodel->results->Size()+1,StressbalanceConvergenceNumStepsEnum,max_nonlinear_iterations));
			InputUpdateFromConstantx(femmodel,converged,ConvergedEnum);
			InputUpdateFromSolutionx(femmodel,ug);		
			break;
		}

		/*Set the matrix entries to zero if we do an other iteration*/
		if(femmodel->loads->numrifts==0){
			Kff->SetZero();
			Kfs->SetZero();
			df->Set(0);
			pf->Set(0);
		}
	}

	/*delete matrices after the iteration loop*/
	if(femmodel->loads->numrifts==0){
		delete Kff; delete pf; delete df; delete Kfs;
	}

	if(VerboseConvergence()) _printf0_("\n   total number of iterations: " << count << "\n");

	/*clean-up*/
	if(conserve_loads){
		delete femmodel->loads;
		int index=femmodel->AnalysisIndex(configuration_type);
		femmodel->loads_list[index]=savedloads;
		femmodel->loads=savedloads;
	}
	delete uf;
	delete ug;
	delete old_uf;
}

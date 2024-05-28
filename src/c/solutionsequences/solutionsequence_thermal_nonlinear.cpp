/*!\file: solutionsequence_thermal_nonlinear.cpp
 * \brief: core of the thermal solution 
 */ 

#include "./solutionsequences.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"

void solutionsequence_thermal_nonlinear(FemModel* femmodel){

	/*solution : */
	Vector<IssmDouble>* tg=NULL; 
	Vector<IssmDouble>* tf=NULL; 
	Vector<IssmDouble>* tf_old=NULL; 
	Vector<IssmDouble>* ys=NULL; 
	IssmDouble melting_offset;

	/*intermediary: */
	Matrix<IssmDouble>* Kff=NULL;
	Matrix<IssmDouble>* Kfs=NULL;
	Vector<IssmDouble>* pf=NULL;
	Vector<IssmDouble>* df=NULL;

	bool converged;
	bool isenthalpy, isdynamicbasalspc;
	int constraints_converged;
	int num_unstable_constraints;
	int count;
	int thermal_penalty_threshold;
	int thermal_maxiter;

	/*parameters:*/
	IssmDouble eps_rel;

	/*Recover parameters: */
	femmodel->parameters->FindParam(&isenthalpy,ThermalIsenthalpyEnum);
	femmodel->parameters->FindParam(&thermal_maxiter,ThermalMaxiterEnum);

	converged=false;
	InputUpdateFromConstantx(femmodel,converged,ConvergedEnum);

	if(isenthalpy){
		femmodel->parameters->FindParam(&isdynamicbasalspc,ThermalIsdynamicbasalspcEnum);
		femmodel->parameters->FindParam(&eps_rel,ThermalReltolEnum);
		femmodel->UpdateConstraintsx();

		//Update the solution to make sure that tf and tf_old are similar (for next step in transient or steadystate)
		GetSolutionFromInputsx(&tg,femmodel);
		Reducevectorgtofx(&tf, tg, femmodel->nodes,femmodel->parameters);
		InputUpdateFromSolutionx(femmodel,tg);
	}
	else{
		femmodel->parameters->FindParam(&thermal_penalty_threshold,ThermalPenaltyThresholdEnum);
		ResetZigzagCounterx(femmodel);
		femmodel->UpdateConstraintsx();
	}

	count=1;

	for(;;){
		delete tf_old;tf_old=tf;

		if(isenthalpy){ 
			SystemMatricesx(&Kff,&Kfs,&pf,&df,NULL,femmodel);
			/*Update old solution, such that sizes of tf_old and tf are comparable*/
			if(isdynamicbasalspc){
				delete tf_old;
				Reducevectorgtofx(&tf_old, tg, femmodel->nodes,femmodel->parameters);
			}
		}
		else SystemMatricesx(&Kff, &Kfs, &pf,&df, &melting_offset,femmodel);
		delete tg;
		CreateNodalConstraintsx(&ys,femmodel->nodes);
		Reduceloadx(pf, Kfs, ys); delete Kfs;
		femmodel->profiler->Start(SOLVER);
		Solverx(&tf, Kff, pf, tf_old, df, femmodel->parameters);
		femmodel->profiler->Stop(SOLVER);
		Mergesolutionfromftogx(&tg, tf,ys,femmodel->nodes,femmodel->parameters); delete ys;
		if(isenthalpy){ 
			convergence(&converged,Kff,pf,tf,tf_old,0.05,eps_rel,NAN); 
			InputUpdateFromConstantx(femmodel,converged,ConvergedEnum);
		}
		delete Kff; delete pf; delete df;
		InputUpdateFromSolutionx(femmodel,tg);
		ConstraintsStatex(&constraints_converged,&num_unstable_constraints,femmodel);
		if(VerboseConvergence()) _printf0_("   number of unstable constraints: " << num_unstable_constraints << "\n");

		if(isenthalpy){ // enthalpy method
			IssmDouble dt;
			femmodel->parameters->FindParam(&dt,TimesteppingTimeStepEnum);

			count++;
			bool max_iteration_state=false;
			if(count>=thermal_maxiter){
				_printf0_("   maximum number of nonlinear iterations (" << thermal_maxiter << ") exceeded\n"); 
				converged=true;
				InputUpdateFromConstantx(femmodel,converged,ConvergedEnum);
				InputUpdateFromSolutionx(femmodel,tg);		
				max_iteration_state=true;
			}
			if(converged==true){
				break;
			}
			else if(dt==0.){
				EnthalpyAnalysis::ComputeBasalMeltingrate(femmodel);
				EnthalpyAnalysis::UpdateBasalConstraints(femmodel);
			}
		}
		else{ // dry ice method
			if(!converged){
				if(num_unstable_constraints<=thermal_penalty_threshold) converged=true;
				if(count>=thermal_maxiter){
					converged=true;
					_printf0_("   maximum number of iterations (" << thermal_maxiter << ") exceeded\n"); 
				}
			}
			count++;
			InputUpdateFromConstantx(femmodel,converged,ConvergedEnum);
			if(converged)break;
		}
	}

	if(isenthalpy){
		if(VerboseConvergence()) _printf0_("\n   total number of iterations: " << count-1 << "\n");
	}
	else{
		InputUpdateFromSolutionx(femmodel,tg);
		femmodel->parameters->SetParam(melting_offset,MeltingOffsetEnum);
	}

	/*Free resources: */
	delete tg;
	delete tf;
	delete tf_old;
}

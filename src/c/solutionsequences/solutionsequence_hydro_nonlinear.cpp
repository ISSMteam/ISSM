/*!\file: solutionsequence_hydro_nonlinear.cpp
 * \brief: core of the hydro solution
 */

#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"
/*FIXME, dirty hack to get the solutionsequence linear needed to compute the slopes*/
#include "../solutionsequences/solutionsequences.h"

void solutionsequence_hydro_nonlinear(FemModel* femmodel, bool* pconv_fail){
	/*solution : */
	Vector<IssmDouble>* ug_sed=NULL;
	Vector<IssmDouble>* uf_sed=NULL;
	Vector<IssmDouble>* uf_sed_sub_iter=NULL;
	Vector<IssmDouble>* ug_sed_main_iter=NULL;
	Vector<IssmDouble>* ug_sed_init=NULL;

	Vector<IssmDouble>* ug_epl=NULL;
	Vector<IssmDouble>* uf_epl=NULL;
	Vector<IssmDouble>* uf_epl_sub_iter=NULL;
	Vector<IssmDouble>* ug_epl_main_iter=NULL;
	Vector<IssmDouble>* ug_epl_init=NULL;

	Vector<IssmDouble>* ys=NULL;
	Vector<IssmDouble>* dug=NULL;

	Matrix<IssmDouble>* Kff=NULL;
	Matrix<IssmDouble>* Kfs=NULL;

	Vector<IssmDouble>* pf=NULL;
	Vector<IssmDouble>* df=NULL;

	HydrologyDCInefficientAnalysis* inefanalysis = NULL;
	HydrologyDCEfficientAnalysis* effanalysis = NULL;

	bool       sedconverged,eplconverged,hydroconverged;
	bool       isefficientlayer;
	bool       sliceadapt;
	int        constraints_converged;
	int        num_unstable_constraints;
	int        sedcount,eplcount,hydrocount;
	int        hydro_maxiter;
	int        epl_fsize,epl_sub_fsize,epl_main_fsize;
	IssmDouble sediment_kmax;
	IssmDouble eps_res,eps_rel,eps_abs;
	IssmDouble ndu_sed,nu_sed;
	IssmDouble ndu_epl,nu_epl;
	IssmDouble ThickCount,L2Count;
	/*Recover parameters: */
	femmodel->SetCurrentConfiguration(HydrologyDCInefficientAnalysisEnum);
	femmodel->parameters->FindParam(&isefficientlayer,HydrologydcIsefficientlayerEnum);
	femmodel->parameters->FindParam(&sliceadapt,HydrologyStepAdaptEnum);
	femmodel->parameters->FindParam(&hydro_maxiter,HydrologydcMaxIterEnum);
	femmodel->parameters->FindParam(&eps_res,StressbalanceRestolEnum);
	femmodel->parameters->FindParam(&eps_rel,HydrologydcRelTolEnum);
	femmodel->parameters->FindParam(&eps_abs,StressbalanceAbstolEnum);
	hydrocount=1;
	hydroconverged=false;
	/*We don't need the outer loop if only one layer is used*/
	if(!isefficientlayer) hydroconverged=true;

	/*{{{*//*Retrieve inputs as the initial state for the non linear iteration*/
	GetBasalSolutionFromInputsx(&ug_sed,femmodel);
	/*Initialize the IDS element mask to exclude frozen nodes*/
	inefanalysis = new HydrologyDCInefficientAnalysis();
	inefanalysis->ElementizeIdsMask(femmodel);

	Reducevectorgtofx(&uf_sed, ug_sed, femmodel->nodes,femmodel->parameters);
	ug_sed_init=ug_sed->Duplicate();
	ug_sed->Copy(ug_sed_init);

	if(isefficientlayer) {
		effanalysis = new HydrologyDCEfficientAnalysis();
		femmodel->SetCurrentConfiguration(HydrologyDCEfficientAnalysisEnum);
		GetBasalSolutionFromInputsx(&ug_epl,femmodel);
		inefanalysis->ElementizeEplMask(femmodel);
		effanalysis->InitZigZagCounter(femmodel);
		ug_epl_init=ug_epl->Duplicate();
		ug_epl->Copy(ug_epl_init);
	}
	/*}}}*/
	/*The real computation starts here, outermost loop is on the two layer system*/
	for(;;){
		sedcount=1;
		eplcount=1;

		/*If there is two layers we need an outer loop value to compute convergence*/
		if(isefficientlayer){
			ug_sed_main_iter=ug_sed->Duplicate();
			ug_sed->Copy(ug_sed_main_iter);
			ug_epl_main_iter=ug_epl->Duplicate();
			ug_epl->Copy(ug_epl_main_iter);
		}
		/*Loop on sediment layer to deal with transfer and head value*/
		femmodel->SetCurrentConfiguration(HydrologyDCInefficientAnalysisEnum);
		ResetZigzagCounterx(femmodel);
		InputUpdateFromConstantx(femmodel,false,ConvergedEnum);
		femmodel->UpdateConstraintsx();

		/*Reset constraint on the ZigZag Lock*/
		ResetConstraintsx(femmodel);

		/*{{{*//*Treating the sediment*/
		femmodel->profiler->Start(SEDLOOP);
		for(;;){
			sedconverged=false;
			uf_sed_sub_iter=uf_sed->Duplicate();_assert_(uf_sed_sub_iter);
			uf_sed->Copy(uf_sed_sub_iter);
			/*{{{*//*Loop on the sediment layer to deal with the penalization*/
			for(;;){
				/*{{{*//*Core of the computation*/
				if(VerboseSolution()) _printf0_("Building Sediment Matrix...\n");

				femmodel->profiler->Start(SEDMatrix);
				SystemMatricesx(&Kff,&Kfs,&pf,&df,&sediment_kmax,femmodel);
				CreateNodalConstraintsx(&ys,femmodel->nodes);
				Reduceloadx(pf,Kfs,ys); delete Kfs;
				delete uf_sed;
				femmodel->profiler->Stop(SEDMatrix);

				femmodel->profiler->Start(SOLVER);
				Solverx(&uf_sed,Kff,pf,uf_sed_sub_iter,df,femmodel->parameters);
				delete df;
				femmodel->profiler->Stop(SOLVER);
				delete ug_sed;
				femmodel->profiler->Start(SEDUpdate);
				Mergesolutionfromftogx(&ug_sed,uf_sed,ys,femmodel->nodes,femmodel->parameters); delete ys;
				InputUpdateFromSolutionx(femmodel,ug_sed);
				ConstraintsStatex(&constraints_converged,&num_unstable_constraints,femmodel);
				femmodel->profiler->Stop(SEDUpdate);
				/*}}}*/
				if (!sedconverged){
					/*First check that all the penalizations are applied*/
					if(VerboseConvergence()) _printf0_("   # Sediment unstable constraints = " << num_unstable_constraints << "\n");
					if(num_unstable_constraints==0) {
						sedconverged = true;
					}
					else{//clean up
						delete Kff;
						delete pf;
					}
					if (sedcount>=hydro_maxiter){
						delete ug_sed;delete uf_sed;delete inefanalysis; delete ug_sed_main_iter;
						if(isefficientlayer)delete ug_epl;delete effanalysis; delete ug_epl_main_iter;
						_error_("   maximum number of Sediment iterations (" << hydro_maxiter << ") exceeded");

					}
				}
				/*Add an iteration and get out of the loop if the penalisation is converged*/
				sedcount++;
				if(sedconverged)break;
			}
			/*}}}*//*End of the sediment penalization loop*/
			sedconverged=false;
			/*Checking convergence on the value of the sediment head*/
			convergence(&sedconverged,Kff,pf,uf_sed,uf_sed_sub_iter,eps_res,eps_rel,eps_abs);
			delete Kff; delete pf;delete uf_sed_sub_iter;
			if(sedconverged){
				femmodel->parameters->SetParam(sediment_kmax,HydrologySedimentKmaxEnum);
				InputUpdateFromConstantx(femmodel,sedconverged,ConvergedEnum);
				InputUpdateFromSolutionx(femmodel,ug_sed);
				InputUpdateFromConstantx(femmodel,sediment_kmax,HydrologySedimentKmaxEnum);
				break;
			}
		}
		femmodel->profiler->Stop(SEDLOOP);
		/*}}}*//*End of the global sediment loop*/
		/*{{{*//*Now dealing with the EPL in the same way*/
		femmodel->profiler->Start(EPLLOOP);
		if(isefficientlayer){
			femmodel->SetCurrentConfiguration(HydrologyDCEfficientAnalysisEnum);
			/*updating mask*/
			if(VerboseSolution()) _printf0_("==updating mask...\n");
			femmodel->HydrologyEPLupdateDomainx(&ThickCount);
			ResetZigzagCounterx(femmodel);
			InputUpdateFromConstantx(femmodel,false,ConvergedEnum);

			for(;;){
				eplconverged=false;
				/*{{{*//*Retrieve the EPL head slopes and compute EPL Thickness*/
				if(VerboseSolution()) _printf0_("computing EPL Head slope...\n");
				femmodel->profiler->Start(EPLMasking);
				femmodel->SetCurrentConfiguration(L2ProjectionEPLAnalysisEnum);
				femmodel->UpdateConstraintsL2ProjectionEPLx(&L2Count);
				femmodel->parameters->SetParam(EplHeadSlopeXEnum,InputToL2ProjectEnum);
				solutionsequence_linear(femmodel);
				femmodel->parameters->SetParam(EplHeadSlopeYEnum,InputToL2ProjectEnum);
				solutionsequence_linear(femmodel);

				femmodel->SetCurrentConfiguration(HydrologyDCEfficientAnalysisEnum);
				effanalysis->ComputeEPLThickness(femmodel);
				//updating mask after the computation of the epl thickness (Allow to close too thin EPL)
				femmodel->HydrologyEPLupdateDomainx(&ThickCount);
				/*}}}*/
				femmodel->profiler->Stop(EPLMasking);
				Reducevectorgtofx(&uf_epl, ug_epl, femmodel->nodes,femmodel->parameters);
				if(VerboseSolution()) _printf0_("Building EPL Matrix...\n");
				uf_epl_sub_iter=uf_epl->Duplicate();_assert_(uf_epl_sub_iter);
				uf_epl->Copy(uf_epl_sub_iter);
				uf_epl->GetSize(&epl_sub_fsize);
				femmodel->profiler->Start(EPLMatrices);
				SystemMatricesx(&Kff,&Kfs,&pf,&df,NULL,femmodel);
				CreateNodalConstraintsx(&ys,femmodel->nodes);
				Reduceloadx(pf,Kfs,ys);
				delete Kfs;delete uf_epl;
				femmodel->profiler->Stop(EPLMatrices);
				femmodel->profiler->Start(SOLVER);
				Solverx(&uf_epl,Kff,pf,uf_epl_sub_iter,df,femmodel->parameters);
				femmodel->profiler->Stop(SOLVER);
				delete df;delete ug_epl;
				femmodel->profiler->Start(EPLUpdate);
				Mergesolutionfromftogx(&ug_epl,uf_epl,ys,femmodel->nodes,femmodel->parameters); delete ys;
				InputUpdateFromSolutionx(femmodel,ug_epl);
				ConstraintsStatex(&constraints_converged,&num_unstable_constraints,femmodel);
				femmodel->profiler->Stop(EPLUpdate);
				uf_epl->GetSize(&epl_fsize);
				if(epl_fsize-epl_sub_fsize==0){
					convergence(&eplconverged,Kff,pf,uf_epl,uf_epl_sub_iter,eps_res,eps_rel,eps_abs);
					delete Kff; delete pf;
				}
				else{
					delete Kff; delete pf;
				}

				if (eplcount>=hydro_maxiter*9/10 && sliceadapt && !eplconverged) {
					if(VerboseSolution()) _printf0_("epl did not converged after "<<eplconverged<<" iteration, we refine the steping\n");
					*pconv_fail = true;
					InputUpdateFromSolutionx(femmodel,ug_epl_init);
					delete ug_epl_init;
					femmodel->SetCurrentConfiguration(HydrologyDCInefficientAnalysisEnum);
					InputUpdateFromSolutionx(femmodel,ug_sed_init);
					delete ug_sed_init;
					break;
				}
				else if (eplcount>=hydro_maxiter){
					delete ug_sed;delete uf_sed;delete inefanalysis;delete ug_sed_main_iter;
					delete ug_epl;delete uf_epl;delete effanalysis;delete ug_epl_main_iter;
					_error_("   maximum number of EPL iterations (" << hydro_maxiter << ") exceeded");

				}
				eplcount++;
				delete uf_epl_sub_iter; delete uf_epl;
				if(eplconverged){
					if(VerboseSolution()) _printf0_("eplconverged...\n");
					effanalysis->ResetCounter(femmodel);
					break;
				}
			}
		}
		femmodel->profiler->Stop(EPLLOOP);
		if(*pconv_fail){
			break;
		}
		/*}}}*//*End of the global EPL loop*/
		/*{{{*//*Now dealing with the convergence of the whole system*/
		if(!hydroconverged){
			//compute norm(du)/norm(u)
			dug=ug_sed_main_iter->Duplicate(); _assert_(dug);
			ug_sed_main_iter->Copy(dug);
			dug->AYPX(ug_sed,-1.0);
			ndu_sed=dug->Norm(NORM_TWO);
			delete dug;
			nu_sed=ug_sed_main_iter->Norm(NORM_TWO);
			delete ug_sed_main_iter;
			if (xIsNan<IssmDouble>(ndu_sed) || xIsNan<IssmDouble>(nu_sed)) _error_("Sed convergence criterion is NaN!");
			if (ndu_sed==0.0 && nu_sed==0.0) nu_sed=1.0e-6; /*Hacking the case where the Sediment is used but empty*/
			dug=ug_epl_main_iter->Duplicate();_assert_(dug);
			ug_epl_main_iter->Copy(dug);
			dug->AYPX(ug_epl,-1.0);
			ndu_epl=dug->Norm(NORM_TWO);
			delete dug;
			nu_epl=ug_epl_main_iter->Norm(NORM_TWO);
			delete ug_epl_main_iter;
			if (xIsNan<IssmDouble>(ndu_epl) || xIsNan<IssmDouble>(nu_epl)) _error_("EPL convergence criterion is NaN!");
			if (ndu_epl==0.0 && nu_epl==0.0) nu_epl=1.0e-6; /*Hacking the case where the EPL is used but empty*/
			if (!xIsNan<IssmDouble>(eps_rel)){
				if ((ndu_epl/nu_epl)<eps_rel && (ndu_sed/nu_sed)<(eps_rel)){
					if (VerboseConvergence()) _printf0_(setw(50) << left << "   Converged after, " << hydrocount << " iterations \n");
					hydroconverged=true;
				}
				else{
					if(VerboseConvergence()) _printf0_(setw(50) << left << "   for iteration:" << hydrocount << " \n");
					if(VerboseConvergence()) _printf0_(setw(50) << left << "   Sediment Convergence criterion:" << ndu_sed/nu_sed*100 << "%, aiming lower than " << eps_rel*100 << " %\n");
					if(VerboseConvergence()) _printf0_(setw(50) << left << "   EPL Convergence criterion:" << ndu_epl/nu_epl*100 << "%, aiming lower than " << eps_rel*100 << " %\n");
					hydroconverged=false;
				}
			}
			else _printf0_(setw(50) << left << "   Convergence criterion:" << ndu_sed/nu_sed*100 << " %\n");

			if (hydrocount>=hydro_maxiter*9/10  && sliceadapt && !hydroconverged) {
				if(VerboseSolution()) _printf0_("hydrology main loop  did not converged after "<<hydrocount<<" iteration, we refine the steping\n");
				*pconv_fail = true;
				InputUpdateFromSolutionx(femmodel,ug_epl_init);
				delete ug_epl_init;
				femmodel->SetCurrentConfiguration(HydrologyDCInefficientAnalysisEnum);
				InputUpdateFromSolutionx(femmodel,ug_sed_init);
				delete ug_sed_init;
				break;
			}
			if (hydrocount>=hydro_maxiter){
				_error_("   maximum number for hydrological global iterations (" << hydro_maxiter << ") exceeded");
				delete ug_sed;delete uf_sed;delete effanalysis;
				delete ug_epl;	delete uf_epl;
			}
		}
		hydrocount++;
		if(hydroconverged)break;
	}
	/*}}}*/
	/*To deal with adaptative stepping we only save results if we are actually converged*/
	if(hydroconverged){
		if(isefficientlayer){
			delete ug_epl_init;
		}
		delete ug_sed_init;
	}
	/*Free resources: */
	delete ug_epl;delete ug_sed;
	delete uf_sed;
	delete inefanalysis;	delete effanalysis;
}

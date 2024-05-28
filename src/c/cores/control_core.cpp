/*!\file: control_core.cpp
 * \brief: core of the control solution 
 */ 

#include "./cores.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"
#include "../solutionsequences/solutionsequences.h"

/*Local prototypes*/
/*{{{*/
IssmDouble FormFunction(IssmDouble* X,void* usr);
IssmDouble FormFunctionGradient(IssmDouble** pG,IssmDouble* X,void* usr);
typedef struct {
	FemModel* femmodel;
	int       nsize;
} AppCtx;
/*}}}*/

void control_core(FemModel* femmodel){/*{{{*/

	/*parameters: */
	int         num_controls,nsize,nsteps;
	int         solution_type;
	bool        isFS,dakota_analysis;
	int        *control_type  = NULL;
	int        *maxiter       = NULL;
	IssmDouble *cm_jump       = NULL;
	IssmDouble *J             = NULL;

	/*Solution and Adjoint core pointer*/
	void (*solutioncore)(FemModel*) = NULL;
	void (*adjointcore)(FemModel*)  = NULL;

	/*Recover parameters used throughout the solution*/
	femmodel->parameters->FindParam(&num_controls,InversionNumControlParametersEnum);
	femmodel->parameters->FindParam(&control_type,NULL,InversionControlParametersEnum);
	femmodel->parameters->FindParam(&nsteps,InversionNstepsEnum);
	femmodel->parameters->FindParam(&maxiter,NULL,InversionMaxiterPerStepEnum);
	femmodel->parameters->FindParam(&cm_jump,NULL,InversionStepThresholdEnum);
	femmodel->parameters->FindParam(&solution_type,SolutionTypeEnum);
	femmodel->parameters->FindParam(&isFS,FlowequationIsFSEnum);
	femmodel->parameters->FindParam(&dakota_analysis,QmuIsdakotaEnum);
	femmodel->parameters->SetParam(false,SaveResultsEnum);

	/*out of solution_type, figure out solution core and adjoint function pointer*/
	CorePointerFromSolutionEnum(&solutioncore,femmodel->parameters,solution_type);
	AdjointCorePointerFromSolutionEnum(&adjointcore,solution_type);

	/*Launch once a complete solution to set up all inputs*/
	if(VerboseControl()) _printf0_("   preparing initial solution\n");
	if(isFS) solutioncore(femmodel);

	/*Get initial guess*/
	Vector<IssmDouble> *Xpetsc = NULL;
	GetVectorFromControlInputsx(&Xpetsc,femmodel->elements,femmodel->nodes,femmodel->vertices,femmodel->loads,femmodel->materials,femmodel->parameters,"value");
	IssmDouble* X0 = Xpetsc->ToMPISerial();
	Xpetsc->GetSize(&nsize);
	delete Xpetsc;

	/*Initialize some of the BrentSearch arguments: */
	OptPars optpars;
	optpars.xmin    = 0; 
	optpars.xmax    = 1;
	optpars.nsteps  = nsteps;
	optpars.nsize   = nsize;
	optpars.maxiter = maxiter;
	optpars.cm_jump = cm_jump;

	/*Initialize function argument*/
	AppCtx usr;
	usr.femmodel = femmodel;
	usr.nsize    = nsize;

	/*Call Brent optimization*/
	BrentSearch(&J,optpars,X0,&FormFunction,&FormFunctionGradient,(void*)&usr);

	if(VerboseControl()) _printf0_("   preparing final solution\n");
	IssmDouble  *XL = NULL;
	IssmDouble  *XU = NULL;
	GetVectorFromControlInputsx(&XL,NULL,femmodel->elements,femmodel->nodes,femmodel->vertices,femmodel->loads,femmodel->materials,femmodel->parameters,"lowerbound");
	GetVectorFromControlInputsx(&XU,NULL,femmodel->elements,femmodel->nodes,femmodel->vertices,femmodel->loads,femmodel->materials,femmodel->parameters,"upperbound");
	for(long i=0;i<nsize;i++){
		if(X0[i]>XU[i]) X0[i]=XU[i];
		if(X0[i]<XL[i]) X0[i]=XL[i];
	}
	xDelete<IssmDouble>(XU);
	xDelete<IssmDouble>(XL);
	SetControlInputsFromVectorx(femmodel,X0);
	femmodel->parameters->SetParam(true,SaveResultsEnum);
	solutioncore(femmodel);

	/*some results not computed by steadystate_core or stressbalance_core: */
	if(!dakota_analysis){ //do not save this if we are running the control core from a qmu run!
		femmodel->OutputControlsx(&femmodel->results);

		#ifdef _HAVE_AD_
		IssmPDouble* J_passive=xNew<IssmPDouble>(nsteps);
		for(int i=0;i<nsteps;i++) J_passive[i]=reCast<IssmPDouble>(J[i]);
		femmodel->results->AddObject(new GenericExternalResult<IssmPDouble*>(femmodel->results->Size()+1,JEnum,J_passive,nsteps,1,0,0));
		xDelete<IssmPDouble>(J_passive);
		#else
		femmodel->results->AddObject(new GenericExternalResult<IssmPDouble*>(femmodel->results->Size()+1,JEnum,J,nsteps,1,0,0));
		#endif
	}

	/*Free resources: */
	xDelete<int>(control_type);
	xDelete<int>(maxiter);
	xDelete<IssmDouble>(cm_jump);
	xDelete<IssmDouble>(J);
	xDelete<IssmDouble>(X0);
}/*}}}*/
IssmDouble FormFunction(IssmDouble* X,void* usrvoid){/*{{{*/

	/*output: */
	IssmDouble J;

	/*parameters: */
	int        solution_type,analysis_type,num_responses;
	bool       conserve_loads = true;
	AppCtx*    usr = (AppCtx*)usrvoid;
	FemModel  *femmodel  = usr->femmodel;
	int        nsize     = usr->nsize;

	/*Recover parameters: */
	femmodel->parameters->FindParam(&analysis_type,AnalysisTypeEnum);
	femmodel->parameters->FindParam(&solution_type,SolutionTypeEnum);
	femmodel->parameters->FindParam(&num_responses,InversionNumCostFunctionsEnum);

	/*Constrain input vector*/
	IssmDouble  *XL = NULL;
	IssmDouble  *XU = NULL;
	GetVectorFromControlInputsx(&XL,NULL,femmodel->elements,femmodel->nodes,femmodel->vertices,femmodel->loads,femmodel->materials,femmodel->parameters,"lowerbound");
	GetVectorFromControlInputsx(&XU,NULL,femmodel->elements,femmodel->nodes,femmodel->vertices,femmodel->loads,femmodel->materials,femmodel->parameters,"upperbound");
	for(long i=0;i<nsize;i++){
		if(X[i]>XU[i]) X[i]=XU[i];
		if(X[i]<XL[i]) X[i]=XL[i];
	}

	/*Update control input*/
	SetControlInputsFromVectorx(femmodel,X);

	/*solve forward: */
	switch(solution_type){
		case SteadystateSolutionEnum:
			femmodel->SetCurrentConfiguration(StressbalanceAnalysisEnum);
			stressbalance_core(femmodel);	//We need a 3D velocity!! (vz is required for the next thermal run)
			break;
		case StressbalanceSolutionEnum:{
			femmodel->SetCurrentConfiguration(StressbalanceAnalysisEnum);

			bool is_schur_cg_solver = false;
			#ifdef _HAVE_PETSC_
			int solver_type;
			PetscOptionsDetermineSolverType(&solver_type);
			if(solver_type==FSSolverEnum) is_schur_cg_solver = true;
			#endif

			if(is_schur_cg_solver){
			 solutionsequence_schurcg(femmodel);
			}else{
			 solutionsequence_nonlinear(femmodel,conserve_loads); 
			}
			}
			 break;
		case BalancethicknessSolutionEnum:
			femmodel->SetCurrentConfiguration(BalancethicknessAnalysisEnum);
			solutionsequence_linear(femmodel); 
			break;
		case BalancethicknessSoftSolutionEnum:
			/*NOTHING*/
			break;
		case Balancethickness2SolutionEnum:
			femmodel->SetCurrentConfiguration(Balancethickness2AnalysisEnum);
			solutionsequence_linear(femmodel); 
			break;
		default:
			_error_("Solution " << EnumToStringx(solution_type) << " not implemented yet");
	}

	/*Compute misfit for this velocity field.*/
	IssmDouble* Jlist = NULL;
	femmodel->CostFunctionx(&J,&Jlist,NULL);
	_printf0_("f(x) = "<<setw(12)<<setprecision(7)<<J<<"  |  ");
	for(int i=0;i<num_responses;i++) _printf0_(" "<<setw(12)<<setprecision(7)<<Jlist[i]);
	_printf0_("\n");

	/*Free resources:*/
	xDelete<IssmDouble>(XU);
	xDelete<IssmDouble>(XL);
	xDelete<IssmDouble>(Jlist);
	return J;
}/*}}}*/
IssmDouble FormFunctionGradient(IssmDouble** pG,IssmDouble* X,void* usrvoid){/*{{{*/

	/*output: */
	IssmDouble J;
	int        temp;

	/*parameters: */
	void (*adjointcore)(FemModel*)=NULL;
	int         solution_type,analysis_type,num_responses,num_controls,numvertices;
	bool        conserve_loads = true;
	IssmDouble *scalar_list    = NULL;
	IssmDouble *Jlist          = NULL;
	IssmDouble *G              = NULL;
	IssmDouble *norm_list      = NULL;
	AppCtx     *usr            = (AppCtx*)usrvoid;
	FemModel   *femmodel       = usr->femmodel;
	int         nsize          = usr->nsize;

	/*Recover parameters: */
	femmodel->parameters->FindParam(&analysis_type,AnalysisTypeEnum);
	femmodel->parameters->FindParam(&solution_type,SolutionTypeEnum);
	femmodel->parameters->FindParam(&num_responses,InversionNumCostFunctionsEnum);
	femmodel->parameters->FindParam(&scalar_list,&temp,&temp,InversionGradientScalingEnum);
	femmodel->parameters->FindParam(&num_controls,InversionNumControlParametersEnum);   _assert_(num_controls);
	numvertices=femmodel->vertices->NumberOfVertices();

	/*Constrain input vector*/
	IssmDouble  *XL = NULL;
	IssmDouble  *XU = NULL;
	GetVectorFromControlInputsx(&XL,NULL,femmodel->elements,femmodel->nodes,femmodel->vertices,femmodel->loads,femmodel->materials,femmodel->parameters,"lowerbound");
	GetVectorFromControlInputsx(&XU,NULL,femmodel->elements,femmodel->nodes,femmodel->vertices,femmodel->loads,femmodel->materials,femmodel->parameters,"upperbound");
	for(long i=0;i<nsize;i++){
		if(X[i]>XU[i]) X[i]=XU[i];
		if(X[i]<XL[i]) X[i]=XL[i];
	}

	/*Update control input*/
	SetControlInputsFromVectorx(femmodel,X);

	/*Compute new temperature at this point*/
	if(solution_type==SteadystateSolutionEnum) steadystate_core(femmodel);

	/*Compute Adjoint*/
	AdjointCorePointerFromSolutionEnum(&adjointcore,solution_type);
	adjointcore(femmodel);

	/*Compute gradient*/
	Gradjx(&G,&norm_list,femmodel->elements,femmodel->nodes,femmodel->vertices,femmodel->loads,femmodel->materials,femmodel->parameters);

	/*Compute scaling factor*/
	IssmDouble scalar = scalar_list[0]/norm_list[0];
	for(int i=1;i<num_controls;i++) scalar=min(scalar,scalar_list[i]/norm_list[i]);

	/*Constrain Gradient*/
	for(int i=0;i<nsize;i++){
		G[i] = scalar*G[i];
	}

	for(long i=0;i<nsize;i++){
		if(X[i]>=XU[i]) G[i]=0.;
		if(X[i]<=XL[i]) G[i]=0.;
	}

	/*Needed for output results*/
	ControlInputSetGradientx(femmodel->elements,femmodel->nodes,femmodel->vertices,femmodel->loads,femmodel->materials,femmodel->parameters,G);

	/*Compute misfit for this velocity field.*/
	femmodel->CostFunctionx(&J,&Jlist,NULL);
	_printf0_("f(x) = "<<setw(12)<<setprecision(7)<<J<<"  |  ");
	for(int i=0;i<num_responses;i++) _printf0_(" "<<setw(12)<<setprecision(7)<<Jlist[i]);
	_printf0_("\n");

	/*Clean-up and return*/
	xDelete<IssmDouble>(XU);
	xDelete<IssmDouble>(XL);
	xDelete<IssmDouble>(norm_list);
	xDelete<IssmDouble>(scalar_list);
	xDelete<IssmDouble>(Jlist);
	*pG = G;
	return J;
}/*}}}*/

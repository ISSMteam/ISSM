/*!\file: controlm1qn3_core.cpp
 * \brief: core of the control solution 
 */ 

#include <config.h>
#include "./cores.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"
#include "../solutionsequences/solutionsequences.h"

#if defined (_HAVE_M1QN3_)
/*m1qn3 prototypes {{{*/
extern "C" void *ctonbe_; // DIS mode : Conversion
extern "C" void *ctcabe_; // DIS mode : Conversion
extern "C" void *euclid_; // Scalar product
typedef void (*SimulFunc) (long* indic,long* n, double* x, double* pf,double* g,long [],float [],void* dzs);
extern "C" void m1qn3_ (void f(long* indic,long* n, double* x, double* pf,double* g,long [],float [],void* dzs),
			void **, void **, void **,
			long *, double [], double *, double [], double*, double *,
			double *, char [], long *, long *, long *, long *, long *, long *, long [], double [], long *,
			long *, long *, long [], float [],void* );

/*Cost function prototype*/
void simul(long* indic,long* n,double* X,double* pf,double* G,long izs[1],float rzs[1],void* dzs);

/*Use struct to provide arguments*/
typedef struct{
	FemModel   * femmodel;
	IssmPDouble* Jlist;
	int         M;
	int         N;
	int*        i;
} m1qn3_struct;
/*}}}*/

void controlm1qn3_core(FemModel* femmodel){/*{{{*/

	/*Intermediaries*/
	long    omode;
	double  f,dxmin,dfmin_frac,gttol; 
	int     maxsteps,maxiter;
	int     intn,num_controls,num_cost_functions,solution_type;
	double *scaling_factors = NULL;
	double *X  = NULL;
	double *G  = NULL;

	/*Get control sizes*/
	int* M = NULL;
	int* N = NULL;
	femmodel->parameters->FindParam(&M,NULL,ControlInputSizeMEnum);
	femmodel->parameters->FindParam(&N,NULL,ControlInputSizeNEnum);

	/*Recover some parameters*/
	femmodel->parameters->FindParam(&solution_type,SolutionTypeEnum);
	femmodel->parameters->FindParam(&num_controls,InversionNumControlParametersEnum);
	femmodel->parameters->FindParam(&num_cost_functions,InversionNumCostFunctionsEnum);
	femmodel->parameters->FindParam(&maxsteps,InversionMaxstepsEnum);
	femmodel->parameters->FindParam(&maxiter,InversionMaxiterEnum);
	femmodel->parameters->FindParamAndMakePassive(&dxmin,InversionDxminEnum);
	femmodel->parameters->FindParamAndMakePassive(&dfmin_frac,InversionDfminFracEnum);
	femmodel->parameters->FindParamAndMakePassive(&gttol,InversionGttolEnum);
	femmodel->parameters->FindParamAndMakePassive(&scaling_factors,NULL,InversionControlScalingFactorsEnum);
	femmodel->parameters->SetParam(false,SaveResultsEnum);

	/*Initialize M1QN3 parameters*/
	if(VerboseControl())_printf0_("   Initialize M1QN3 parameters\n");
	SimulFunc costfuncion  = &simul;    /*Cost function address*/
	void**    prosca       = &euclid_;  /*Dot product function (euclid is the default)*/
	char      normtype[]   = "dfn";     /*Norm type: dfn = scalar product defined by prosca*/
	long      izs[5];                   /*Arrays used by m1qn3 subroutines*/
	long      iz[5];                    /*Integer m1qn3 working array of size 5*/
	float     rzs[1];                   /*Arrays used by m1qn3 subroutines*/
	long      impres       = 0;         /*verbosity level*/
	long      imode[3]     = {0};       /*scaling and starting mode, 0 by default*/
	long      indic        = 4;         /*compute f and g*/
	long      reverse      = 0;         /*reverse or direct mode*/
	long      io           = 6;         /*Channel number for the output*/

	/*Optimization criterions (need to be cast to long for m1qn3)*/
	long niter = long(maxsteps); /*Maximum number of iterations*/
	long nsim  = long(maxiter);  /*Maximum number of function calls*/

	/*Get initial guess*/
	GetPassiveVectorFromControlInputsx(&X,&intn,femmodel->elements,femmodel->nodes,femmodel->vertices,femmodel->loads,femmodel->materials,femmodel->parameters,"value");

	/*Get problem dimension and initialize gradient and initial guess*/
	long n = long(intn);
	G = xNew<double>(n);

	/*Scale control for M1QN3*/
	int offset = 0;
	for(int c=0;c<num_controls;c++){
		for(int i=0;i<M[c]*N[c];i++){
			int index = offset+i;
			X[index] = X[index]/reCast<IssmPDouble>(scaling_factors[c]);
		}
		offset += M[c]*N[c];
	}

	/*Allocate m1qn3 working arrays (see documentation)*/
	long      m   = 100;
	long      ndz = 4*n+m*(2*n+1);
	double*   dz  = xNew<double>(ndz);

	if(VerboseControl())_printf0_("   Computing initial solution\n");
	InversionStatsHeader(num_cost_functions);

	/*Prepare structure for m1qn3*/
	m1qn3_struct mystruct;
	mystruct.femmodel = femmodel;
	mystruct.M        = maxiter;
	mystruct.N        = num_cost_functions+1;
	mystruct.Jlist    = xNewZeroInit<IssmPDouble>(mystruct.M*mystruct.N);
	mystruct.i        = xNewZeroInit<int>(1);

	/*Initialize Gradient and cost function of M1QN3*/
	indic = 4; /*gradient required*/
	simul(&indic,&n,X,&f,G,izs,rzs,(void*)&mystruct);

	/*Estimation of the expected decrease in f during the first iteration*/
	if(dfmin_frac==0.) dfmin_frac=1.;
	double df1=dfmin_frac*f;

	/*Call M1QN3 solver*/
	m1qn3_(costfuncion,prosca,&ctonbe_,&ctcabe_,
				&n,X,&f,G,&dxmin,&df1,
				&gttol,normtype,&impres,&io,imode,&omode,&niter,&nsim,iz,dz,&ndz,
				&reverse,&indic,izs,rzs,(void*)&mystruct);

	/*Print exit flag*/
	InversionStatsFooter(num_cost_functions);
	switch(int(omode)){
		case 0:  _printf0_("   Stop requested (indic = 0)\n"); break;
		case 1:  _printf0_("   Convergence reached (gradient satisfies stopping criterion)\n"); break;
		case 2:  _printf0_("   Bad initialization\n"); break;
		case 3:  _printf0_("   Line search failure\n"); break;
		case 4:  _printf0_("   Maximum number of iterations exceeded\n");break;
		case 5:  _printf0_("   Maximum number of function calls exceeded\n"); break;
		case 6:  _printf0_("   stopped on dxmin during line search\n"); break;
		case 7:  _printf0_("   <g,d> > 0  or  <y,s> <0\n"); break;
		default: _printf0_("   Unknown end condition\n");
	}

	/*Constrain solution vector*/
	double  *XL = NULL;
	double  *XU = NULL;
	GetPassiveVectorFromControlInputsx(&XL,NULL,femmodel->elements,femmodel->nodes,femmodel->vertices,femmodel->loads,femmodel->materials,femmodel->parameters,"lowerbound");
	GetPassiveVectorFromControlInputsx(&XU,NULL,femmodel->elements,femmodel->nodes,femmodel->vertices,femmodel->loads,femmodel->materials,femmodel->parameters,"upperbound");

	offset = 0;
	for(int c=0;c<num_controls;c++){
		for(int i=0;i<M[c]*N[c];i++){
			int index = offset+i;
			X[index] = X[index]*scaling_factors[c];
			if(X[index]>XU[index]) X[index]=XU[index];
			if(X[index]<XL[index]) X[index]=XL[index];
		}
		offset += M[c]*N[c];
	}

	/*Set X as our new control (need to recast)*/
	#ifdef _HAVE_AD_
	IssmDouble* aX=xNew<IssmDouble>(intn);
	IssmDouble* aG=xNew<IssmDouble>(intn);
	for(int i=0;i<intn;i++) {
		aX[i] = reCast<IssmDouble>(X[i]); 
		aG[i] = reCast<IssmDouble>(G[i]);
	}
	ControlInputSetGradientx(femmodel->elements,femmodel->nodes,femmodel->vertices,femmodel->loads,femmodel->materials,femmodel->parameters,aG);
	SetControlInputsFromVectorx(femmodel,aX);
	xDelete(aX);
	xDelete(aG);
	#else
	SetControlInputsFromVectorx(femmodel,X);
	ControlInputSetGradientx(femmodel->elements,femmodel->nodes,femmodel->vertices,femmodel->loads,femmodel->materials,femmodel->parameters,G);
	#endif

	femmodel->OutputControlsx(&femmodel->results);
	femmodel->results->AddObject(new GenericExternalResult<double*>(femmodel->results->Size()+1,JEnum,mystruct.Jlist,(*mystruct.i),mystruct.N,0,0));
	femmodel->results->AddObject(new GenericExternalResult<int>(femmodel->results->Size()+1,InversionStopFlagEnum,int(omode)));

	/*Finalize*/
	if(VerboseControl()) _printf0_("   preparing final solution\n");
	femmodel->parameters->SetParam(true,SaveResultsEnum);
	void (*solutioncore)(FemModel*)=NULL;
	CorePointerFromSolutionEnum(&solutioncore,femmodel->parameters,solution_type);
	solutioncore(femmodel);

	/*Clean-up and return*/
	xDelete<int>(M);
	xDelete<int>(N);
	xDelete<double>(G);
	xDelete<double>(X);
	xDelete<double>(dz);
	xDelete<double>(XU);
	xDelete<double>(XL);
	xDelete<double>(scaling_factors);
	xDelete<double>(mystruct.Jlist);
	xDelete<int>(mystruct.i);
}/*}}}*/
void simul(long* indic,long* n,double* X,double* pf,double* G,long izs[1],float rzs[1],void* dzs){/*{{{*/

	/*Recover Arguments*/
	m1qn3_struct *input_struct = (m1qn3_struct*)dzs;
	FemModel     *femmodel     = input_struct->femmodel;
	IssmPDouble  *Jlist        = input_struct->Jlist;
	int           JlistM       = input_struct->M;
	int           JlistN       = input_struct->N;
	int          *Jlisti       = input_struct->i;

	/*Recover some parameters*/
	int num_responses,num_controls,solution_type;
	double* scaling_factors = NULL;
	int* M = NULL;
	int* N = NULL;
	femmodel->parameters->FindParam(&num_responses,InversionNumCostFunctionsEnum);
	femmodel->parameters->FindParam(&N,NULL,ControlInputSizeNEnum);
	femmodel->parameters->FindParam(&M,NULL,ControlInputSizeMEnum);
	femmodel->parameters->FindParam(&num_controls,InversionNumControlParametersEnum);
	femmodel->parameters->FindParamAndMakePassive(&scaling_factors,NULL,InversionControlScalingFactorsEnum);
	femmodel->parameters->FindParam(&solution_type,SolutionTypeEnum);

	/*Constrain input vector and update controls*/
	double *XL = NULL;
	double *XU = NULL;
	GetPassiveVectorFromControlInputsx(&XL,NULL,femmodel->elements,femmodel->nodes,femmodel->vertices,femmodel->loads,femmodel->materials,femmodel->parameters,"lowerbound");
	GetPassiveVectorFromControlInputsx(&XU,NULL,femmodel->elements,femmodel->nodes,femmodel->vertices,femmodel->loads,femmodel->materials,femmodel->parameters,"upperbound");

	int offset = 0;
	for(int c=0;c<num_controls;c++){
		for(int i=0;i<M[c]*N[c];i++){
			int index = offset+i;
			X[index] = X[index]*scaling_factors[c];
			if(X[index]>XU[index]) X[index]=XU[index];
			if(X[index]<XL[index]) X[index]=XL[index];
		}
		offset += M[c]*N[c];
	}

	#ifdef _HAVE_AD_
	IssmDouble* aX=xNew<IssmDouble>(*n);
	for(int i=0;i<*n;i++) aX[i] = reCast<IssmDouble>(X[i]); 
	SetControlInputsFromVectorx(femmodel,aX);
	xDelete(aX);
	#else
	SetControlInputsFromVectorx(femmodel,X);
	#endif

	/*Compute solution and adjoint*/
	void (*solutioncore)(FemModel*)=NULL;
	void (*adjointcore)(FemModel*)=NULL;
	CorePointerFromSolutionEnum(&solutioncore,femmodel->parameters,solution_type);
	solutioncore(femmodel);

	/*Check size of Jlist to avoid crashes*/
	_assert_((*Jlisti)<JlistM);
	_assert_(JlistN==num_responses+1);

	/*Compute objective function*/
	IssmDouble* Jtemp = NULL;
	IssmDouble  J;
	femmodel->CostFunctionx(&J,&Jtemp,NULL);
	*pf = reCast<double>(J);

	/*Record cost function values and delete Jtemp*/
	for(int i=0;i<num_responses;i++) Jlist[(*Jlisti)*JlistN+i] = reCast<double>(Jtemp[i]);
	Jlist[(*Jlisti)*JlistN+num_responses] = *pf;
	xDelete<IssmDouble>(Jtemp);

	if(*indic==0){
		/*dry run, no gradient required*/
		InversionStatsIter( (*Jlisti)+1, *pf, NAN, &Jlist[(*Jlisti)*JlistN], num_responses);

		*Jlisti = (*Jlisti) +1;
		xDelete<double>(XU);
		xDelete<double>(XL);
		return;
	}

	/*Compute Adjoint*/
	AdjointCorePointerFromSolutionEnum(&adjointcore,solution_type);
	adjointcore(femmodel);

	/*Compute gradient*/
	IssmDouble* G2 = NULL;
	Gradjx(&G2,NULL,femmodel->elements,femmodel->nodes,femmodel->vertices,femmodel->loads,femmodel->materials,femmodel->parameters);
	for(long i=0;i<*n;i++) G[i] = -reCast<double>(G2[i]);
	xDelete<IssmDouble>(G2);

	/*Constrain Gradient*/
	IssmDouble  Gnorm = 0.;
	offset = 0;
	for(int c=0;c<num_controls;c++){
		for(int i=0;i<M[c]*N[c];i++){
			int index = offset+i;
			if(X[index]>=XU[index]) G[index]=0.;
			if(X[index]<=XL[index]) G[index]=0.;
			G[index] = G[index]*scaling_factors[c];
			X[index] = X[index]/scaling_factors[c];
			Gnorm += G[index]*G[index];
		}
		offset += M[c]*N[c];
	}
	Gnorm = sqrt(Gnorm);

	/*Print info*/
	InversionStatsIter( (*Jlisti)+1, *pf, reCast<double>(Gnorm), &Jlist[(*Jlisti)*JlistN], num_responses);

	/*Clean-up and return*/
	*Jlisti = (*Jlisti) +1;
	xDelete<int>(M);
	xDelete<int>(N);
	xDelete<double>(XU);
	xDelete<double>(XL);
	xDelete<double>(scaling_factors);
}/*}}}*/

#else
void controlm1qn3_core(FemModel* femmodel){_error_("M1QN3 not installed");}
#endif //_HAVE_M1QN3_

/*!\file: controladm1qn3_core.cpp
 * \brief: core of the control solution
 */
#include <ctime>
#include <config.h>
#include "./cores.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"
#include "../solutionsequences/solutionsequences.h"

#if _HAVE_CODIPACK_
#include "../toolkits/codipack/CoDiPackGlobal.h"
#include <fenv.h>
double transient_ad(FemModel* femmodel, double* G,double* Jlist);
#endif

#if defined (_HAVE_M1QN3_) && defined(_HAVE_AD_)
/*m1qn3 prototypes {{{*/
extern "C" void *ctonbe_; // DIS mode : Conversion
extern "C" void *ctcabe_; // DIS mode : Conversion
extern "C" void *euclid_; // Scalar product
typedef void (*SimulFunc) (long* indic,long* n, double* x,double* pf,double* g,long [],float [],void* dzs);
extern "C" void m1qn3_ (void f(long* indic,long* n, double* x,double* pf,double* g,long [],float [],void* dzs),
			void **, void **, void **,
			long *, double [],double *, double[], double*, double *,
			double *, char [], long *, long *, long *, long *, long *, long *, long [],double [], long *,
			long *, long *, long [], float [],void* );

/*Use struct to provide arguments*/
typedef struct{
	FemModel   * femmodel;
	IssmPDouble* Jlist;
	int          M;
	int          N;
	int*         i;
} m1qn3_struct;
/*}}}*/

/*m1qm3 functions*/
void simul_starttrace(FemModel* femmodel){/*{{{*/

	#if defined(_HAVE_ADOLC_)
	/*Retrive ADOLC parameters*/
	IssmDouble gcTriggerRatio;
	IssmDouble gcTriggerMaxSize;
	IssmDouble obufsize;
	IssmDouble lbufsize;
	IssmDouble cbufsize;
	IssmDouble tbufsize;
	femmodel->parameters->FindParam(&gcTriggerRatio,AutodiffGcTriggerRatioEnum);
	femmodel->parameters->FindParam(&gcTriggerMaxSize,AutodiffGcTriggerMaxSizeEnum);
	femmodel->parameters->FindParam(&obufsize,AutodiffObufsizeEnum);
	femmodel->parameters->FindParam(&lbufsize,AutodiffLbufsizeEnum);
	femmodel->parameters->FindParam(&cbufsize,AutodiffCbufsizeEnum);
	femmodel->parameters->FindParam(&tbufsize,AutodiffTbufsizeEnum);

	/*Set garbage collection parameters: */
	setStoreManagerControl(reCast<IssmPDouble>(gcTriggerRatio),reCast<size_t>(gcTriggerMaxSize));

	/*Start trace: */
	int skipFileDeletion=1;
	int keepTaylors=1;
	int my_rank=IssmComm::GetRank();
	trace_on(my_rank,keepTaylors,reCast<size_t>(obufsize),reCast<size_t>(lbufsize),reCast<size_t>(cbufsize),reCast<size_t>(tbufsize),skipFileDeletion);

	#elif defined(_HAVE_CODIPACK_)

		//fprintf(stderr, "*** Codipack IoModel::StartTrace\n");
		codi_global.start();
	#else
	_error_("not implemented");
	#endif
}/*}}}*/
void simul_stoptrace(){/*{{{*/

	#if defined(_HAVE_ADOLC_)
	trace_off();
	if(VerboseAutodiff()){ /*{{{*/

		#ifdef _HAVE_ADOLC_
		int my_rank=IssmComm::GetRank();
		size_t  tape_stats[15];
		tapestats(my_rank,tape_stats); //reading of tape statistics
		int commSize=IssmComm::GetSize();
		int *sstats=new int[7];
		sstats[0]=tape_stats[NUM_OPERATIONS];
		sstats[1]=tape_stats[OP_FILE_ACCESS];
		sstats[2]=tape_stats[NUM_LOCATIONS];
		sstats[3]=tape_stats[LOC_FILE_ACCESS];
		sstats[4]=tape_stats[NUM_VALUES];
		sstats[5]=tape_stats[VAL_FILE_ACCESS];
		sstats[6]=tape_stats[TAY_STACK_SIZE];
		int *rstats=NULL;
		if (my_rank==0) rstats=new int[commSize*7];
		ISSM_MPI_Gather(sstats,7,ISSM_MPI_INT,rstats,7,ISSM_MPI_INT,0,IssmComm::GetComm());
		if (my_rank==0) {
			int offset=50;
			int rOffset=(commSize/10)+1;
			_printf_("   ADOLC statistics: \n");
			_printf_("     "<<setw(offset)<<left<<"#independents: " <<setw(12)<<right<<tape_stats[NUM_INDEPENDENTS] << "\n");
			_printf_("     "<<setw(offset)<<left<<"#dependents: " <<setw(12)<<right<<tape_stats[NUM_DEPENDENTS] << "\n");
			_printf_("     "<<setw(offset)<<left<<"max #live active variables: " <<setw(12)<<right<<tape_stats[NUM_MAX_LIVES] << "\n");
			_printf_("     operations: entry size "<< sizeof(unsigned char) << " Bytes \n");
			_printf_("     "<<setw(offset)<<left<<"  #entries in buffer (AutodiffObufsizeEnum) " <<setw(12)<<right<<tape_stats[OP_BUFFER_SIZE] << "\n");
			for (int r=0;r<commSize;++r)
			 _printf_("       ["<<setw(rOffset)<<right<<r<<"]"<<setw(offset-rOffset-4)<<left<<" #entries total" <<setw(12)<<right<<rstats[r*7+0] << (rstats[r*7+1]?" ->file":"") << "\n");
			_printf_("     locations: entry size " << sizeof(locint) << " Bytes\n");
			_printf_("     "<<setw(offset)<<left<<"  #entries in buffer (AutodiffLbufsizeEnum) " <<setw(12)<<right<<tape_stats[LOC_BUFFER_SIZE] << "\n");
			for (int r=0;r<commSize;++r)
			 _printf_("       ["<<setw(rOffset)<<right<<r<<"]"<<setw(offset-rOffset-4)<<left<<" #entries total" <<setw(12)<<right<<rstats[r*7+2] << (rstats[r*7+3]?" ->file":"") << "\n");
			_printf_("     constant values: entry size " << sizeof(double) << " Bytes\n");
			_printf_("     "<<setw(offset)<<left<<"  #entries in buffer (AutodiffCbufsizeEnum) " <<setw(12)<<right<<tape_stats[VAL_BUFFER_SIZE] << "\n");
			for (int r=0;r<commSize;++r)
			 _printf_("       ["<<setw(rOffset)<<right<<r<<"]"<<setw(offset-rOffset-4)<<left<<" #entries total" <<setw(12)<<right<<rstats[r*7+4] << (rstats[r*7+5]?" ->file":"") << "\n");
			_printf_("     Taylor stack: entry size " << sizeof(revreal) << " Bytes\n");
			_printf_("     "<<setw(offset)<<left<<"  #entries in buffer (AutodiffTbufsizeEnum) " <<setw(12)<<right<<tape_stats[TAY_BUFFER_SIZE] << "\n");
			for (int r=0;r<commSize;++r)
			 _printf_("       ["<<setw(rOffset)<<right<<r<<"]"<<setw(offset-rOffset-4)<<left<<" #entries total" <<setw(12)<<right<<rstats[r*7+6] << (rstats[r*7+6]>tape_stats[TAY_BUFFER_SIZE]?" ->file":"") << "\n");
			delete []rstats;
		}
		delete [] sstats;
		#endif
	} /*}}}*/

	#elif defined(_HAVE_CODIPACK_)

	codi_global.stop();
	if(VerboseAutodiff()){
		int my_rank=IssmComm::GetRank();
		if(my_rank == 0) {
			codi_global.print(std::cout);
		}
	}
	#else
	_error_("not implemented");
	#endif
}/*}}}*/
void simul_ad(long* indic,long* n,double* X,double* pf,double* G,long izs[1],float rzs[1],void* dzs){/*{{{*/

	/*Get rank*/
	int my_rank=IssmComm::GetRank();

	/*Recover Arguments*/
	m1qn3_struct *input_struct = (m1qn3_struct*)dzs;

	FemModel* femmodel = input_struct->femmodel;
	int num_responses,num_controls,solution_type;
	femmodel->parameters->FindParam(&solution_type,SolutionTypeEnum);

	/*we need to make sure we do not modify femmodel at each iteration, make a copy*/
	femmodel = input_struct->femmodel->copy();

	IssmPDouble*  Jlist  = input_struct->Jlist;
	int           JlistM = input_struct->M;
	int           JlistN = input_struct->N;
	int*          Jlisti = input_struct->i;
	int           intn   = (int)*n;

	/*Recover some parameters*/
	double *scaling_factors = NULL;
	int    *M = NULL;
	int    *N = NULL;
	int    *control_enum    = NULL;
	int     checkpoint_frequency;
	femmodel->parameters->FindParam(&num_responses,InversionNumCostFunctionsEnum);
	femmodel->parameters->FindParam(&num_controls,InversionNumControlParametersEnum);
	femmodel->parameters->FindParamAndMakePassive(&scaling_factors,NULL,InversionControlScalingFactorsEnum);
	femmodel->parameters->FindParam(&N,NULL,ControlInputSizeNEnum);
	femmodel->parameters->FindParam(&M,NULL,ControlInputSizeMEnum);
	femmodel->parameters->FindParam(&control_enum,NULL,InversionControlParametersEnum);
	femmodel->parameters->FindParam(&checkpoint_frequency,SettingsCheckpointFrequencyEnum);

	/*Constrain input vector and update controls*/
	double  *XL = NULL;
	double  *XU = NULL;
	GetPassiveVectorFromControlInputsx(&XL,NULL,femmodel->elements,femmodel->nodes,femmodel->vertices,femmodel->loads,femmodel->materials,femmodel->parameters,"lowerbound");
	GetPassiveVectorFromControlInputsx(&XU,NULL,femmodel->elements,femmodel->nodes,femmodel->vertices,femmodel->loads,femmodel->materials,femmodel->parameters,"upperbound");

	int offset = 0;
	for(int c=0;c<num_controls;c++){
		for(int i=0;i<M[c]*N[c];i++){
			int index = offset+i;
			X[index] = X[index]*scaling_factors[c];
			if(X[index]>XU[index]) X[index]=XU[index];
			if(X[index]<XL[index]) X[index]=XL[index];
			_assert_(!xIsNan(X[index]));
		}
		offset += M[c]*N[c];
	}

	/*Special case: do we need to run AD with checkpointing?*/
	#ifdef _HAVE_CODIPACK_
	if(checkpoint_frequency && solution_type == TransientSolutionEnum){
		SetControlInputsFromVectorx(femmodel,X);
		*pf = transient_ad(femmodel, G, &Jlist[(*Jlisti)*JlistN]);
	}
	else
	#endif
	  {

		/*Start Tracing*/
		simul_starttrace(femmodel);
		/*Set X as our new control input and as INDEPENDENT*/
		#ifdef _HAVE_AD_
		IssmDouble* aX=xNew<IssmDouble>(intn,"t");
		#else
		IssmDouble* aX=xNew<IssmDouble>(intn);
		#endif

		if(my_rank==0){
			for(int i=0;i<intn;i++){
				#if defined(_HAVE_ADOLC_)
					aX[i]<<=X[i];
				#elif defined(_HAVE_CODIPACK_)
					aX[i]=X[i];
					codi_global.registerInput(aX[i]);
				#else
					_error_("not suppoted");
				#endif
			}
		}

		ISSM_MPI_Bcast(aX,intn,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());
		SetControlInputsFromVectorx(femmodel,aX);
		xDelete<IssmDouble>(aX);

		/*Compute solution (forward)*/
		void (*solutioncore)(FemModel*)=NULL;
		CorePointerFromSolutionEnum(&solutioncore,femmodel->parameters,solution_type);
		solutioncore(femmodel);

		/*Get Dependents*/
		int          num_dependents;
		IssmPDouble *dependents;
		IssmDouble   J = 0.;
		DataSet     *dependent_objects = ((DataSetParam*)femmodel->parameters->FindParamObject(AutodiffDependentObjectsEnum))->value;
		femmodel->parameters->FindParam(&num_dependents,AutodiffNumDependentsEnum);

		/*Go through our dependent variables, and compute the response:*/
		dependents=xNew<IssmPDouble>(num_dependents);
		int i=-1;
		for(Object* & object:dependent_objects->objects){
			i++;
			DependentObject* dep=xDynamicCast<DependentObject*>(object);

			/*Get cost function for this dependent*/
			dep->RecordResponsex(femmodel);
			IssmDouble output_value = dep->GetValue();
			dependents[i] = output_value.getValue();
			#if defined(_HAVE_ADOLC_)
			output_value>>=dependents[i];
			#endif

			J+=output_value;
		}

		#if defined(_HAVE_CODIPACK_)
		// TODO: Registration of output values is more fine grained for ADOL-c.
		codi_global.registerOutput(J);
		#endif

		/*Turning off trace tape*/
		simul_stoptrace();

		/*intermediary: */
		int          num_independents=intn;
		IssmPDouble *aWeightVector=NULL;
		IssmPDouble *weightVectorTimesJac=NULL;
		IssmPDouble *totalgradient=xNewZeroInit<IssmPDouble>(num_independents);

		/*if no dependents, no point in running a driver: */
		if(!(num_dependents*num_independents)) _error_("this is not allowed");

		/*for adolc to run in parallel, we 0 out on rank~=0. But we still keep track of num_dependents:*/
		int num_dependents_old   = num_dependents;
		int num_independents_old = num_independents;

		#if defined(_HAVE_ADOLC_)
		/*Get gradient for ADOLC {{{*/
		if(my_rank!=0){
			num_dependents   = 0;
			num_independents = 0;
		}

		/*get the EDF pointer:*/
		ext_diff_fct *anEDF_for_solverx_p=xDynamicCast<GenericParam<Adolc_edf> * >(femmodel->parameters->FindParamObject(AdolcParamEnum))->GetParameterValue().myEDF_for_solverx_p;

		/* these are always needed regardless of the interpreter */
		anEDF_for_solverx_p->dp_x=xNew<double>(anEDF_for_solverx_p->max_n);
		anEDF_for_solverx_p->dp_y=xNew<double>(anEDF_for_solverx_p->max_m);

		/* Ok, now we are going to call the fos_reverse in a loop on the index, from 0 to num_dependents, so
		 * as to generate num_dependents gradients: */
		for(int aDepIndex=0;aDepIndex<num_dependents_old;aDepIndex++){

			/*initialize direction index in the weights vector: */
			aWeightVector=xNewZeroInit<IssmPDouble>(num_dependents);
			if (my_rank==0) aWeightVector[aDepIndex]=1.;

			/*initialize output gradient: */
			weightVectorTimesJac=xNew<IssmPDouble>(num_independents);

			/*set the forward method function pointer: */
			#ifdef _HAVE_GSL_
			anEDF_for_solverx_p->fos_reverse=EDF_fos_reverse_for_solverx;
			#endif
			#ifdef _HAVE_MUMPS_
			anEDF_for_solverx_p->fos_reverse_iArr=fos_reverse_mumpsSolveEDF;
			#endif

			anEDF_for_solverx_p->dp_U=xNew<IssmPDouble>(anEDF_for_solverx_p->max_m);
			anEDF_for_solverx_p->dp_Z=xNew<IssmPDouble>(anEDF_for_solverx_p->max_n);

			/*call driver: */
			fos_reverse(my_rank,num_dependents,num_independents, aWeightVector, weightVectorTimesJac );

			/*Add to totalgradient: */
			if(my_rank==0) for(int i=0;i<num_independents;i++) {
				totalgradient[i]+=weightVectorTimesJac[i];
			}

			/*free resources :*/
			xDelete(weightVectorTimesJac);
			xDelete(aWeightVector);
		}
		/*}}}*/
		#elif defined(_HAVE_CODIPACK_)
		/*Get gradient for CoDiPack{{{*/
		if(VerboseAutodiff())_printf0_("   CoDiPack fos_reverse\n");
		if(my_rank==0) codi_global.setGradient(0, 1.0);
		codi_global.evaluate();
		codi_global.getFullGradient(totalgradient, num_independents);

		/*Clear tape*/
		codi_global.clear();
		/*}}}*/
		#else
		_error_("not suppoted");
		#endif

		/*Broadcast gradient to other ranks*/
		ISSM_MPI_Bcast(totalgradient,num_independents_old,ISSM_MPI_PDOUBLE,0,IssmComm::GetComm());

		/*Check size of Jlist to avoid crashes*/
		_assert_((*Jlisti)<JlistM);
		_assert_(JlistN==num_responses+1);

		/*Compute objective function and broadcast it to other cpus*/
		*pf = reCast<double>(J);
		ISSM_MPI_Bcast(pf,1,ISSM_MPI_PDOUBLE,0,IssmComm::GetComm());

		/*Record cost function values and delete Jtemp*/
		for(int i=0;i<num_responses;i++) Jlist[(*Jlisti)*JlistN+i] = dependents[i];
		Jlist[(*Jlisti)*JlistN+num_responses] = reCast<IssmPDouble>(J);

		if(*indic==0){
			/*dry run, no gradient required*/
			InversionStatsIter( (*Jlisti)+1, *pf, NAN, &Jlist[(*Jlisti)*JlistN], num_responses);

			*Jlisti = (*Jlisti) +1;
			xDelete<double>(XU);
			xDelete<double>(XL);
			return;
		}

		/*Compute gradient*/
		for(long i=0;i<num_independents_old;i++) G[i] = totalgradient[i];

		xDelete<IssmPDouble>(dependents);
		xDelete<IssmPDouble>(totalgradient);
	  }

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
	_assert_(!xIsNan(Gnorm));
	_assert_(!xIsInf(Gnorm));

	/*Print info*/
	InversionStatsIter( (*Jlisti)+1, *pf, reCast<double>(Gnorm), &Jlist[(*Jlisti)*JlistN], num_responses);

	/*Clean-up and return*/
	delete femmodel;
	*Jlisti = (*Jlisti) +1;
	xDelete<double>(XU);
	xDelete<double>(XL);
	xDelete<int>(control_enum);
	xDelete<int>(M);
	xDelete<int>(N);
	xDelete<double>(scaling_factors);
}/*}}}*/
void controladm1qn3_core(FemModel* femmodel){/*{{{*/

	/*Intermediaries*/
	long    omode;
	double  f,dxmin,dfmin_frac,gttol;
	int     maxsteps,maxiter;
	int     intn ,num_controls,num_cost_functions,solution_type;
	double *scaling_factors = NULL;
	double *X               = NULL;
	double *G               = NULL;
	int    *N               = NULL;
	int    *M               = NULL;
	int    *control_enum;

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
	femmodel->parameters->FindParam(&control_enum,NULL,InversionControlParametersEnum);
	femmodel->parameters->SetParam(false,SaveResultsEnum);
	femmodel->parameters->FindParam(&N,NULL,ControlInputSizeNEnum);
   femmodel->parameters->FindParam(&M,NULL,ControlInputSizeMEnum);

	/*Initialize M1QN3 parameters*/
	if(VerboseControl())_printf0_("   Initialize M1QN3 parameters\n");
	SimulFunc simul_ptr    = &simul_ad; /*Cost function address*/
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

	/*Optimization criterions*/
	long niter = long(maxsteps); /*Maximum number of iterations*/
	long nsim  = long(maxiter);/*Maximum number of function calls*/

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
         X[index] = X[index]/scaling_factors[c];
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
	simul_ad(&indic,&n,X,&f,G,izs,rzs,(void*)&mystruct);

	/*Estimation of the expected decrease in f during the first iteration*/
	if(dfmin_frac==0.) dfmin_frac=1.;
	double df1=dfmin_frac*f;

	/*Call M1QN3 solver*/
	m1qn3_(simul_ptr,prosca,&ctonbe_,&ctcabe_,
				&n,X,&f,G,&dxmin,&df1,
				&gttol,normtype,&impres,&io,imode,&omode,&niter,&nsim,iz,dz,&ndz,
				&reverse,&indic,izs,rzs,(void*)&mystruct);

	/*Print exit flag*/
	InversionStatsFooter(num_cost_functions);
	_printf0_("   Exit code "<<int(omode));
	switch(int(omode)){
		case 0:  _printf0_(": Stop requested (indic = 0)\n"); break;
		case 1:  _printf0_(": Convergence reached (gradient satisfies stopping criterion)\n"); break;
		case 2:  _printf0_(": Bad initialization\n"); break;
		case 3:  _printf0_(": Line search failure\n"); break;
		case 4:  _printf0_(": Maximum number of iterations exceeded\n");break;
		case 5:  _printf0_(": Maximum number of function calls exceeded\n"); break;
		case 6:  _printf0_(": stopped on dxmin during line search\n"); break;
		case 7:  _printf0_(": <g,d> > 0  or  <y,s> <0\n"); break;
		default: _printf0_(": Unknown end condition\n");
	}

	/*Constrain solution vector*/
	double  *XL = NULL;
	double  *XU = NULL;
	GetPassiveVectorFromControlInputsx(&XL,NULL,femmodel->elements,femmodel->nodes,femmodel->vertices,femmodel->loads,femmodel->materials,femmodel->parameters,"lowerbound");
	GetPassiveVectorFromControlInputsx(&XU,NULL,femmodel->elements,femmodel->nodes,femmodel->vertices,femmodel->loads,femmodel->materials,femmodel->parameters,"upperbound");

   offset = 0;
   for (int c=0;c<num_controls;c++){
      for(int i=0;i<M[c]*N[c];i++){
         int index = offset+i;
         X[index] = X[index]*scaling_factors[c];
         if(X[index]>XU[index]) X[index]=XU[index];
         if(X[index]<XL[index]) X[index]=XL[index];
      }
      offset += M[c]*N[c];
   }

	/*Set X as our new control*/
	IssmDouble* aX=xNew<IssmDouble>(intn);
	IssmDouble* aG=xNew<IssmDouble>(intn);

	for(int i=0;i<intn;i++) {
		aX[i] = reCast<IssmDouble>(X[i]);
		aG[i] = reCast<IssmDouble>(G[i]);
	}

	ControlInputSetGradientx(femmodel->elements,femmodel->nodes,femmodel->vertices,femmodel->loads,femmodel->materials,femmodel->parameters,aG);
	SetControlInputsFromVectorx(femmodel,aX);
	xDelete(aX);

	if (solution_type == TransientSolutionEnum){
		int step = 1;
		femmodel->parameters->SetParam(step,StepEnum);
		femmodel->results->AddObject(new GenericExternalResult<IssmPDouble*>(femmodel->results->Size()+1,JEnum,mystruct.Jlist,(*mystruct.i),mystruct.N));

		int offset = 0;
		for(int i=0;i<num_controls;i++){

			/*Disect results*/
			GenericExternalResult<IssmPDouble*>* G_output = new GenericExternalResult<IssmPDouble*>(femmodel->results->Size()+1,Gradient1Enum+i,&G[offset],N[i],M[i]);
			GenericExternalResult<IssmPDouble*>* X_output = new GenericExternalResult<IssmPDouble*>(femmodel->results->Size()+1,control_enum[i],&X[offset],N[i],M[i]);

			/*transpose for consistency with MATLAB's formating*/
			G_output->Transpose();
			X_output->Transpose();

			/*Add to results*/
			femmodel->results->AddObject(G_output);
			femmodel->results->AddObject(X_output);

			offset += N[i]*M[i];
		}
	}
	else{
		//FIXME: merge with code above?
		femmodel->results->AddObject(new GenericExternalResult<IssmPDouble*>(femmodel->results->Size()+1,JEnum,mystruct.Jlist,(*mystruct.i),mystruct.N,0,0));
		femmodel->OutputControlsx(&femmodel->results);
	}
	femmodel->results->AddObject(new GenericExternalResult<int>(femmodel->results->Size()+1,InversionStopFlagEnum,int(omode)));

	xDelete(aG);

	/*Add last cost function to results*/

	/*Finalize*/
	if(VerboseControl()) _printf0_("   preparing final solution\n");
	femmodel->parameters->SetParam(true,SaveResultsEnum);
	femmodel->parameters->SetParam(0,SettingsCheckpointFrequencyEnum);
	void (*solutioncore)(FemModel*)=NULL;
	CorePointerFromSolutionEnum(&solutioncore,femmodel->parameters,solution_type);
	solutioncore(femmodel);

	/*Clean-up and return*/
	xDelete<double>(G);
	xDelete<double>(X);
	xDelete<double>(dz);
	xDelete<double>(XU);
	xDelete<double>(XL);
	xDelete<double>(scaling_factors);
	xDelete<IssmPDouble>(mystruct.Jlist);
	xDelete<int>(mystruct.i);
	xDelete<int>(control_enum);
	xDelete<int>(M);
	xDelete<int>(N);
}/*}}}*/

#else
void controladm1qn3_core(FemModel* femmodel){_error_("M1QN3 or ADOLC/CoDiPack not installed");}
#endif //_HAVE_M1QN3_

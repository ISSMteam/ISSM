/*!\file: controlvalidation_core.cpp
 * \brief: core of the control solution
 */

#include "./cores.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"
#include "../toolkits/codipack/CoDiPackGlobal.h"

#ifdef _HAVE_AD_
void simul_starttrace2(FemModel* femmodel){/*{{{*/

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
void simul_stoptrace2(){/*{{{*/

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

	#if _CODIPACK_MAJOR_>=2
	auto& tape_codi = IssmDouble::getTape();
	#elif _CODIPACK_MAJOR_==1
	auto& tape_codi = IssmDouble::getGlobalTape();
	#else
	#error "_CODIPACK_MAJOR_ not supported"
	#endif

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
#endif

void controlvalidation_core(FemModel* femmodel){

	int         solution_type,n;
	int         num_responses;
	IssmDouble  j0,j;
	IssmDouble  Ialpha,exponent,alpha;
	IssmDouble* scaling_factors = NULL;
	IssmDouble* jlist = NULL;
	int my_rank=IssmComm::GetRank();

	/*Recover parameters used throughout the solution*/
	femmodel->parameters->FindParam(&solution_type,SolutionTypeEnum);
	femmodel->parameters->SetParam(false,SaveResultsEnum);
	femmodel->parameters->FindParam(&num_responses,InversionNumCostFunctionsEnum);
	femmodel->parameters->FindParam(&scaling_factors,NULL,InversionControlScalingFactorsEnum);

	/*Get initial guess*/
	IssmPDouble* X0 = NULL;
	GetPassiveVectorFromControlInputsx(&X0,&n,femmodel->elements,femmodel->nodes,femmodel->vertices,femmodel->loads,femmodel->materials,femmodel->parameters,"value");

	/*Allocate vectors*/
	IssmDouble*  X = xNew<IssmDouble>(n);
	IssmPDouble* G = xNew<IssmPDouble>(n);

	/*out of solution_type, figure out solution core and adjoint function pointer*/
	void (*solutioncore)(FemModel*)=NULL;
	CorePointerFromSolutionEnum(&solutioncore,femmodel->parameters,solution_type);

	#if defined(_HAVE_ADOLC_)
	/*{{{*/
	IssmDouble* aX=xNew<IssmDouble>(n);
	if(my_rank==0){
		for(int i=0;i<n;i++){
			aX[i]<<=X0[i];
		}
	}
	_error_("not implemented yet...");
	/*}}}*/
	#elif defined(_HAVE_CODIPACK_)
	/*{{{*/
	simul_starttrace2(femmodel);
	IssmDouble* aX=xNew<IssmDouble>(n);

	codi_global.start();
	if(my_rank==0){
		for (int i=0;i<n;i++) {
			aX[i]=X0[i];
			codi_global.registerInput(aX[i]);
		}
	}
	SetControlInputsFromVectorx(femmodel,aX);
	xDelete(aX);

	if(VerboseControl()) _printf0_("   Compute Initial cost function\n");
	solutioncore(femmodel);

	/*Get Dependents*/
	DataSet     *dependent_objects = ((DataSetParam*)femmodel->parameters->FindParamObject(AutodiffDependentObjectsEnum))->value;
	IssmDouble	J=0.;

	/*Go through our dependent variables, and compute the response:*/
	int i=-1;
	for(Object* & object:dependent_objects->objects){
		DependentObject* dep=xDynamicCast<DependentObject*>(object);
		i++;
		dep->RecordResponsex(femmodel);
		IssmDouble output_value = dep->GetValue();

		_printf0_("=== output ="<<output_value<<" \n");
		if(my_rank==0) {
			codi_global.registerOutput(output_value);
			J+=output_value;
		}
	}
	j0 = J;
	_printf0_("Initial cost function J(x) = "<<setw(12)<<setprecision(7)<<j0<<"\n");
	_assert_(j0>0.);
	simul_stoptrace2();
	/*initialize direction index in the weights vector: */
	if(my_rank==0){
		codi_global.setGradient(0, 1.0); // TODO: This is different form J. Also: Only one output is seeded and not all.
	}
	codi_global.evaluate();

	/*Get gradient for this dependent */
	codi_global.getFullGradient(G, n);

	/*Clear tape*/
	codi_global.clear();
/*}}}*/
	#else
	/*{{{*/
	void (*adjointcore)(FemModel*)  = NULL;
	AdjointCorePointerFromSolutionEnum(&adjointcore,solution_type);

	if(VerboseControl()) _printf0_("   Compute Initial solution\n");
	solutioncore(femmodel);
	if(VerboseControl()) _printf0_("   Compute Adjoint\n");
	adjointcore(femmodel);

	if(VerboseControl()) _printf0_("   Compute Initial cost function\n");
	femmodel->CostFunctionx(&j0,&jlist,NULL);
	_printf0_("Initial cost function J(x) = "<<setw(12)<<setprecision(7)<<j0<<"\n");
	xDelete<IssmDouble>(jlist);

	if(VerboseControl()) _printf0_("   Compute Gradient\n");
	Gradjx(&G,NULL,femmodel->elements,femmodel->nodes,femmodel->vertices,femmodel->loads,femmodel->materials,femmodel->parameters);
	for(int i=0;i<n;i++) G[i] = -G[i];
	/*}}}*/
	#endif

	/*Allocate output*/
	int num = 26;
	IssmDouble* output = xNew<IssmDouble>(2*num);

	/*Start loop*/
	_printf0_("       alpha      Ialpha \n");
	_printf0_("_________________________\n");
	for(int m=0;m<num;m++){

		/*Create new vector*/
		alpha    = pow(2.,-m);
		for(int i=0;i<n;i++) X[i] = X0[i] + alpha*scaling_factors[0];

		/*Calculate j(k+alpha delta k) */
		SetControlInputsFromVectorx(femmodel,X);
		solutioncore(femmodel);

		#if defined(_HAVE_CODIPACK_)
		j=0.;
		for(Object* & object:dependent_objects->objects){
			DependentObject* dep=xDynamicCast<DependentObject*>(object);
			dep->RecordResponsex(femmodel);
			IssmDouble output_value = dep->GetValue();
			j+=output_value;
		}
		#else
		femmodel->CostFunctionx(&j,NULL,NULL);
		#endif

		IssmDouble Den = 0.;
		for(int i=0;i<n;i++) Den += alpha* G[i] * scaling_factors[0];
		Ialpha = fabs((j - j0)/Den - 1.);
		_assert_(fabs(Den)>0.);

		_printf0_(" " << setw(11) << setprecision (5)<<alpha<<" " << setw(11) << setprecision (5)<<Ialpha<<"\n");
		output[m*2+0] = alpha;
		output[m*2+1] = Ialpha;
	}

	/*output*/
	#ifdef _HAVE_AD_
	IssmPDouble* J_passive=xNew<IssmPDouble>(2*num);
	for(int i=0;i<2*num;i++) J_passive[i]=reCast<IssmPDouble>(output[i]);
	femmodel->results->AddObject(new GenericExternalResult<IssmPDouble*>(femmodel->results->Size()+1,JEnum,J_passive,num,2,0,0));
	xDelete<IssmPDouble>(J_passive);
	IssmDouble* aG=xNew<IssmDouble>(n);
	for(int i=0;i<n;i++) aG[i] = G[i];
	ControlInputSetGradientx(femmodel->elements,femmodel->nodes,femmodel->vertices,femmodel->loads,femmodel->materials,femmodel->parameters,aG);
	xDelete<IssmDouble>(aG);
	#else
	femmodel->results->AddObject(new GenericExternalResult<IssmPDouble*>(femmodel->results->Size()+1,JEnum,output,num,2,0,0));
	ControlInputSetGradientx(femmodel->elements,femmodel->nodes,femmodel->vertices,femmodel->loads,femmodel->materials,femmodel->parameters,G);
	#endif
	femmodel->OutputControlsx(&femmodel->results);

	/*Clean up and return*/
	xDelete<IssmDouble>(output);
	xDelete<IssmPDouble>(G);
	xDelete<IssmDouble>(X);
	xDelete<double>(X0);
	xDelete<IssmDouble>(scaling_factors);
}

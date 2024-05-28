/*!\file ad_core
 * \brief: compute outputs from the AD mode,  using our dependents and independents, and drivers available in Adolc.
 */

/*Includes: {{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include <set>
#include "./cores.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"
#include "../solutionsequences/solutionsequences.h"
/*}}}*/

#ifdef _HAVE_CODIPACK_
CoDi_global codi_global = {};
#endif
void ad_core(FemModel* femmodel){

	/*diverse: */
	int     i;
	int     dummy;
	int     num_dependents=0;
	int     num_independents=0;
	bool    isautodiff,iscontrol;
	char   *driver           = NULL;
	size_t  tape_stats[15];

	/*state variables: */
	IssmDouble *axp = NULL;
	double     *xp  = NULL;
	int my_rank=IssmComm::GetRank();

	/*AD mode on?: */
	femmodel->parameters->FindParam(&isautodiff,AutodiffIsautodiffEnum);
	femmodel->parameters->FindParam(&iscontrol,InversionIscontrolEnum);

	if(isautodiff && !iscontrol){

		#if defined(_HAVE_ADOLC_)
			if(VerboseAutodiff())_printf0_("   start ad core\n");

			/*First, stop tracing: */
			trace_off();

			/*Print tape statistics so that user can kill this run if something is off already:*/
			if(VerboseAutodiff()){ /*{{{*/
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
			} /*}}}*/

			/*retrieve parameters: */
			femmodel->parameters->FindParam(&num_dependents,AutodiffNumDependentsEnum);
			femmodel->parameters->FindParam(&num_independents,AutodiffNumIndependentsEnum);

			/*if no dependents, no point in running a driver: */
			if(!(num_dependents*num_independents)) return;

			/*for adolc to run in parallel, we 0 out on rank~=0:*/
			if (my_rank!=0){
				num_dependents=0; num_independents=0;
			}

			/*retrieve state variable: */
			femmodel->parameters->FindParam(&axp,&dummy,AutodiffXpEnum);

			/* driver argument */
			xp=xNew<double>(num_independents);
			for(i=0;i<num_independents;i++){
				xp[i]=reCast<double,IssmDouble>(axp[i]);
			}

			/*get the EDF pointer:*/
			ext_diff_fct *anEDF_for_solverx_p=xDynamicCast<GenericParam<Adolc_edf> * >(femmodel->parameters->FindParamObject(AdolcParamEnum))->GetParameterValue().myEDF_for_solverx_p;

			/*Branch according to AD driver: */
			femmodel->parameters->FindParam(&driver,AutodiffDriverEnum);

			if (strcmp(driver,"fos_forward")==0){ /*{{{*/

				int     anIndepIndex;
				double *tangentDir         = NULL;
				double *jacTimesTangentDir = NULL;
				double *theOutput          = NULL;

				/*retrieve direction index: */
				femmodel->parameters->FindParam(&anIndepIndex,AutodiffFosForwardIndexEnum);

				if (anIndepIndex<0 || anIndepIndex>=num_independents) _error_("index value for AutodiffFosForwardIndexEnum should be in [0,num_independents-1]");

				tangentDir=xNewZeroInit<double>(num_independents);
				tangentDir[anIndepIndex]=1.0;

				jacTimesTangentDir=xNew<double>(num_dependents);
				theOutput=xNew<double>(num_dependents);

				/*set the forward method function pointer: */
#ifdef _HAVE_GSL_
				anEDF_for_solverx_p->fos_forward=EDF_fos_forward_for_solverx;
#endif

				/*call driver: */
				fos_forward(my_rank,num_dependents,num_independents, 0, xp, tangentDir, theOutput, jacTimesTangentDir );

				/*add to results*/
				femmodel->results->AddObject(new GenericExternalResult<IssmPDouble*>(femmodel->results->Size()+1,AutodiffJacobianEnum,jacTimesTangentDir,num_dependents,1,0,0.0));

				/*free resources :*/
				xDelete(theOutput);
				xDelete(jacTimesTangentDir);
				xDelete(tangentDir);
			} /*}}}*/
			else if ((strcmp(driver,"fov_forward")==0) || (strcmp(driver,"fov_forward_all")==0)){ /*{{{*/

				int      tangentDirNum;
				int      dummy;
				int     *indepIndices  = NULL;
				double **jacTimesSeed  = NULL;
				double **seed          = NULL;
				double  *theOutput     = NULL;
				std::set<unsigned int> anIndexSet;

				/*retrieve directions:*/
				if (strcmp(driver,"fov_forward_all")==0){
					tangentDirNum=num_independents;
					indepIndices=xNewZeroInit<int>(tangentDirNum);
					for(i=0;i<num_independents;i++)indepIndices[i]=1;
				}
				else{
					femmodel->parameters->FindParam(&indepIndices,&tangentDirNum,&dummy,AutodiffFovForwardIndicesEnum);
				}

				/*Some checks: */
				if (tangentDirNum<1 || tangentDirNum>num_independents) _error_("tangentDirNum should be in [1,num_independents]");

				/* full Jacobian or Jacobian projection:*/
				jacTimesSeed=xNew<double>(num_dependents,tangentDirNum);

				/*set the forward method function pointers: */
#ifdef _HAVE_GSL_
				anEDF_for_solverx_p->fov_forward=EDF_fov_forward_for_solverx;
#endif
				// anEDF_for_solverx_p->fov_reverse=EDF_fov_reverse_for_solverx;

				/*seed matrix: */
				seed=xNewZeroInit<double>(num_independents,tangentDirNum);

				/*collect indices in a set to prevent accidental duplicates as long as we don't do compression:*/
				for (int i=0; i<tangentDirNum; ++i) {
					/* make sure the index is in range*/
					if (indepIndices[i]>num_independents) {
						_error_("indepIndices values must be in [0,num_independents-1]");
					}
					if (anIndexSet.find(indepIndices[i])!=anIndexSet.end()) {
						_error_("duplicate indepIndices values are not allowed until we implement Jacobian decompression");
					}
					anIndexSet.insert(indepIndices[i]);
					/* now populate the seed matrix from the set of independent indices;
					 * simple setup with a single 1.0 per column and at most a single 1.0 per row*/
					seed[indepIndices[i]][i]=1.0;
				}

				/*allocate output: */
				theOutput=xNew<double>(num_dependents);

				/*call driver: */
				fov_forward(my_rank,num_dependents,num_independents, tangentDirNum, xp, seed, theOutput, jacTimesSeed );
				/*Free resources: */
				xDelete(theOutput);
				xDelete(indepIndices);
				xDelete(seed);

				/*add to results: */
				femmodel->results->AddObject(new GenericExternalResult<IssmPDouble*>(femmodel->results->Size()+1,AutodiffJacobianEnum,*jacTimesSeed,num_dependents*tangentDirNum,1,0,0.0));

				/*Free resources: */
				xDelete(jacTimesSeed);
				xDelete(indepIndices);
			} /*}}}*/
			else if (strcmp(driver,"fos_reverse")==0) { /*{{{*/

				int     aDepIndex=0;
				double *aWeightVector=NULL;
				double *weightVectorTimesJac=NULL;

				/*retrieve direction index: */
				femmodel->parameters->FindParam(&aDepIndex,AutodiffFosReverseIndexEnum);
				aWeightVector=xNewZeroInit<double>(num_dependents);
				if (my_rank==0) {
					if (aDepIndex<0 || aDepIndex>=num_dependents) _error_("index value for AutodiffFosReverseIndexEnum should be in [0,num_dependents-1]");
					aWeightVector[aDepIndex]=1.0;
				}
				weightVectorTimesJac=xNew<double>(num_independents);

				/*set the forward method function pointer: */
#ifdef _HAVE_GSL_
				anEDF_for_solverx_p->fos_reverse=EDF_fos_reverse_for_solverx;
#endif
#ifdef _HAVE_MUMPS_
				anEDF_for_solverx_p->fos_reverse_iArr=fos_reverse_mumpsSolveEDF;
#endif

				/*call driver: */
				fos_reverse(my_rank,num_dependents,num_independents, aWeightVector, weightVectorTimesJac );

				/*add to results*/
				femmodel->results->AddObject(new GenericExternalResult<IssmPDouble*>(femmodel->results->Size()+1,AutodiffJacobianEnum,weightVectorTimesJac,num_independents,1,0,0.0));

				/*free resources :*/
				xDelete(weightVectorTimesJac);
				xDelete(aWeightVector);
			} /*}}}*/
			else if ((strcmp(driver,"fov_reverse")==0) || (strcmp(driver,"fov_reverse_all")==0)){ /*{{{*/

				int* depIndices=NULL;
				int weightNum;
				int dummy;
				double **weightsTimesJac=NULL;
				double **weights=NULL;
				std::set<unsigned int> anIndexSet;

				/*retrieve directions:*/
				if (strcmp(driver,"fov_reverse_all")==0){
					weightNum=num_dependents;
					depIndices=xNewZeroInit<int>(weightNum);
					for(i=0;i<num_dependents;i++)depIndices[i]=1;
				}
				else{
					femmodel->parameters->FindParam(&depIndices,&weightNum,&dummy,AutodiffFovForwardIndicesEnum);
				}

				/*Some checks: */
				if (weightNum<1 || weightNum>num_dependents) _error_("tangentDirNum should be in [1,num_dependents]");

				/* full Jacobian or Jacobian projection:*/
				weightsTimesJac=xNew<double>(weightNum,num_independents);

				/*set the forward method function pointers: */
				#ifdef _HAVE_GSL_
				anEDF_for_solverx_p->fov_reverse=EDF_fov_reverse_for_solverx;
				#endif

				/*seed matrix: */
				weights=xNewZeroInit<double>(weightNum,num_dependents);

				/*collect indices in a set to prevent accidental duplicates as long as we don't do compression:*/
				for (int i=0; i<weightNum; ++i) {
					/* make sure the index is in range*/
					if (depIndices[i]>num_dependents) {
						_error_("depIndices values must be in [0,num_dependents-1]");
					}
					if (anIndexSet.find(depIndices[i])!=anIndexSet.end()) {
						_error_("duplicate depIndices values are not allowed until we implement Jacobian decompression");
					}
					anIndexSet.insert(depIndices[i]);
					/* now populate the seed matrix from the set of independent indices;
					 * simple setup with a single 1.0 per column and at most a single 1.0 per row*/
					weights[depIndices[i]][i]=1.0;
				}

				/*call driver: */
				fov_reverse(my_rank,num_dependents,num_independents, weightNum, weights, weightsTimesJac );

				/*add to results: */
				femmodel->results->AddObject(new GenericExternalResult<IssmPDouble*>(femmodel->results->Size()+1,AutodiffJacobianEnum,*weightsTimesJac,weightNum*num_independents,1,0,0.0));

				/*Free resources: */
				xDelete(weights);
				xDelete(weightsTimesJac);
				xDelete(depIndices);
			} /*}}}*/
			else _error_("driver: " << driver << " not yet supported!");

			if(VerboseAutodiff())_printf0_("   end AD core\n");

			/*Free resources: */
			xDelete(xp);
			xDelete(axp); 
			xDelete(driver);

			#elif defined(_HAVE_CODIPACK_)
			if(VerboseAutodiff())_printf0_("   start CoDiPack ad core\n");

			/*First, stop tracing: */
			#if _CODIPACK_MAJOR_==2
			auto& tape_codi = IssmDouble::getTape();
			#elif _CODIPACK_MAJOR_==1
			auto& tape_codi = IssmDouble::getGlobalTape();
			#else
			#error "_CODIPACK_MAJOR_ not supported"
			#endif

			tape_codi.setPassive();

			if(VerboseAutodiff()){ /*{{{*/
				if(my_rank == 0) {
					// FIXME codi "just because" for now
					tape_codi.printStatistics(std::cout);
					codi_global.print(std::cout);
				}
			}

			/*retrieve parameters: */
			femmodel->parameters->FindParam(&num_dependents,AutodiffNumDependentsEnum);
			femmodel->parameters->FindParam(&num_independents,AutodiffNumIndependentsEnum);

			/*if no dependents, no point in running a driver: */
			if(!(num_dependents*num_independents)) return;

			/*Branch according to AD driver: */
			femmodel->parameters->FindParam(&driver,AutodiffDriverEnum);
			if(VerboseAutodiff())_printf0_("   driver: " << driver << "\n");

			if (strcmp(driver,"fos_reverse")==0) { /*{{{*/
				if(VerboseAutodiff())_printf0_("   CoDiPack fos_reverse\n");
				int     aDepIndex=0;
				double *weightVectorTimesJac=NULL;

				/*retrieve direction index: */
				femmodel->parameters->FindParam(&aDepIndex,AutodiffFosReverseIndexEnum);
				if (my_rank==0) {
					if (aDepIndex<0 || aDepIndex>=num_dependents
							|| codi_global.output_indices.size() <= aDepIndex){
						_error_("index value for AutodiffFosReverseIndexEnum should be in [0,num_dependents-1]");
					}
					tape_codi.setGradient(codi_global.output_indices[aDepIndex], 1.0);
				}

				tape_codi.evaluate();

				weightVectorTimesJac=xNew<double>(num_independents);
				/*call driver: */
				auto in_size = codi_global.input_indices.size();
				for(size_t i = 0; i < in_size; ++i) {
					weightVectorTimesJac[i] = tape_codi.getGradient(codi_global.input_indices[i]);
				}

				/*add to results*/
				femmodel->results->AddObject(new GenericExternalResult<IssmPDouble*>(femmodel->results->Size()+1,AutodiffJacobianEnum,weightVectorTimesJac,num_independents,1,0,0.0));

				/*free resources :*/
				xDelete(weightVectorTimesJac);
			} /*}}}*/
			else _error_("driver: " << driver << " not yet supported!");

			if(VerboseAutodiff())_printf0_("   end CoDiPack ad core\n");
			xDelete(driver);
		#else
			_error_("Should not be requesting AD drivers when an AD library is not available!");
		#endif
	}
}

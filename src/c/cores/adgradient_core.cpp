/*!\file adgradient_core
 * \brief: compute gradient for all scalar depenendents, then sum them up as output. This relies mainly on the fos_reverse 
 * driver, hence the incapacity to merge this with ad_core.cpp.
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

void adgradient_core(FemModel* femmodel){

	/*diverse: */
	int     i;
	int     dummy;
	int     num_dependents=0;
	int     num_dependents_old=0;
	int     num_independents=0;
	bool    isautodiff       = false;
	int     aDepIndex=0;
	int     my_rank=IssmComm::GetRank();

	/*state variables: */
	IssmDouble *axp = NULL;
	IssmPDouble     *xp  = NULL;

	/*intermediary: */
	IssmPDouble *aWeightVector=NULL;
	IssmPDouble *weightVectorTimesJac=NULL;

	/*output: */
	IssmPDouble *totalgradient=NULL;

	/*AD mode on?: */
	femmodel->parameters->FindParam(&isautodiff,AutodiffIsautodiffEnum);

	if(isautodiff){

		#if defined(_HAVE_ADOLC_)
			if(VerboseAutodiff())_printf0_("   start ad core\n"); 

			/*First, stop tracing: */
			trace_off();

			if(VerboseAutodiff()){ /*{{{*/
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
		} /*}}}*/

			/*retrieve parameters: */
			femmodel->parameters->FindParam(&num_dependents,AutodiffNumDependentsEnum);
			femmodel->parameters->FindParam(&num_independents,AutodiffNumIndependentsEnum);

			/*if no dependents, no point in running a driver: */
			if(!(num_dependents*num_independents)) return;

			/*for adolc to run in parallel, we 0 out on rank~=0. But we still keep track of num_dependents:*/
			num_dependents_old=num_dependents;
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

			/* these are always needed regardless of the interpreter */
			anEDF_for_solverx_p->dp_x=xNew<double>(anEDF_for_solverx_p->max_n);
			anEDF_for_solverx_p->dp_y=xNew<double>(anEDF_for_solverx_p->max_m);

			/* Ok, now we are going to call the fos_reverse in a loop on the index, from 0 to num_dependents, so 
			 * as to generate num_dependents gradients: */

			/*Initialize outputs: */
			totalgradient=xNewZeroInit<IssmPDouble>(num_independents);

			for(aDepIndex=0;aDepIndex<num_dependents_old;aDepIndex++){

				/*initialize direction index in the weights vector: */
				aWeightVector=xNewZeroInit<IssmPDouble>(num_dependents);
				if (my_rank==0) aWeightVector[aDepIndex]=1.0;

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
				if(my_rank==0)for(i=0;i<num_independents;i++)totalgradient[i]+=weightVectorTimesJac[i];

				/*free resources :*/
				xDelete(weightVectorTimesJac);
				xDelete(aWeightVector);
			}

			/*add totalgradient to results*/
			femmodel->results->AddObject(new GenericExternalResult<IssmPDouble*>(femmodel->results->Size()+1,AutodiffJacobianEnum,totalgradient,num_independents,1,1,0.0));

			if(VerboseAutodiff())_printf0_("   end ad core\n");

			/* delete the allocated space for the parameters and Free resources:{{{*/
			xDelete(anEDF_for_solverx_p->dp_x);
			xDelete(anEDF_for_solverx_p->dp_X);
			xDelete(anEDF_for_solverx_p->dpp_X);
			xDelete(anEDF_for_solverx_p->dp_y);
			xDelete(anEDF_for_solverx_p->dp_Y);
			xDelete(anEDF_for_solverx_p->dpp_Y);
			xDelete(anEDF_for_solverx_p->dp_U);
			xDelete(anEDF_for_solverx_p->dpp_U);
			xDelete(anEDF_for_solverx_p->dp_Z);
			xDelete(anEDF_for_solverx_p->dpp_Z);
			xDelete(xp);
			xDelete(totalgradient);
			xDelete(axp); /*}}}*/

		#elif defined(_HAVE_CODIPACK_)
			fprintf(stderr, "*** Codipack adgradient_core()\n");
		#else
			_error_("Should not be requesting AD drivers when an AD library is not available!");
		#endif
	}
}

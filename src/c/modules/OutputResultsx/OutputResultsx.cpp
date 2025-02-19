/*!\file:  OutputResultsx.cpp
 * \brief: go through our finite elements, and see what results they have stored within. 
 * Then output them into serialized patch arrays, and dump to disk.
 */ 

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include <stdio.h>
#include "./OutputResultsx.h"
#include "../../shared/io/io.h"
#include "../../classes/classes.h"

void OutputResultsx(FemModel* femmodel){

	int         my_rank;
	FILE       *fid                     = NULL;
	char       *outputfilename          = NULL;
	char        outputfilename2[1000];        //easier to convert an integer with sprintf
	bool        io_gather;
	int         solutiontype;
	char*       solutiontypestring      = NULL;

	/*recover my_rank:*/
	my_rank=IssmComm::GetRank();

	/*If we are running dakota, do we want to output?*/
	bool dakota_analysis;
	femmodel->parameters->FindParam(&dakota_analysis,QmuIsdakotaEnum);
	if(dakota_analysis){
		bool dakota_output;
		femmodel->parameters->FindParam(&dakota_output,QmuOutputEnum);
		if(!dakota_output) return; 
	}

	/*Results do not include the type of solution being run	. In parallel, we output results to a filename, 
	 *therefore, we need to include the solutiontype into the filename: */
	if(my_rank==0){
		femmodel->parameters->FindParam(&solutiontype,SolutionTypeEnum);
		EnumToStringx(&solutiontypestring,solutiontype);
		femmodel->results->AddResult(new GenericExternalResult<char*>(femmodel->results->Size()+1,SolutionTypeEnum,solutiontypestring));
		xDelete<char>(solutiontypestring);
	}

#ifdef _HAVE_JAVASCRIPT_
	femmodel->parameters->FindParam(&fid,OutputFilePointerEnum);
#else

	/*Now, open file for writing*/
	_assert_(!femmodel->parameters->Exist(OutputFilePointerEnum));
	femmodel->parameters->FindParam(&outputfilename,OutputFileNameEnum);
	femmodel->parameters->FindParam(&io_gather,SettingsIoGatherEnum);

	if(io_gather){
		/*Just open the file for output on cpu 0. We are gathering the data on cpu 0 from all other cpus: */
		if(!dakota_analysis){
			if(my_rank==0) fid=pfopen0(outputfilename ,"ab+");
		}
		else{
			if(my_rank==0){
				/*a little bit complicated. Either statistic computations are requested, which means we 
				 * put our outbin files in subidirectories with numbers, or we don't, and we dump our 
				 * outbins directly in the current directory:*/
				int currEvalId ;
				int nfilesperdirectory;
				bool statistics=false;
				char* root=NULL;
				char* modelname=NULL;

				femmodel->parameters->FindParam(&currEvalId,QmuCurrEvalIdEnum);
				femmodel->parameters->FindParam(&statistics,QmuStatisticsEnum);

				if(statistics){
					femmodel->parameters->FindParam(&nfilesperdirectory,QmuNfilesPerDirectoryEnum);
					femmodel->parameters->FindParam(&root,RootPathEnum);
					femmodel->parameters->FindParam(&modelname,ModelnameEnum);
					snprintf(outputfilename2, sizeof(outputfilename2),"%s/%i/%s.outbin.%i",root,(int)(floor((currEvalId-1)/nfilesperdirectory)+1),modelname,currEvalId);
				}
				else{
					snprintf(outputfilename2, sizeof(outputfilename2),"%s.%i",outputfilename,currEvalId);
				}
				fid=pfopen0(outputfilename2,"ab+");
			}
		}
	}
	else{
		/*We are opening different  files for output on all cpus. Append the  rank to the filename, and open: */
		snprintf(outputfilename2, sizeof(outputfilename2),"%s.%i",outputfilename,my_rank);
		fid=pfopen(outputfilename2 ,"ab+");
	}

	/*Add file pointer in parameters for further calls to OutputResultsx: */
	femmodel->parameters->SetParam(fid,OutputFilePointerEnum);
#endif

	/*Write results to disk: */
	femmodel->results->Write(femmodel->parameters);

#ifdef _HAVE_JAVASCRIPT_
	/*Delete and reinitialize results, in parallel: */
	femmodel->results->clear();
#else
	femmodel->parameters->Delete(OutputFilePointerEnum);

	/*Delete and reinitialize results, in parallel: */
	femmodel->results->clear();

	/*Close output file? :*/
	if(io_gather){
		if(!dakota_analysis){
			if(my_rank==0) pfclose(fid,outputfilename);
		}
		else{
			if(my_rank==0) pfclose(fid,outputfilename2);
		}
	}
	else pfclose(fid,outputfilename2);
#endif

	/*Clean up and return*/
	xDelete<char>(outputfilename);
}

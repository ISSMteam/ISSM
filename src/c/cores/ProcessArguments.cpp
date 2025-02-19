/*!\file:  ProcessArguments.cpp
 * \brief: process arguments
 */ 

#include <stdio.h>
#include <cstring>

#include "../shared/shared.h"

void ProcessArguments(int* solution_type,char** pbinfilename,char** poutbinfilename,char** ptoolkitsfilename,char** plockfilename,char** prestartfilename, char** prootpath, char** pmodelname, int argc,char **argv){

	/*Check input arguments*/
	if(argc<2) _error_("Usage error: no solution requested");
	if(argc<3) _error_("Usage error: missing execution directory");
	if(argc<4) _error_("Usage error: missing model name");

	/*Get some arguments*/
	*solution_type = StringToEnumx(argv[1]);
	char* rootpatharg = argv[2];
	char* modelname   = xNew<char>(strlen(argv[3])+1); 
	xMemCpy<char>(modelname,argv[3],strlen(argv[3])+1);

	/*Recover myrank and length of string "my_rank" */
	int my_rank     = IssmComm::GetRank();
	int rank_length = (my_rank == 0 ? 1 : (int)(log10(static_cast<double>(my_rank))+1)); 

	/*Create rootpath from argument*/
	int   rootpath_len = strlen(rootpatharg)+2;
	char* rootpath = xNew<char>(rootpath_len);
	snprintf(rootpath, rootpath_len,"%s/",rootpatharg);

	/*Create all file paths*/
	int base_length = strlen(rootpath)+strlen(modelname);
	int binfilename_len      = base_length+strlen(".bin")     +1;
	int outbinfilename_len   = base_length+strlen(".outbin")  +1;
	int toolkitsfilename_len = base_length+strlen(".toolkits")+1;
	int lockfilename_len     = base_length+strlen(".lock")    +1;
	int restartfilename_len  = base_length+strlen("_rank")+rank_length+strlen(".rst")+1;
	char* binfilename      = xNew<char>(binfilename_len);      snprintf(binfilename,      binfilename_len,     "%s%s%s",rootpath,modelname,".bin");
	char* outbinfilename   = xNew<char>(outbinfilename_len);   snprintf(outbinfilename,   outbinfilename_len,  "%s%s%s",rootpath,modelname,".outbin");
	char* toolkitsfilename = xNew<char>(toolkitsfilename_len); snprintf(toolkitsfilename, toolkitsfilename_len,"%s%s%s",rootpath,modelname,".toolkits");
	char* lockfilename     = xNew<char>(lockfilename_len);     snprintf(lockfilename,     lockfilename_len,    "%s%s%s",rootpath,modelname,".lock");
	char* restartfilename  = xNew<char>(restartfilename_len);  snprintf(restartfilename,  restartfilename_len, "%s%s%s%i%s",rootpath,modelname,"_rank",my_rank,".rst");

	/*Clean up and assign output pointer*/
	*pbinfilename      = binfilename;
	*poutbinfilename   = outbinfilename;
	*ptoolkitsfilename = toolkitsfilename;
	*plockfilename     = lockfilename;
	*prestartfilename  = restartfilename;
	*prootpath         = rootpath;
	*pmodelname        = modelname;
}

/*!\file:  pfopen.cpp
 * \brief fopen wrapper for parallel solution
 */ 

#ifdef HAVE_CONFIG_H
#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include <stdio.h>
#include "../Print/Print.h"
#include "../Comm/IssmComm.h"
#include "../../Exceptions/exceptions.h"

FILE* pfopen0(char* filename,const char* format){

	FILE* fid=NULL;

	/*recover my_rank:*/
	int my_rank  = IssmComm::GetRank();
	if(my_rank) _error_("This function should only be called by cpu 0");

	/*Open handle to data on disk*/
	fid = fopen(filename,format);
	if(fid==NULL)_error_("could not open file " << filename << " for binary reading or writing");

	return fid;
}
FILE* pfopen(char* filename,const char* format,bool errorout){

	FILE* fid=NULL;

	/*recover my_rank:*/
	int my_rank  = IssmComm::GetRank();
	int num_proc = IssmComm::GetSize();

	/*Open handle to data on disk (one by one to avoid errors)*/
	for(int i=0;i<num_proc;i++){
		if(my_rank==i) fid = fopen(filename,format);
		ISSM_MPI_Barrier(IssmComm::GetComm());
	}
	if(errorout && fid==NULL)_error_("could not open file " << filename << " for binary reading or writing");

	return fid;
}

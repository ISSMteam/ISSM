/*!\file:  WriteLockFile.cpp
 * \brief
 */ 
#include <cstdio>
#include "../../Exceptions/exceptions.h"
#include "../Comm/IssmComm.h"
#include "../Print/Print.h"
#include <cstdio>

void WriteLockFile(char* filename){

	/*recover my_rank:*/
	int my_rank=IssmComm::GetRank();

	/* output: */
	FILE* fid=NULL;

	/* Open lock file and write 1 into it: */
	if(my_rank==0){
		fid=fopen(filename,"w");
		if(fid==NULL) _error_("error message: could not open lock file " << filename);

		/*Close file: */
		if(fclose(fid)!=0) _error_("could not close lock file " << filename);
	}

}	

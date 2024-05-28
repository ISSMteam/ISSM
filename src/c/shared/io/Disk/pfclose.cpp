/*!\file:  pfclose.cpp
 * \brief fclose wrapper for parallel solution
 */ 

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include <stdio.h>
#include "../../shared.h"

void pfclose(FILE* fid,char* filename){

	/*Close file handle: */
	_assert_(fid);
	if(fclose(fid)!=0)_error_("could not close file " << filename);
}

/*!\file OptionUtilities.cpp
 * \brief: implementation of the options utilities
 */

/*Headers:*/
/*{{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "../classes.h"
#include "../../shared/shared.h"
/*}}}*/

int StringFromSize(char* cstr, int* size, int ndims) {/*{{{*/

	sprintf(&cstr[0],"[");
	for(int i=0; i<ndims-1; i++) sprintf(&cstr[strlen(cstr)],"%dx",size[i]);
	sprintf(&cstr[strlen(cstr)],"%d]",size[ndims-1]);

	return(0);
}/*}}}*/

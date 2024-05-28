/* \file ApiPrintf.c:
 * \brief: API specific symbols from libISSMCore that we need to resolve here
 */

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./javascriptio.h"

/*Javascript printf i/o: */
void ApiPrintf(const char* string){

	/*use mexPrintf in matlab: */
	printf(string);
	return;
}

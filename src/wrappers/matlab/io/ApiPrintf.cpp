/* \file ApiPrintf.c:
 * \brief: API specific symbols from libISSMCore that we need to resolve here
 */

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./matlabio.h"

/*Matlab printf i/o: */
void ApiPrintf(const char* string){

	/*use mexPrintf in matlab: */
	mexPrintf(string); /*Messes up percentages (like "interpolation progress: 100.00%")*/
	//printf("%s",string);
	return;
}

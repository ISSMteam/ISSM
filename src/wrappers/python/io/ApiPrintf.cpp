/* \file ApiPrintf.c:
 * \brief: pyton api specific symbols which are unresolved from libISSMCore.a
 */

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./pythonio.h"
#include "../../c/shared/shared.h"

/*Python printf i/o: */
void ApiPrintf(const char* string){

	/*use printf: */
	printf("%s",string);
	return;
}

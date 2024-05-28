/*!\file: issm.h
 * \brief prototype wrapper for issm.h
 */ 

#ifndef _ISSM_H_
#define _ISSM_H_

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./globals.h" //only include this header file once!
#include "../shared/shared.h"
#include "../classes/classes.h"
#include "../toolkits/toolkits.h"
#include "../cores/cores.h"
#include "../solutionsequences/solutionsequences.h"
#include "../modules/modules.h"

/*Environment*/
ISSM_MPI_Comm EnvironmentInit(int argc,char** argv);
void EnvironmentFinalize(void);

#endif //ifndef _ISSM_H_

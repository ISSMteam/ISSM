/*!\file:  IssmToolkitUtils.h
 * \brief routines used throughout the ISSM toolkit
 */ 

#ifndef _ISSM_TOOLKIT_UTILS_H_
#define _ISSM_TOOLKIT_UTILS_H_

/*Headers:*/
/*{{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

/*}}}*/

int IssmMatTypeFromToolkitOptions(void);
int IssmVecTypeFromToolkitOptions(void);
int IssmSolverTypeFromToolkitOptions(void);

#endif //#ifndef _ISSMMAT_H_

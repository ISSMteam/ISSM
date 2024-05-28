/*!\file:  OutputResultsx.h
 * \brief header file for outputing results
 */ 

#ifndef _OUTPUTRESULTSX_H
#define _OUTPUTRESULTSX_H

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "../../classes/classes.h"

void OutputResultsx(FemModel* femmodel);

#endif  /* _OUTPUTRESULTS_H */

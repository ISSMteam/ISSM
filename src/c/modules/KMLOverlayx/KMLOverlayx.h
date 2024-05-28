/*!\file:  KMLOverlayx.h
 * \brief header file for kml file overlay routines.
 */ 

#ifndef _KMLOVERLAYX_H
#define _KMLOVERLAYX_H

#include <float.h>    /*  DBL_MAX  */
#include "../../classes/classes.h"

/* local prototypes: */
void KMLOverlayx(int* ierror,
				 double* lataxis, double* longaxis,
				 int nimages, char** pimages,
				 FILE* fid);

#endif  /* _KMLOVERLAYX_H */

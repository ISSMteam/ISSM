/*!\file:  Kml2Expx.h
 * \brief header file for kml to exp conversion routines.
 */ 

#ifndef _KML2EXPX_H
#define _KML2EXPX_H

#include <float.h>    /*  DBL_MAX  */
#include "../../classes/classes.h"

/* local prototypes: */
int Kml2Expx(char* filkml,char* filexp,
			 int sgn);
int Kml2Expx(char* filkml,char* filexp,
			 int sgn,double cm,double sp);

#endif  /* _KML2EXPX_H */

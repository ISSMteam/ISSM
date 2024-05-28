/*!\file:  Shp2Kmlx.h
 * \brief header file for shp to kml conversion routines.
 */ 

#ifndef _SHP2KMLX_H
#define _SHP2KMLX_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#ifdef _HAVE_SHAPELIB_ //only works if Shapelib library has been compiled in.

#include "shapefil.h"

#endif

#include "../../classes/classes.h"

/* local prototypes: */
int Shp2Kmlx(char* filshp,char* filkml, int sgn);
int Shp2Kmlx(char* filshp,char* filkml, int sgn,double cm,double sp);

#endif  /* _SHP2KMLX_H */

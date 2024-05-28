/*!\file:  ElementCoordinatesx.h
 */ 

#ifndef _ELEMENT_COORDINATESX_H
#define _ELEMENT_COORDINATESX_H

#include "../../classes/classes.h"

/* local prototypes: */
void ElementCoordinatesx( IssmDouble** pxe, IssmDouble** pye, IssmDouble** pze,IssmDouble** pareae, Elements* elements,bool spherical=false);
void ElementCoordinatesx( IssmDouble** plonge, IssmDouble** plate, IssmDouble** pareae, Elements* elements);

#endif  /* _ELEMENT_COORDINATESX_H */

/*!\file:  SurfaceLogVxVyMisfitx.h
 * \brief header file for inverse methods misfit computation
 */ 

#ifndef _SURFACELOGVXVYMISFITX_H
#define _SURFACELOGVXVYMISFITX_H

#include "../../classes/classes.h"

/* local prototypes: */
void SurfaceLogVxVyMisfitx( IssmDouble* pJ, Elements* elements,Nodes* nodes, Vertices* vertices, Loads* loads, Materials* materials, Parameters* parameters);
IssmDouble SurfaceLogVxVyMisfit(Element* element);

#endif

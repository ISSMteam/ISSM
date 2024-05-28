/*!\file:  SurfaceLogVelMisfitx.h
 * \brief header file for inverse methods misfit computation
 */ 

#ifndef _SURFACELOGVELMISFITX_H
#define _SURFACELOGVELMISFITX_H

#include "../../classes/classes.h"

/* local prototypes: */
void SurfaceLogVelMisfitx( IssmDouble* pJ, Elements* elements,Nodes* nodes, Vertices* vertices, Loads* loads, Materials* materials, Parameters* parameters);
IssmDouble SurfaceLogVelMisfit(Element* element);

#endif

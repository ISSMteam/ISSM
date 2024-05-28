/*!\file:  SurfaceRelVelMisfitx.h
 * \brief header file for inverse methods misfit computation
 */ 

#ifndef _SURFACERELVELMISFITX_H
#define _SURFACERELVELMISFITX_H

#include "../../classes/classes.h"

/* local prototypes: */
void SurfaceRelVelMisfitx( IssmDouble* pJ, Elements* elements,Nodes* nodes, Vertices* vertices, Loads* loads, Materials* materials, Parameters* parameters);
IssmDouble SurfaceRelVelMisfit(Element* element);

#endif

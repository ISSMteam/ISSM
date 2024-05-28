/*!\file:  SurfaceAbsVelMisfitx.h
 * \brief header file for inverse methods misfit computation
 */ 

#ifndef _SURFACEABSVELMISFITX_H
#define _SURFACEABSVELMISFITX_H

#include "../../classes/classes.h"

/* local prototypes: */
void SurfaceAbsVelMisfitx( IssmDouble* pJ, Elements* elements,Nodes* nodes, Vertices* vertices, Loads* loads, Materials* materials, Parameters* parameters);
IssmDouble SurfaceAbsVelMisfit(Element* element);

#endif

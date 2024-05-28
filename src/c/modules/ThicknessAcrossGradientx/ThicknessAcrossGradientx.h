/*!\file:  ThicknessAcrossGradientx.h
 * \brief header file for inverse methods misfit computation
 */ 

#ifndef _THICKNESSACROSSGRADIENT_H
#define _THICKNESSACROSSGRADIENT_H

#include "../../classes/classes.h"

/* local prototypes: */
void ThicknessAcrossGradientx( IssmDouble* pJ, Elements* elements,Nodes* nodes, Vertices* vertices, Loads* loads, Materials* materials, Parameters* parameters);
IssmDouble ThicknessAcrossGradient(Element* element);

#endif

/*!\file:  ThicknessAlongGradientx.h
 * \brief header file for inverse methods misfit computation
 */ 

#ifndef _THICKNESSALONGGRADIENT_H
#define _THICKNESSALONGGRADIENT_H

#include "../../classes/classes.h"

/* local prototypes: */
void ThicknessAlongGradientx( IssmDouble* pJ, Elements* elements,Nodes* nodes, Vertices* vertices, Loads* loads, Materials* materials, Parameters* parameters);
IssmDouble ThicknessAlongGradient(Element* element);

#endif

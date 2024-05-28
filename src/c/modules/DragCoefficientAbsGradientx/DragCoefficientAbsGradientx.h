/*!\file:  DragCoefficientAbsGradientx.h
 * \brief header file for inverse methods misfit computation
 */ 

#ifndef _DRAGCOEFFABSGRADX_H
#define _DRAGCOEFFABSGRADX_H

#include "../../classes/classes.h"

/* local prototypes: */
void DragCoefficientAbsGradientx( IssmDouble* pJ, Elements* elements,Nodes* nodes, Vertices* vertices, Loads* loads, Materials* materials, Parameters* parameters);
IssmDouble DragCoefficientAbsGradient(Element* element);

#endif

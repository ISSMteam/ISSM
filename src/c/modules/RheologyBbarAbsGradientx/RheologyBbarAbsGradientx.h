/*!\file:  RheologyBbarAbsGradientx.h
 * \brief header file for inverse methods misfit computation
 */ 

#ifndef _RHEOLOGYBBARGRADIENTX_H
#define _RHEOLOGYBBARGRADIENTX_H

#include "../../classes/classes.h"

/* local prototypes: */
void RheologyBbarAbsGradientx( IssmDouble* pJ, Elements* elements,Nodes* nodes, Vertices* vertices, Loads* loads, Materials* materials, Parameters* parameters);
IssmDouble RheologyBbarAbsGradient(Element* element);

#endif

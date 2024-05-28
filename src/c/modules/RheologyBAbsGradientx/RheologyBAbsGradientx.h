/*!\file:  RheologyBAbsGradientx.h
 * \brief header file for inverse methods misfit computation
 */ 

#ifndef _RHEOLOGYBGRADIENTX_H
#define _RHEOLOGYBGRADIENTX_H

#include "../../classes/classes.h"

/* local prototypes: */
void RheologyBAbsGradientx( IssmDouble* pJ, Elements* elements,Nodes* nodes, Vertices* vertices, Loads* loads, Materials* materials, Parameters* parameters);
IssmDouble RheologyBAbsGradient(Element* element);

void RheologyBInitialguessMisfitx( IssmDouble* pJ, Elements* elements,Nodes* nodes, Vertices* vertices, Loads* loads, Materials* materials, Parameters* parameters);
IssmDouble RheologyBInitialguessMisfit(Element* element);

#endif

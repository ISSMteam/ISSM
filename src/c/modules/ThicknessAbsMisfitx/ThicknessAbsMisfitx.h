/*!\file:  ThicknessAbsMisfitx.h
 * \brief header file for inverse methods misfit computation
 */ 

#ifndef _THICKNESSABSMISFITX_H
#define _THICKNESSABSMISFITX_H

#include "../../classes/classes.h"

/* local prototypes: */
void ThicknessAbsMisfitx( IssmDouble* pJ, Elements* elements,Nodes* nodes, Vertices* vertices, Loads* loads, Materials* materials, Parameters* parameters);
IssmDouble ThicknessAbsMisfit(Element* element);

#endif

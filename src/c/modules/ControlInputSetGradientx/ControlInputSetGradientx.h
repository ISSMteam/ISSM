/*!\file:  ControlInputSetGradientx.h
 */ 

#ifndef _CONTROLINPUTSSETGRADIENTX_H
#define _CONTROLINPUTSSETGRADIENTX_H

#include "../../classes/classes.h"

void	ControlInputSetGradientx(Elements* elements,Nodes* nodes, Vertices* vertices,Loads* loads, Materials* materials,  Parameters* parameters,IssmDouble* gradient);
void	ControlInputSetGradientx(Elements* elements,Nodes* nodes, Vertices* vertices,Loads* loads, Materials* materials,  Parameters* parameters,Vector<IssmDouble>* gradient);

#endif

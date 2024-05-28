/*!\file:  Gradjx.h
 * \brief header file for inverse methods gradient computation
 */ 

#ifndef _GRADJX_H
#define _GRADJX_H

#include "../../classes/classes.h"
#include "../../analyses/analyses.h"

/* local prototypes: */
void Gradjx(Vector<IssmDouble>** pgrad_g,IssmDouble** pgrad_norm,Elements* elements,Nodes* nodes, Vertices* vertices, Loads* loads, Materials* materials,Parameters* parameters);
void Gradjx(IssmDouble** pgrad_g,IssmDouble** pgrad_norm,Elements* elements,Nodes* nodes, Vertices* vertices, Loads* loads, Materials* materials,Parameters* parameters);

#endif  /* _GRADJX_H */

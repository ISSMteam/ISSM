/*!\file:  NodalValuex.h
 * \brief header file for NodalValuex
 */ 

#ifndef _NODALVALUEX_H
#define _NODALVALUEX_H

#include "../../classes/classes.h"

/* local prototypes: */
void NodalValuex( IssmDouble* pnodalvalue, int natureofdataenum,Elements* elements,Nodes* nodes, Vertices* vertices, Loads* loads, Materials* materials, Parameters* parameters);

#endif  /* _NODALVALUEX_H */

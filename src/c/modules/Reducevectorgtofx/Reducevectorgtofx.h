/*!\file:  Reducevectorgtofx.h
 * \brief reduce petsc vector from g set to f set (free dofs), using the nodeset partitioning 
 * vectors.
 */ 

#ifndef _REDUCEVECTORGTOFX_H
#define _REDUCEVECTORGTOFX_H

#include "../../classes/classes.h"

/* local prototypes: */
void Reducevectorgtofx(Vector<IssmDouble>** puf, Vector<IssmDouble>* ug, Nodes* nodes,Parameters* parameters);

#endif  /* _REDUCEVECTORGTOFX_H */
